library(tximport)
library(edgeR)
library(iCOBRA)
library(dplyr)
library(Hmisc)
library(BiocParallel)
library(DESeq2)
library(DEXSeq)
library(biomaRt)
library(Biostrings)
library(scales)
source("/Users/koenvandenberge/Dropbox/PhD/Research/transcriptLevel/helpFunctionsCharlotte.R")
muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
           "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
           "#90C987","#CAEDAB","#777777")
meta <- data.frame(sample = paste0("sample", c("A1", "A2", "A3", "B1", "B2", "B3")),
                   condition = c("A", "A", "A", "B", "B", "B"),
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$sample


### get GTF file
library(rtracklayer)
gtfWholeGenome="/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.gtf"
gtff <- import(gtfWholeGenome)
gtff2 <- subset(gtff, seqnames == "1")
export(gtff2, "/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.chr1.gtf", format = "gtf")

## generate tx 2 gene map
cdna_fasta <-"/Users/koenvandenberge/Dropbox/edgeR_zeroinflation/simulated_data_bulk_tx/sim2_human/reference_files/Homo_sapiens.GRCh37.71.cdna.chr1.fa"
gtf <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.chr1.gtf"
feature_lengths_file <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/feature_lengths.Rdata"
tx_gene_file <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/tx_gene_map.Rdata"
calc_lengths_mapping(gtf = gtf, cdna_fasta = cdna_fasta, 
                       feature_lengths_file = feature_lengths_file, 
                       tx_gene_file = tx_gene_file) 
load(feature_lengths_file)
load(tx_gene_file) #gene2tx
tx2gene = gene2tx[,c(2,1)]

truth_tx_file <- "/Users/koenvandenberge/Dropbox/edgeR_zeroinflation/simulated_data_bulk_tx/sim2_human/truth_transcript.txt"
truth_tx <- read.delim(truth_tx_file, header = TRUE, as.is = TRUE)


### derive salmon transcript counts
salmon_basedir <- "/Users/koenvandenberge/Dropbox/edgeR_zeroinflation/simulated_data_bulk_tx/sim2_human/salmon"

salmon_files <- list.files(salmon_basedir, pattern = "sample", full.names = TRUE)
salmon_files <- salmon_files[file.info(salmon_files)$isdir]
salmon_files <- paste0(salmon_files, "/quant.sf")
salmon_files <- salmon_files[file.exists(salmon_files)]
names(salmon_files) <- basename(gsub("/quant.sf", "", salmon_files))
txi_salmonsimplesum <- tximport(files = salmon_files, type = "salmon", txIn = TRUE,
                                txOut = FALSE, countsFromAbundance = "no", 
                                tx2gene = tx2gene)

txi_salmonscaledtpm <- tximport(files = salmon_files, type = "salmon", txIn = TRUE,
                                txOut = FALSE, countsFromAbundance = "scaledTPM", 
                                 tx2gene = tx2gene)

txi_salmontx <- tximport(files = salmon_files, type = "salmon", txIn = TRUE,
                         txOut = TRUE, countsFromAbundance = "no", tx2gene = tx2gene)

salmon_quant <- list(geneCOUNT_sal_simplesum = txi_salmonsimplesum$counts,
                     geneCOUNT_sal_scaledTPM = txi_salmonscaledtpm$counts,
                     avetxlength = txi_salmonsimplesum$length,
                     geneTPM_sal = txi_salmonsimplesum$abundance,
                     txTPM_sal = txi_salmontx$abundance,
                     txCOUNT_sal = txi_salmontx$counts,
                     txi_salmonsimplesum = txi_salmonsimplesum,
                     txi_salmonscaledtpm = txi_salmonscaledtpm,
                     txi_salmontx = txi_salmontx)
stopifnot(all(colnames(salmon_quant$txCOUNT_sal) == rownames(meta)))

#### DEXSEQ on salmon transcript counts
## discard genes with only one transcript
genesWithOneTx <- names(table(tx2gene$gene))[table(tx2gene$gene)==1]
txFromGenesWithOneTx <- tx2gene$transcript[match(genesWithOneTx,tx2gene$gene)]
txCount <- round(salmon_quant$txCOUNT_sal)
#avoid NA p-values
txCount <- txCount[!rownames(txCount)%in%txFromGenesWithOneTx,]
txCount <- txCount[!rowSums(txCount)==0,]

## regular DEXSeq analysis
dxd_sal <- DEXSeqDataSet(countData = txCount, 
                         sampleData = meta, 
                         design = ~sample + exon + condition:exon,
                         featureID = rownames(txCount),
                         groupID = tx2gene$gene[match(rownames(txCount),
                                                      tx2gene$transcript)])
dxd_sal <- estimateSizeFactors(dxd_sal)
dxd_sal <- estimateDispersions(dxd_sal)
dxd_sal <- testForDEU(dxd_sal)
dxr_sal <- DEXSeqResults(dxd_sal) #still couple NA pvals but not a lot
hist(dxr_sal$pvalue)
qval_dtu_salmon <- perGeneQValue(dxr_sal)

## stage-wise DEXSeq analysis
significantGenes <- names(qval_dtu_salmon)[which(qval_dtu_salmon<.05)]
alphaAdjusted <- 0.05*length(significantGenes)/length(qval_dtu_salmon)
genesStageII <- dxr_sal$groupID[dxr_sal$groupID%in%significantGenes]
uniqueGenesStageII <- unique(genesStageII)
txStageII <- dxr_sal$featureID[dxr_sal$groupID%in%significantGenes]
pvalStageII <- dxr_sal$pvalue[dxr_sal$featureID%in%txStageII]
pvalGeneList <- list()
for(i in 1:length(uniqueGenesStageII)){
    id <- which(genesStageII==uniqueGenesStageII[i])
    pvalHlp <- pvalStageII[id]
    names(pvalHlp) <- txStageII[id]
    pvalGeneList[[i]] <- pvalHlp
}
padjGeneList <- lapply(pvalGeneList,function(x) p.adjust(x,method="holm"))


## ROC gene-level analysis
truth_gene_file <- "/Users/koenvandenberge/Dropbox/edgeR_zeroinflation/simulated_data_bulk_tx/sim2_human/truth_gene.txt"
truth_gene <- read.delim(truth_gene_file, header = TRUE, as.is = TRUE, row.names = 1)
cobra <- COBRAData(padj = data.frame(salmon_dexseq = qval_dtu_salmon,
                                     row.names = names(qval_dtu_salmon),
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth_gene, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "diffisouse")
cobraplotGene <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green","steelblue"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplotGene)
plot_roc(cobraplotGene,xaxisrange=c(0,0.1))

## ROC transcript-level analysis
rownames(truth_tx) <- truth_tx$transcript
cobra <- COBRAData(padj = data.frame(salmon_dexseq = p.adjust(dxr_sal$pvalue,"BH"),
                                     row.names = dxr_sal$featureID,
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth_tx, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "status")
cobraplotTx <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green","steelblue"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplotTx)
plot_roc(cobraplotTx,xaxisrange=c(0,0.1))


### combine in one plot
library(scales)
plot(x=cobraplotGene@fdrtprcurve$FDR,y=cobraplotGene@fdrtprcurve$TPR, type="l",col=2,lwd=1, xlim=c(0,0.6), ylab="True Positive Rate", xlab="False Discovery Rate")
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col=2, pch="o", cex=1.2)
lines(x=cobraplotTx@fdrtprcurve$FDR,y=cobraplotTx@fdrtprcurve$TPR,col=4,lwd=1)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col=4, pch="o", cex=1.2)
abline(v=c(0.01,0.05,0.1,seq(0.1,1,.1)), col=alpha("grey",.8), lty=2)
legend("topleft",c("gene-level","transcript-level"),col=c(2,4),lty=1, bty="n")








########
## ROC transcript-level analysis through stage-wise testing
## since the number of p-values changes for a different alpha through stage-wise testing and the stage II has an adjusted alpha, I will write a function
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
pvalDexSeq <- dxr_sal$pvalue
names(pvalDexSeq) <- dxr_sal$featureID
pvalDexSeq[is.na(pvalDexSeq)] <- 1
padjDexSeq <- p.adjust(pvalDexSeq,method="BH")
tprSW <- vector(length=length(pvalSeq))
fprSW <- vector(length=length(pvalSeq))
fdrSW <- vector(length=length(pvalSeq))
for(id in 1:length(pvalSeq)){
    print(id)
    ### DEXSeq stage-wise
    significantGenes <- names(qval_dtu_salmon)[which(qval_dtu_salmon<=pvalSeq[id])]
    if(length(significantGenes)==0){
	tprSW[id]=0
	fprSW[id]=0
	fdrSW[id]=0
	next
    } else {
    alphaAdjusted <- pvalSeq[id]*length(significantGenes)/length(qval_dtu_salmon)
    genesStageII <- dxr_sal$groupID[dxr_sal$groupID%in%significantGenes]
    uniqueGenesStageII <- unique(genesStageII)
    txStageII <- dxr_sal$featureID[dxr_sal$groupID%in%significantGenes]
    pvalStageII <- dxr_sal$pvalue[dxr_sal$featureID%in%txStageII]
    pvalGeneList <- list()
    for(i in 1:length(uniqueGenesStageII)){
        idHlp <- which(genesStageII==uniqueGenesStageII[i])
        pvalHlp <- pvalStageII[idHlp]
        names(pvalHlp) <- txStageII[idHlp]
        pvalGeneList[[i]] <- pvalHlp
    }
    padjGeneList <- lapply(pvalGeneList,function(x) p.adjust(x,method="holm"))
    padj <- unlist(padjGeneList)
    positives <- names(padj[padj<=alphaAdjusted])
    tprSW[id] <- mean(truth_tx$transcript[truth_tx$status==1]%in%positives)
    truePosId <- truth_tx$status[match(positives,truth_tx$transcript)]
    fprSW[id] <- mean(truth_tx$transcript[truth_tx$status==0]%in%positives)
    fdrSW[id] <- 1-mean(truePosId)
    }
}


### DEXSeq regular
evalDexSeqRegular <- t(sapply(pvalSeq,function(alpha){
	   positives <- names(padjDexSeq[which(padjDexSeq<alpha)])
	   tprRegular <- mean(truth_tx$transcript[truth_tx$status==1]%in%positives)
	   fprRegular <- mean(truth_tx$transcript[truth_tx$status==0]%in%positives)
	   truePosId <- truth_tx$status[match(positives,truth_tx$transcript)]
	   fdrRegular <- 1-mean(truePosId)
	   return(cbind(tprRegular=tprRegular, fprRegular=fprRegular, fdrRegular=fdrRegular))
				   }))
colnames(evalDexSeqRegular) <- c("tpr","fpr","fdr")

## ROC curve
plot(x=fprSW,y=tprSW,type="n", xlim=c(0,0.15))
lines(x=fprSW,y=tprSW,col=2,lwd=2)
points(x=fprSW[516],y=tprSW[516],pch=3,col=2)
lines(x=evalDexSeqRegular[,"fpr"],y=evalDexSeqRegular[,"tpr"],col=3,lwd=2)
points(x=evalDexSeqRegular[516,"fpr"],y=evalDexSeqRegular[516,"tpr"],pch=3,col=3)

## FDR-TPR curve
plot(x=fdrSW,y=tprSW,type="n")
lines(x=fdrSW,y=tprSW,col=2,lwd=2)
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],pch="o",col=2)
lines(x=evalDexSeqRegular[,"fdr"],y=evalDexSeqRegular[,"tpr"],col=3,lwd=2)
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],col=3,pch="o")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)

## there is no big difference in this simulation study, because there is no huge multiple testing correction: the number of genes is only 1800.







