### simulation Drosophila from Genome Biology paper Soneson 2016
baseDir <- "/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/diff_splicing_comparison_drosophila/"
files=list.files(baseDir)
fastaFiles <- files[grep(x=files,pattern=".fq.gz")]
meta=read.delim("/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/E-MTAB-3766.sdrf.txt",header=TRUE) #also includes fragment length and SD for kallisto
names(fastaFiles) <- (meta$Array.Data.File[1:12])
kallistoIndex="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/KallistoIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel"
txConversionFile="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/KallistoIndex/TranscriptID_conversion.txt"
kallistoDir="/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/diff_splicing_comparison_drosophila/quantifications/kallisto/"
truthFile="/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/Dm_truth.txt"
library(DEXSeq)
library(iCOBRA)
library(scales)
library(dplyr)

### kallisto quantification
sample <- unlist(lapply(strsplit(names(fastaFiles),split="_"),function(x) x[3]))
fileNames <- unlist(lapply(strsplit(names(fastaFiles),split=".",fixed=TRUE),function(x) x[1]))
fileNames <- substr(x=fileNames,1,11)[seq(1,12,2)]
#for(i in 1:6){
#    pairedFasta <- fastaFiles[sample==i]
#    cmd <- paste0("kallisto quant -i ",kallistoIndex,
#	      " -o ",kallistoDir,fileNames[i],
#	      " -b 30 ",
#	      baseDir,pairedFasta[1]," ",baseDir,pairedFasta[2])
#    message(cmd)
#    system(cmd)
#}

### get kallisto results
files2 <- list.files(kallistoDir)
sampleDirs <- files2[grep(x=files2,pattern="Dm_sample_[1-6]")]
dir=sampleDirs[1]
hlp=read.table(paste0(kallistoDir,"/",dir,"/abundance.tsv"), header=TRUE)
data <- as.data.frame(sapply(sampleDirs,function(dir) read.table(paste0(kallistoDir,"/",dir,"/abundance.tsv"), header=TRUE)[,"est_counts"]), row.names=hlp$target_id)
kal2tx=read.table(txConversionFile)
colnames(kal2tx) <- c("kallisto","transcript")
#rownames(data) <- kal2tx$transcript[match(kal2tx$kallisto,rownames(data))]
rownames(data) <- kal2tx$transcript[match(kal2tx$kallisto,rownames(data))]
#data <- data[!rowSums(data)==0,]

## the truth_tx file is incorrect for tx level evaluation: all tx from a gene get a differential splicing status, but actually only two of them should have.
truth <- read.table(truthFile,header=TRUE)
truth_gene <- truth[,c(1,3)]
truth_gene <- truth_gene[!duplicated(truth$gene_id),]
rownames(truth_gene) <- truth_gene[,1]

truth_tx <- truth[,2:3]
truth_tx$transcript_id <- as.character(truth_tx$transcript_id)
rownames(truth_tx) <- truth_tx[,1]

tx2gene <- truth[,1:2]
tx2gene$gene_id <- as.character(tx2gene$gene_id)
tx2gene$transcript_id <- as.character(tx2gene$transcript_id)

simFilesDir <- "/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/"
simFiles=list.files(simFilesDir)
firstSampleFiles <- simFiles[grep(x=simFiles,pattern="Dm_sample[1-3]")]
sample1Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[1]),header=TRUE)
sample2Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[2]),header=TRUE)
sample3Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[3]),header=TRUE)

genesWithOneTx <- names(table(sample1Tpm$gene_id)[table(sample1Tpm$gene_id)==1])
sample1Tpm <- sample1Tpm[!sample1Tpm$gene_id%in%genesWithOneTx,]
hlp=group_by(sample1Tpm,by=gene_id)
highestExpressedTxPerGene <- dplyr::summarize(hlp,max1=transcript_id[order(TPM,decreasing=TRUE)[1]],max2=transcript_id[order(TPM,decreasing=TRUE)[2]])
head(highestExpressedTxPerGene)

#consistent across samples?
sample2Tpm <- sample2Tpm[!sample2Tpm$gene_id%in%genesWithOneTx,]
hlp=group_by(sample2Tpm,by=gene_id)
highestExpressedTxPerGene2 <- dplyr::summarize(hlp,max1=transcript_id[order(TPM,decreasing=TRUE)[1]],max2=transcript_id[order(TPM,decreasing=TRUE)[2]])
head(highestExpressedTxPerGene2)
#yes, consistent across samples

truthTx <- data.frame(transcript=truth_tx[,1],tx_ds_status=0,row.names=rownames(truth_tx))
simulatedDTUGenes <- rownames(truth_gene[truth_gene$gene_ds_status==1,])
length(simulatedDTUGenes) #should be 1000
highestExpressedTxForSignificantGenes <- highestExpressedTxPerGene[highestExpressedTxPerGene$by%in%simulatedDTUGenes,]
dtuTx <- as.character(unlist(c(as.data.frame(highestExpressedTxForSignificantGenes)[,2:3])))
truthTx[dtuTx,"tx_ds_status"]=1


### DEXSeq analysis
genesWithOneTx <- names(table(tx2gene$gene))[table(tx2gene$gene)==1]
txFromGenesWithOneTx <- tx2gene$transcript[match(genesWithOneTx,tx2gene$gene)]
txCount <- round(data)
#avoid NA p-values
txCount <- txCount[!rownames(txCount)%in%txFromGenesWithOneTx,]
txCount <- txCount[!rowSums(txCount)==0,]

geneTx <- tx2gene$gene_id[match(rownames(txCount),tx2gene$transcript_id)]
sampleData <- data.frame(condition=factor(rep(c("A","B"),each=3)))
dxd <- DEXSeqDataSet(countData = txCount, 
                         sampleData = sampleData, 
                         design = ~ sample + exon + condition:exon,
                         featureID = rownames(txCount),
                         groupID = geneTx)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxr <- DEXSeqResults(dxd)
qvalDxr <- perGeneQValue(dxr)

## stage-wise DEXSeq analysis
significantGenes <- names(qvalDxr)[which(qvalDxr<=.05)]
alphaAdjusted <- 0.05*length(significantGenes)/length(qvalDxr)
genesStageII <- dxr$groupID[dxr$groupID%in%significantGenes]
uniqueGenesStageII <- unique(genesStageII)
txStageII <- dxr$featureID[dxr$groupID%in%significantGenes]
pvalStageII <- dxr$pvalue[dxr$featureID%in%txStageII]
pvalGeneList <- list()
for(i in 1:length(uniqueGenesStageII)){
    id <- which(genesStageII==uniqueGenesStageII[i])
    pvalHlp <- pvalStageII[id]
    names(pvalHlp) <- txStageII[id]
    pvalGeneList[[i]] <- pvalHlp
}
padjGeneList <- lapply(pvalGeneList,function(x) p.adjust(x,method="holm"))

sum(p.adjust(dxr$pvalue,method="BH")<.05,na.rm=TRUE) ; sum(unlist(padjGeneList)<alphaAdjusted)

## stage-wise DEXSeq analysis
significantGenes <- names(qvalDxr)[which(qvalDxr<=.05)]
alphaAdjusted <- 0.05*length(significantGenes)/length(qvalDxr)
genesStageII <- dxr$groupID[dxr$groupID%in%significantGenes]
uniqueGenesStageII <- unique(genesStageII)
txStageII <- dxr$featureID[dxr$groupID%in%significantGenes]
pvalStageII <- dxr$pvalue[dxr$featureID%in%txStageII]
pvalGeneList <- list()
for(i in 1:length(uniqueGenesStageII)){
    id <- which(genesStageII==uniqueGenesStageII[i])
    pvalHlp <- pvalStageII[id]
    names(pvalHlp) <- txStageII[id]
    pvalGeneList[[i]] <- pvalHlp
}
padjGeneList <- lapply(pvalGeneList,function(x) p.adjust(x,method="holm"))

sum(p.adjust(dxr$pvalue,method="BH")<.05,na.rm=TRUE) ; sum(unlist(padjGeneList)<alphaAdjusted)


### gene-level ROC
truth_gene <- truth[,c(1,3)]
truth_gene <- truth_gene[!duplicated(truth$gene_id),]
rownames(truth_gene) <- truth_gene[,1]
cobra <- COBRAData(padj = data.frame(kallisto_dexseq = qvalDxr,
                                     row.names = names(qvalDxr),
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth_gene, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "gene_ds_status")
cobraplotGene <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplotGene)
plot_roc(cobraplotGene)

### transcript-level ROC
padjTxDexSeq <- p.adjust(dxr$pvalue,"BH")
truth_tx <- truth[,2:3]
truth_tx$transcript_id <- as.character(truth_tx$transcript_id)
rownames(truth_tx) <- truth_tx[,1]
cobra <- COBRAData(padj = data.frame(kallisto_dexseq = padjTxDexSeq,
                                     row.names = dxr$featureID,
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truthTx, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "tx_ds_status")
cobraplotTx <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplotTx)
plot_roc(cobraplotTx)

### combine in one plot
library(scales)
plot(x=cobraplotGene@fdrtprcurve$FDR,y=cobraplotGene@fdrtprcurve$TPR, type="l",col=2,lwd=1, xlim=c(0,0.8), ylab="True Positive Rate", xlab="False Discovery Rate")
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col=2, pch="o", cex=1.2)
lines(x=cobraplotTx@fdrtprcurve$FDR,y=cobraplotTx@fdrtprcurve$TPR,col=4,lwd=1)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col=4, pch="o", cex=1.2)
abline(v=c(0.01,0.05,0.1,seq(0.1,1,.1)), col=alpha("grey",.8), lty=2)
legend("topleft",c("gene-level","transcript-level"),col=c(2,4),lty=1, bty="n")



### ROC stage-wise testing on transcript level: loop
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
pvalDexSeq <- dxr$pvalue
names(pvalDexSeq) <- dxr$featureID
pvalDexSeq[is.na(pvalDexSeq)] <- 1
padjDexSeq <- p.adjust(pvalDexSeq,method="BH")
tprSW <- vector(length=length(pvalSeq))
fprSW <- vector(length=length(pvalSeq))
fdrSW <- vector(length=length(pvalSeq))
for(id in 1:length(pvalSeq)){
    #print(id)
    ### DEXSeq stage-wise
    significantGenes <- names(qvalDxr)[which(qvalDxr<=pvalSeq[id])]
    if(length(significantGenes)==0){
	tprSW[id]=0
	fprSW[id]=0
	fdrSW[id]=0
	next
    } else {
    alphaAdjusted <- pvalSeq[id]*length(significantGenes)/length(qvalDxr)
    
    genesStageII <- dxr$groupID[dxr$groupID%in%significantGenes]
    uniqueGenesStageII <- unique(genesStageII)
    txStageII <- dxr$featureID[dxr$groupID%in%significantGenes]
    pvalStageII <- dxr$pvalue[dxr$featureID%in%txStageII]
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
    #tprSW[id] <- mean(truth_tx$transcript_id[truth_tx$gene_ds_status==1]%in%positives)
    tprSW[id] <- mean(truthTx$transcript[truthTx$tx_ds_status==1]%in%positives)  
    #fprSW[id] <- mean(truth_tx$transcript_id[truth_tx$gene_ds_status==0]%in%positives)
    fprSW[id] <- mean(truthTx$transcript[truthTx$tx_ds_status==0]%in%positives)
   #truePosId <- truth_tx$gene_ds_status[match(positives,truth_tx$transcript_id)]
    truePosId <- truthTx$tx_ds_status[match(positives,truthTx$transcript)]
    fdrSW[id] <- 1-mean(truePosId)
    }
}

### DEXSeq regular
evalDexSeqRegular <- t(sapply(pvalSeq,function(alpha){
	   positives <- names(padjDexSeq[which(padjDexSeq<alpha)])
	   #tprRegular <- mean(truth_tx$transcript_id[truth_tx$gene_ds_status==1]%in%positives)
	   tprRegular <- mean(truthTx$transcript[truthTx$tx_ds_status==1]%in%positives)	   
	   #fprRegular <- mean(truth_tx$transcript_id[truth_tx$gene_ds_status==0]%in%positives)
	   fprRegular <- mean(truthTx$transcript[truthTx$tx_ds_status==0]%in%positives)	   
	   #truePosId <- truth_tx$gene_ds_status[match(positives,truth_tx$transcript_id)]
	   truePosId <- truthTx$tx_ds_status[match(positives,truthTx$transcript)]	   
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
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3)

## FDR-TPR curve
plot(x=fdrSW,y=tprSW,type="n")
lines(x=fdrSW,y=tprSW,col=2,lwd=2)
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],pch="o",col=2)
lines(x=evalDexSeqRegular[,"fdr"],y=evalDexSeqRegular[,"tpr"],col=3,lwd=2)
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],col=3,pch="o")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3)


### diffSplice in edgeR
condition=factor(rep(c("A","B"),each=3))
design=model.matrix(~condition)

genesData=data.frame(GeneID=geneTx, TxID=rownames(txCount))
d=DGEList(txCount, genes=genesData)
d=calcNormFactors(d)
d=estimateDisp(d, design)
plotBCV(d)
fit=glmFit(d,design)
sp=diffSpliceDGE(fit,coef=2, geneid="GeneID", exonid="TxID")
hist(sp$exon.p.value)
geneQDS <- p.adjust(sp$gene.p.value,"BH")
names(geneQDS) <- rownames(sp$gene.p.value)
sum(geneQDS<.05) ; sum(qvalDxr<.05) #more genes for DexSeq
pTxDex <- p.adjust(dxr$pvalue,"BH")
pTxDex <- pTxDex[!is.na(pTxDex)]
sum(p.adjust(sp$exon.p.value,"BH")<.05) ; sum(p.adjust(pTxDex,"BH")<.05) #more transcripts for DEXSeq

## fdr-tpr
cobra <- COBRAData(padj = data.frame(kallisto_edgeR = geneQDS,
                                     row.names = names(geneQDS),
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(object_to_extend=cobra,
		   padj = data.frame(kallisto_DEXSeq = qvalDxr,
				     row.names = names(qvalDxr),
				     stringsAsFactors=FALSE))
cobra <- COBRAData(truth = truth_gene, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "gene_ds_status")
cobraplotGene <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplotGene)
plot_roc(cobraplotGene)


