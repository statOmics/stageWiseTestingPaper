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
#library(rtracklayer)
## download the GTF file and save it.
#gtfWholeGenome="/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.gtf"
#gtff <- import(gtfWholeGenome)
#gtff2 <- subset(gtff, seqnames == "1")
#export(gtff2, "/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.chr1.gtf", format = "gtf")

## generate tx 2 gene map
cdna_fasta <-"/Users/koenvandenberge/Dropbox/edgeR_zeroinflation/simulated_data_bulk_tx/sim2_human/reference_files/Homo_sapiens.GRCh37.71.cdna.chr1.fa"
gtf <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/Homo_sapiens.GRCh37.71.chr1.gtf"
feature_lengths_file <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/feature_lengths.Rdata"
tx_gene_file <- "/Users/koenvandenberge/PhD_Data/dtu/annotation/tx_gene_map.Rdata"
#calc_lengths_mapping(gtf = gtf, cdna_fasta = cdna_fasta, 
#                       feature_lengths_file = feature_lengths_file, 
#                       tx_gene_file = tx_gene_file) 
load(feature_lengths_file)
load(tx_gene_file) #gene2tx
tx2gene = gene2tx[,c(2,1)]

truth_tx_file <- "/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper/stageWiseTestingPaper/DTU_DTE/simulation/introduction/truth_transcript.txt"
truth_tx <- read.delim(truth_tx_file, header = TRUE, as.is = TRUE)


### derive salmon transcript counts
salmon_basedir <- "/Users/koenvandenberge/PhD_Data/dtu/sim2_human/salmon"

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
## discard genes with only one transcript and not expressed transcripts
txCount <- ceiling(salmon_quant$txCOUNT_sal)
txCount <- txCount[!rowSums(txCount)==0,]
genesWithOneTx <- names(table(tx2gene$gene))[table(tx2gene$gene)==1]
txFromGenesWithOneTx <- tx2gene$transcript[match(genesWithOneTx,tx2gene$gene)]
txCount <- txCount[!rownames(txCount)%in%txFromGenesWithOneTx,]

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
dxr_sal <- DEXSeqResults(dxd_sal)
hist(dxr_sal$pvalue)
qval_dtu_salmon <- perGeneQValue(dxr_sal)

## ROC gene-level analysis
truth_gene_file <- "/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper/stageWiseTestingPaper/DTU_DTE/simulation/introduction/truth_gene.txt"
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
par(bty="l", cex.axis=1.5, cex.lab=1.5, mar=c(5,4.5,4,1))
plot(x=cobraplotGene@fdrtprcurve$FDR,y=cobraplotGene@fdrtprcurve$TPR, type="l",col="black",lwd=2, xlim=c(0,0.7), ylab="True Positive Rate", xlab="False Discovery Proportion")
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col="white", pch=19, cex=1.2)
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col="black", pch="o", cex=1.2)
lines(x=cobraplotTx@fdrtprcurve$FDR,y=cobraplotTx@fdrtprcurve$TPR,col="red",lwd=2)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col="white", pch=19, cex=1.2)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col="red", pch="o", cex=1.2)
abline(v=c(0.01,0.05,0.1,seq(0.1,1,.1)), col=alpha("grey",.8), lty=2)
legend("bottomright",c("gene-level","transcript-level"),col=c("black","red"),lty=1, bty="n", cex=1.5)

