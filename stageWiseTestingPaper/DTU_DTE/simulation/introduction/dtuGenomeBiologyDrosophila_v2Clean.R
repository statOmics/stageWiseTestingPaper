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

### check 2 most abundant tx
simFilesDir <- "/Users/koenvandenberge/PhD_Data/dtu/E-MTAB-3766/"
simFiles=list.files(simFilesDir)
firstSampleFiles <- simFiles[grep(x=simFiles,pattern="Dm_sample[1-3]")]
sample1Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[1]),header=TRUE)
sample2Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[2]),header=TRUE)
sample3Tpm <- read.table(paste0(simFilesDir,firstSampleFiles[3]),header=TRUE)

group1Tpm <- sample1Tpm
group1Tpm$TPM2 <- sample2Tpm$TPM
group1Tpm$TPM3 <- sample3Tpm$TPM
genesWithOneTx <- names(table(sample1Tpm$gene_id)[table(sample1Tpm$gene_id)==1])
group1Tpm <- group1Tpm[!group1Tpm$gene_id%in%genesWithOneTx,]
hlp=group_by(group1Tpm,by=gene_id)
highestExpressedTxPerGene <- dplyr::summarize(hlp, max1=transcript_id[order(apply(cbind(TPM,TPM2,TPM3),1,mean),decreasing=TRUE)[1]], max2=transcript_id[order(apply(cbind(TPM,TPM2,TPM3),1,mean),decreasing=TRUE)[2]])
head(highestExpressedTxPerGene)

truthTx <- data.frame(transcript=truth_tx[,1],tx_ds_status=0,row.names=rownames(truth_tx))
simulatedDTUGenes <- rownames(truth_gene[truth_gene$gene_ds_status==1,])
length(simulatedDTUGenes) #should be 1000
highestExpressedTxForSignificantGenes <- highestExpressedTxPerGene[highestExpressedTxPerGene$by%in%simulatedDTUGenes,]
dtuTx <- as.character(unlist(c(as.data.frame(highestExpressedTxForSignificantGenes)[,2:3])))
truthTx[dtuTx,"tx_ds_status"]=1


### DEXSeq analysis
genesWithOneTx <- names(table(tx2gene$gene))[table(tx2gene$gene)==1]
txFromGenesWithOneTx <- tx2gene$transcript[match(genesWithOneTx,tx2gene$gene)]
txCount <- ceiling(data)
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


### gene-level
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

### transcript-level
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
par(bty="l", cex.axis=1.5, cex.lab=1.5, mar=c(5,4.5,4,1))
plot(x=cobraplotGene@fdrtprcurve$FDR,y=cobraplotGene@fdrtprcurve$TPR,col=2,lwd=1, xlim=c(0,0.8), ylab="True Positive Rate", xlab="False Discovery Proportion", type="n")
abline(v=c(0.01,0.05,0.1,seq(0.1,1,.1)), col=alpha("grey",.8), lty=2)
lines(x=cobraplotGene@fdrtprcurve$FDR,y=cobraplotGene@fdrtprcurve$TPR,col="black",lwd=2)
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col="white", pch=19, cex=1.2)
points(x=cobraplotGene@fdrtpr$FDR,y=cobraplotGene@fdrtpr$TPR, col="black", pch="o", cex=1.2)
lines(x=cobraplotTx@fdrtprcurve$FDR,y=cobraplotTx@fdrtprcurve$TPR,col="red",lwd=2)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col="white", pch=19, cex=1.2)
points(x=cobraplotTx@fdrtpr$FDR,y=cobraplotTx@fdrtpr$TPR, col="red", pch="o", cex=1.2)
legend("bottomright",c("gene-level","transcript-level"),col=c("black","red"),lty=1, bty="n", cex=1.5)


