### simulation Human from Genome Biology paper Soneson 2016
baseDir <- "/Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb/hsapiens/diffexpression/non_null_simulation_dge/non_null_simulation/"
files=list.files(baseDir,recursive=TRUE)
fastaFiles <- files[grep(x=files,pattern=".fq")]
sampleNum <- as.numeric(unlist(lapply(strsplit(fastaFiles,split="_"),function(x) x[3])))
fastaFiles <- fastaFiles[order(sampleNum)]
names(fastaFiles) <- rep(paste0("Hs_sample_",1:10),each=2)
kallistoIndex="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/hsapiens/reference_files/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly"
txConversionFile="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/hsapiens/reference_files/KallistoIndex/TranscriptID_conversion.txt"
kallistoDir=paste0(baseDir,"quantifications/kallisto/")
truthFile=paste0(baseDir,"3_truth/simulation_details.txt")
library(edgeR)
library(iCOBRA)
library(scales)
library(dplyr)
library(readr)
library(tximport)


### kallisto quantification
sample <- as.numeric(unlist(lapply(strsplit(names(fastaFiles),split="_"),function(x) x[3])))
fileNames <- unlist(lapply(strsplit(names(fastaFiles),split=".",fixed=TRUE),function(x) x[1]))[seq(1,20,2)]
#for(i in 1:10){
#    pairedFasta <- fastaFiles[sample==i]
#    cmd <- paste0("kallisto quant -i ",kallistoIndex,
#	      " -o ",kallistoDir,fileNames[i],
#	      " -b 30 ",
#	      baseDir,pairedFasta[1]," ",baseDir,pairedFasta[2])
#    message(cmd)
#    system(cmd)
#}

## read in truth and gene names
truth=read.table(truthFile,header=TRUE)
tx2gene <- truth[,c("transcript_id","gene_id")]
tx2gene$gene_id <- as.character(tx2gene$gene_id)
tx2gene$transcript_id <- as.character(tx2gene$transcript_id)
kal2tx=read.table(txConversionFile)
colnames(kal2tx) <- c("kallisto","transcript")


### get kallisto results
files2 <- list.files(kallistoDir)
sampleDirs <- files2[grep(x=files2,pattern="Hs_sample_[123456789]")]
#sort them numerically
sampleDirs <- sampleDirs[order(as.numeric(unlist(lapply(strsplit(sampleDirs,split="_"),function(x) x[3]))))]
txData <- tximport(files=paste0(kallistoDir,"/",sampleDirs,"/abundance.tsv"), type="kallisto", txIn=TRUE, txOut=TRUE, reader=read_tsv, tx2gene=tx2gene)
rownames(txData$counts) <- kal2tx$transcript[match(rownames(txData$counts),kal2tx$kallisto)]
rownames(txData$abundance) <- kal2tx$transcript[match(rownames(txData$abundance),kal2tx$kallisto)]
rownames(txData$length) <- kal2tx$transcript[match(rownames(txData$length),kal2tx$kallisto)]
geneData <- summarizeToGene(txData,tx2gene)

### DTE analysis
condition=factor(rep(0:1,each=5))
design=model.matrix(~condition)
library(edgeR)
alpha=0.05

################################
######## gene level test #######
################################
## tx level analysis edgeR using perGeneQValue to aggregate p-values.
dTx=DGEList(txData$counts)
dTx=calcNormFactors(dTx)
dTx=estimateGLMCommonDisp(dTx,design)
dTx=estimateGLMTrendedDisp(dTx,design)
dTx=estimateGLMTagwiseDisp(dTx,design)
plotBCV(dTx)
fitTx=glmFit(dTx,design)
lrtTx=glmLRT(fitTx,coef=2)
padjTx=p.adjust(lrtTx$table$PValue,"BH")
sum(padjTx<.05) #nr of tx found
genesAll=tx2gene$gene_id[match(rownames(txData$counts),tx2gene$transcript_id)]
significantGenesTx=unique(genesAll[padjTx<.05])
significantTxTx <- rownames(lrtTx)[padjTx<.05]
## using perGeneQValue
source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/perGeneQValue_kvdb.R")
library(Biobase)
pvals=lrtTx$table$PValue
object=list()
object$groupID=genesAll
qvals=perGeneQValue_kvdb(object=object,pvals=pvals)
significantGenesQval=names(qvals)[qvals<=alpha]

#stage-wise analysis with stageR
library(stageR)
pConfirmation=matrix(pvals,ncol=1,dimnames=list(rownames(lrtTx),"transcript"))
stageRObj = stageRTx(pScreen=qvals, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj = stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
head(getSignificantGenes(stageRObj))
head(getSignificantTx(stageRObj))
head(getAdjustedPValues(stageRObj, order=TRUE))


### TPR
#I will combine DTU and DTE status because both result in DTE
truthAll = data.frame(gene_id=truth$gene_id, transcript_id=truth$transcript_id, status=as.numeric(truth$transcript_ds_status | truth$gene_de_status), statusGene=as.numeric(truth$gene_ds_status | truth$gene_de_status))
#compare gene level test on gene level and tx level test on tx level.
mean(truthAll$gene_id[truthAll$status==1]%in%significantGenesQval)
mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxTx)

### FDR
mean(!significantGenesQval%in%truthAll$gene_id[truthAll$status==1])
mean(!significantTxTx%in%truthAll$transcript_id[truthAll$status==1])

### FDR-TPR curve
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
performanceData=data.frame(tprQval=NA, tprTxTx=NA, tprTxSW=NA, fdrQval=NA, fdrTx=NA, fdrTxSW=NA)
stageRObjSim = stageRTx(pScreen=qvals, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)

for(i in 1:length(pvalSeq)){
    print(i)
	   alpha=pvalSeq[i]
	   significantGenesQval <- names(qvals)[qvals<=alpha]
	   significantTxTx <- rownames(lrtTx)[padjTx<=alpha]
	   ## stage-wise follow up
	   	stageRObjRes = stageWiseAdjustment(stageRObjSim, method="dte", alpha=alpha)
		significantTxSW = rownames(getSignificantTx(stageRObjRes))
	  

	   tprQval=mean(truthAll$gene_id[truthAll$status==1]%in%significantGenesQval)
	   tprTxTx=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxTx)
	   tprTxSW=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxSW)	  

	   fdrQval=mean(!significantGenesQval%in%truthAll$gene_id[truthAll$status==1])
	   #ofdrQval=mean(hlpOfdr$ofdrFalse)
	   fdrTxTx=mean(!significantTxTx%in%truthAll$transcript_id[truthAll$status==1])
	   fdrTxSW=mean(!significantTxSW%in%truthAll$transcript_id[truthAll$status==1])

	   performanceData[i,] <- c(tprQval, tprTxTx, tprTxSW, fdrQval, fdrTxTx, fdrTxSW)
}
performanceDataHuman=performanceData
save(performanceDataHuman,file="~/performaneDataHuman.rda")

### FDR-TPR curve
par(bty="l", mar=c(5,4.5,4,1))
plot(x=performanceData$fdrQval, y=performanceData$tprQval, type="l", ylab="True Positive Rate", xlab="False Discovery Proportion", ylim=c(0,1), lwd=2, main="", cex.lab=1.5, cex.axis=1.5)
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)],col="black")
lines(x=performanceData$fdrTx, y=performanceData$tprTxTx, col=2, lwd=2)
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)],col=2)
lines(x=performanceData$fdrTxSW, y=performanceData$tprTxSW, col=3, lwd=2)
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)],col=3)
legend("bottomright",c("gene-level","tx-level","tx-level stage-wise"),lty=1,col=1:3, cex=1.25, bty="n")


