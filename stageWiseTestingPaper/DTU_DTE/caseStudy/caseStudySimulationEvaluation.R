###############################################################
#### kallisto data: RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. Downloaded from http://lair.berkeley.edu/eswaran12/?
###############################################################
#####################################################
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
dataMessy <- read.csv(file="kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
dataMessy <- dataMessy[,c("target_id","est_counts","sample")]
dataClean <- spread(dataMessy,key=sample,value=est_counts)
rm(dataMessy)
rownames(dataClean) <- dataClean[,"target_id"]
dataClean <- dataClean[,-1]
library(biomaRt)
mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host='mar2016.archive.ensembl.org', dataset="hsapiens_gene_ensembl")
tx2gene = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id"), filters="ensembl_transcript_id", values=rownames(dataClean), mart=mart, bmHeader=TRUE, uniqueRows=TRUE)
colnames(tx2gene) <- c("Ensembl.Transcript.ID","Ensembl.Gene.ID")


## metadata
metaData <- read.table("sampleDataRelationship.txt",header=TRUE,sep="\t")
assays <- metaData$Assay.Name
runs <- as.character(metaData$Comment.ENA_RUN.)[seq(1,length(assays),2)]
samples=gsub(x=assays,pattern="_[1-2]",replacement="")[seq(1,length(assays),2)]
patient=factor(sapply(samples,function(x) substr(x,1,nchar(x)-1)))
condition=factor(sapply(samples,function(x) substr(x,nchar(x),nchar(x))))
dataClean <- dataClean[,match(runs,colnames(dataClean))] #same ordering as metadata
sampleData <- data.frame(condition=condition,patient=patient)
rownames(sampleData)=colnames(dataClean)


### clean up
dataClean <- ceiling(dataClean)
## remove tx without gene match
dataClean <- dataClean[!is.na(match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)),]
## more stringent filtering to guarantee some expression for swapping tx.
#dataClean <- dataClean[!rowSums(dataClean)==0,]
dataClean <- dataClean[rowSums(dataClean>5)>=6,]
## remove genes with only one tx
geneTable <- table(as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]))
genesWithOneTx <- names(geneTable)[geneTable==1]
txFromGenesWithOneTx <- tx2gene$Ensembl.Transcript.ID[tx2gene$Ensembl.Gene.ID%in%genesWithOneTx]
dataClean <- dataClean[!rownames(dataClean)%in%as.character(txFromGenesWithOneTx),]

txGeneData = as.data.frame(cbind(rownames(dataClean),as.character(tx2gene$Ensembl.Transcript.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]),as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)])))
colnames(txGeneData)=c("tx","transcript","gene")
barplot(table(table(txGeneData$gene)), main="Distribution of number of tx per gene")

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

#randomly assign disease and control samples stratified within patient.
set.seed(45)
for(i in seq(1,28,2)) condition[i:(i+1)] = sample(condition[i:(i+1)])
#select 1000 genes. within a gene, sample a nr of tx to swap from a Bi(n=1,p=1/3).
#get the rows of the tx and shuffle them for inducing DTE/DTU.
library(dplyr)
deGenes <- sample(unique(txGeneData$gene),1000,replace=FALSE)
deGeneTx <- txGeneData[txGeneData$gene%in%deGenes,] %>% 
    group_by(gene) %>% 
    filter(transcript%in%sample(transcript,max(2,rbinom(size=n(),n=1,prob=1/3)))) %>%
    mutate(row=match(transcript,rownames(dataClean))) %>%
    mutate(rowSampled={
	       samp = sample(row)
	       while(any(samp==row)) samp=sample(row)
	       samp
})

#swap expression between tx for tumor samples
dataSwap <- dataClean
dataSwap[deGeneTx$row,condition=="T"] <- dataClean[deGeneTx$rowSampled,condition=="T"]

#### DTE analysis
library(edgeR)
d=DGEList(dataSwap)
d=calcNormFactors(d)
design=model.matrix(~patient+condition)
d=estimateGLMCommonDisp(d,design)
d=estimateGLMTrendedDisp(d,design)
d=estimateGLMTagwiseDisp(d,design)
plotBCV(d)
fit=glmFit(d,design)
lrt=glmLRT(fit,coef="conditionT")
hist(lrt$table$PValue)
padjTx=p.adjust(lrt$table$PValue,"BH")

### DTE: stage-wise testing
source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/perGeneQValue_kvdb.R")
library(Biobase)
library(stageR)
pvals=lrt$table$PValue
genesAll=tx2gene$Ensembl.Gene.ID[match(rownames(d$counts),tx2gene$Ensembl.Transcript.ID)]
object=list()
object$groupID=genesAll
qvals=perGeneQValue_kvdb(object=object,pvals=pvals)
pConfirmation=matrix(pvals,ncol=1,dimnames=list(rownames(lrt),"transcript"))
pScreen=qvals

stageRObj = stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj = stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
dim(getSignificantGenes(stageRObj))
dim(getSignificantTx(stageRObj))


### Performance
truthAll=data.frame(gene_id=txGeneData$gene, transcript_id=txGeneData$transcript, status=0)
truthAll$status[deGeneTx$row]=1
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
performanceData=data.frame(tprQval=NA, tprTxTx=NA, tprTxSW=NA, fdrQval=NA, fdrTx=NA, fdrTxSW=NA)
stageRObjSim = stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
for(i in 1:length(pvalSeq)){
    print(i)
	   alpha=pvalSeq[i]
	   significantGenesQval <- names(qvals)[qvals<=alpha]
	   significantTxTx <- rownames(lrt)[padjTx<=alpha]
	   ## stage-wise follow up
		stageRObjRes = stageWiseAdjustment(stageRObjSim, method="dte", alpha=alpha)
		significantTxSW = rownames(getSignificantTx(stageRObjRes))
	  
	   tprQval=mean(truthAll$gene_id[truthAll$status==1]%in%significantGenesQval)
	   tprTxTx=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxTx)
	   tprTxSW=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxSW)	  

	   fdrQval=mean(!significantGenesQval%in%truthAll$gene_id[truthAll$status==1])
	   fdrTxTx=mean(!significantTxTx%in%truthAll$transcript_id[truthAll$status==1])
	   fdrTxSW=mean(!significantTxSW%in%truthAll$transcript_id[truthAll$status==1])

	   performanceData[i,] <- c(tprQval, tprTxTx, tprTxSW, fdrQval, fdrTxTx, fdrTxSW)
}
save(performanceData,file="performanceDataDTECaseStudySimulation.rda")


### FDR-TPR curve: dteResultsRealDataSimulation.pdf
par(mar=c(4,5,4,1), cex.axis=1.5)
plot(x=performanceData$fdrQval, y=performanceData$tprQval, type="l", ylab="True Positive Rate", xlab="False Discovery Proportion", ylim=c(0,1), lwd=3, bty="l", main="", cex.lab=1.5, col="darkseagreen")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)],col="darkseagreen")
lines(x=performanceData$fdrTx, y=performanceData$tprTxTx, col="steelblue", lwd=3)
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)],col="steelblue")
lines(x=performanceData$fdrTxSW, y=performanceData$tprTxSW, col="orange", lwd=3)
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)],col="orange")
legend("bottomright",c("gene-level","transcript-level","transcript-level stage-wise"),lty=1,col=c("darkseagreen","steelblue","orange"), bty="n", cex=1.25)


### DTU analysis
library(DEXSeq)
sampleData$condition=condition
dxd <- DEXSeqDataSet(countData = dataSwap,
                         sampleData = sampleData,
                         design = ~ sample + exon + patient + condition:exon,
                         featureID = rownames(dataSwap),
                         groupID = as.character(txGeneData$gene))
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, fitType="local")
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon + patient)
#save(dxd,file="dxdCaseStudySimulationEvaluation.RData")
load("dxdCaseStudySimulationEvaluation.RData")
plotDispEsts(dxd)
dxr <- DEXSeqResults(dxd)
hist(dxr$pvalue)
qvalDxr <- perGeneQValue(dxr)


## stage-wise DEXSeq analysis
library(stageR)
pScreen=qvalDxr
pScreen[is.na(pScreen)]=1 #filtered p-values
pConfirmation <- matrix(dxr$pvalue,ncol=1,dimnames=list(dxr$featureID,"transcript"))
tx2gene <- as.data.frame(cbind(dxr$featureID,dxr$groupID))
rowsNotFiltered=tx2gene[,2]%in%names(qvalDxr)
pConfirmation=matrix(pConfirmation[rowsNotFiltered,],ncol=1,dimnames=list(dxr$featureID[rowsNotFiltered],"transcript"))
tx2gene <- tx2gene[rowsNotFiltered,]

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
dim(getSignificantGenes(stageRObj))
dim(getSignificantTx(stageRObj)) 



### performance stage-wise testing on transcript level: loop
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
pvalDexSeq <- dxr$pvalue
names(pvalDexSeq) <- dxr$featureID
pvalDexSeq[is.na(pvalDexSeq)] <- 1
padjDexSeq <- p.adjust(pvalDexSeq,method="BH")
tprSW <- vector(length=length(pvalSeq))
fprSW <- vector(length=length(pvalSeq))
fdrSW <- vector(length=length(pvalSeq))
stageRObjSim <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
truth=data.frame(gene_id=txGeneData$gene, transcript_id=txGeneData$transcript, status=0)
truth$status[deGeneTx$row]=1
truth_tx=truth

for(id in 1:length(pvalSeq)){
    print(id)
    ### DEXSeq stage-wise
    significantGenes <- names(qvalDxr)[which(qvalDxr<=pvalSeq[id])]
    if(length(significantGenes)==0){
	tprSW[id]=0
	fprSW[id]=0
	fdrSW[id]=0
	next
    } else {
      stageRRes <- stageWiseAdjustment(stageRObjSim, method="dtu", alpha=pvalSeq[id])
      positives <- rownames(getSignificantTx(stageRRes))

    tprSW[id] <- mean(truth_tx$transcript_id[truth_tx$status==1]%in%positives)  
    fprSW[id] <- mean(truth_tx$transcript_id[truth_tx$status==0]%in%positives)
    truePosId <- truth_tx$status[match(positives,truth_tx$transcript_id)]
    fdrSW[id] <- 1-mean(truePosId)
    }
}
evalDtuTxSW=cbind(tprSW,fprSW,fdrSW)
#save(evalDtuTxSW,file="evalDtuTxSW.rda")

### DEXSeq regular
evalDexSeqRegular <- t(sapply(pvalSeq,function(alpha){
	   positives <- names(padjDexSeq[which(padjDexSeq<alpha)])
	   tprRegular <- mean(truth_tx$transcript_id[truth_tx$status==1]%in%positives)	   
	   fprRegular <- mean(truth_tx$transcript_id[truth_tx$status==0]%in%positives)	   
	   truePosId <- truth_tx$status[match(positives,truth_tx$transcript_id)]	   
	   fdrRegular <- 1-mean(truePosId)
	   return(cbind(tprRegular=tprRegular, fprRegular=fprRegular, fdrRegular=fdrRegular))
				   }))
colnames(evalDexSeqRegular) <- c("tpr","fpr","fdr")


### gene-level 
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
truth=data.frame(gene_id=txGeneData$gene, transcript_id=txGeneData$transcript, status=0)
truth[truth$gene_id%in%deGeneTx$gene,"status"]=1
truth_gene=truth
positiveGenes <- unique(as.character(truth_gene$gene_id[truth_gene$status==1]))
negativeGenes <- unique(as.character(truth_gene$gene_id[truth_gene$status==0]))
evalGeneLevelSW <- t(sapply(pvalSeq, function(alpha){
  significantGenes <- names(qvalDxr)[qvalDxr<=alpha]
  tpr <- mean(positiveGenes%in%significantGenes)
  fpr <- mean(negativeGenes%in%significantGenes)
  fdr <- mean(significantGenes%in%negativeGenes)
  return(cbind(tpr,fpr,fdr))
}))
colnames(evalGeneLevelSW) <- c("tpr","fpr","fdr")


## FDR-TPR curve

### combine gene level, tx level and stage-wise analysis in one FDR-TPR plot.
#gene level
par(mar=c(5,4,4,1)+0.2, cex.axis=1.5)
plot(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"], type="n", xlab="False Discovery Proportion", ylab="True Positive Rate", ylim=c(0.2,1), bty="l", main="", cex.lab=1.5, xlim=c(0,1))
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.8),lty=2)
lines(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"], col="darkseagreen", lwd=3)
points(x=evalGeneLevelSW[c(508,516,526),"fdr"],y=evalGeneLevelSW[c(508,516,526),"tpr"],pch=19,col="white")
points(x=evalGeneLevelSW[c(508,516,526),"fdr"],y=evalGeneLevelSW[c(508,516,526),"tpr"],col="darkseagreen")
#transcript level
lines(x=c(evalDexSeqRegular[,"fdr"],0.9672017121),y=c(evalDexSeqRegular[,"tpr"],1),col="steelblue",lwd=3)
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],pch=19,col="white")
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],col="steelblue")
#stage-wise transcript level
lines(x=fdrSW,y=tprSW, lwd=3, col="orange")
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],pch=19,col="white")
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],col="orange")
legend("bottomright",c("gene-level","transcript-level","transcript-level stage-wise"),lty=1,col=c("darkseagreen","steelblue","orange"), bty="n", cex=1.5, lwd=2)








