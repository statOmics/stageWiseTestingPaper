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
adjustShafferDTE <- function(p){
    p <- sort(p)
    n <- length(p)
    pAdj <- vector(length=n)
    adjustments <- c(n-1) #Shaffer. if only tx, p-values is zero.
    if(n>1){
	adjustments2 <- c(n-(2:n)+1) #Holm
	adjustments <- c(adjustments,adjustments2)
    }
    pAdj <- p*adjustments    
    pAdj[pAdj>1] <- 1
    # check monotone increase of adjusted p-values
    if(any(diff(pAdj)<0)){
	id <- which(diff(pAdj)<0)
	pAdj[id+1] <- pAdj[id]
    }
    return(pAdj)
}


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
## edgeR using gene level data: fails to pick up DTU genes with equal output between conditions.
d=DGEList(geneData$counts)
normMat <- geneData$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(d$counts/normMat)) + log(colSums(d$counts/normMat))
d$offset <- t(t(log(normMat))+o)
d=estimateGLMCommonDisp(d,design)
d=estimateGLMTrendedDisp(d,design)
d=estimateGLMTagwiseDisp(d,design)
plotBCV(d)
fit=glmFit(d,design)
lrt=glmLRT(fit,coef=2)
significantGenes=rownames(lrt)[p.adjust(lrt$table$PValue,"BH")<alpha]
length(significantGenes)

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
source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/perGeneQValue_kvdb.R"
library(Biobase)
pvals=lrtTx$table$PValue
object=list()
object$groupID=genesAll
qvals=perGeneQValue_kvdb(object=object,pvals=pvals)
significantGenesQval=names(qvals)[qvals<=alpha]
length(significantGenesQval)
## stage-wise follow up
genesUnique <- unique(genesAll)
pvalList=list()
for(i in 1:length(genesUnique)){
    id=which(genesAll==genesUnique[i])
    pvalList[[i]]=pvals[id]
    names(pvalList[[i]]) <- rownames(lrtTx)[id]
}
names(pvalList) <- genesUnique
pvalList=pvalList[significantGenesQval] #only select selected genes
padjList=lapply(pvalList, function(x) p.adjust(x,method="holm"))
## faster follow up for loop
hlp=group_by(data.frame(pvals=pvals,transcript=rownames(lrtTx)),by=factor(genesAll))
hlp=mutate(hlp,padj=p.adjust(pvals,"holm"))


## using Fisher's method
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE) #from Mike Love
dat=data.frame(gene=genesAll,pval=pvals)
hlp=group_by(dat,by=genesAll)
geneFisherP=summarize(hlp,fisherP=fishersMethod(pval))
geneFisherPadj <- p.adjust(geneFisherP$fisherP)
significantGenesFisher <- geneFisherP$by[geneFisherPadj<.05]
length(significantGenesFisher)

### number of genes found
length(significantGenes) #edgeR gene level data
length(significantGenesQval) #tx-level analysis, perGeneQValue aggregation
length(significantGenesFisher) #tx-level analysis, Fisher's method aggregation
length(significantGenesTx) #tx-level analysis: how many genes

### TPR
#I will combine DTU and DTE status because both result in DTE
truthAll = data.frame(gene_id=truth$gene_id, transcript_id=truth$transcript_id, status=as.numeric(truth$transcript_ds_status | truth$gene_de_status), statusGene=as.numeric(truth$gene_ds_status | truth$gene_de_status))
#as in F1000 paper: compare gene level test on gene level and tx level test on tx level.
mean(truthAll$gene_id[truthAll$status==1]%in%significantGenesQval)
mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxTx)

### FDR
mean(!significantGenesQval%in%truthAll$gene_id[truthAll$status==1])
mean(!significantTxTx%in%truthAll$transcript_id[truthAll$status==1])

### FDR-TPR curve
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
performanceData=data.frame(tprQval=NA, tprTxTx=NA, tprTxSW=NA, fdrQval=NA, fdrTx=NA, fdrTxSW=NA, ofdrQval=NA)
for(i in 1:length(pvalSeq)){
    print(i)
	   alpha=pvalSeq[i]
	   significantGenesQval <- names(qvals)[qvals<=alpha]
	   significantTxTx <- rownames(lrtTx)[padjTx<=alpha]
	   ## stage-wise follow up
	   alphaAdjusted=alpha*length(significantGenesQval)/length(genesUnique)	   
	   selectedGenesId=genesAll%in%significantGenesQval
	   genesAllLoop=genesAll[selectedGenesId]
	   pvalsLoop=pvals[selectedGenesId]
	   txnamesLoop=rownames(lrtTx)[selectedGenesId]
genesUniqueLoop=unique(genesAll[selectedGenesId])
pvalList=list()
for(j in 1:length(genesUniqueLoop)){
    id=which(genesAllLoop==genesUniqueLoop[j])
    pvalList[[j]]=pvalsLoop[id]
    names(pvalList[[j]]) <- txnamesLoop[id]
}
names(pvalList) <- genesUniqueLoop
pvalList=pvalList[significantGenesQval] #only select selected genes
padjList=lapply(pvalList, function(x) p.adjust(x,method="holm"))
padjListShaffer=lapply(pvalList, function(x) adjustShafferDTE(x))
significantTxSW <- unlist(lapply(strsplit(names(unlist(padjListShaffer)[unlist(padjListShaffer)<=alphaAdjusted]),split=".",fixed=TRUE), function(x) x[2]))
	   #hlp=group_by(data.frame(pvals=pvalsLoop,transcript=txnamesLoop),by=factor(genesAllLoop))
	   #hlp=mutate(hlp,padj=p.adjust(pvals,"holm"), falseGene=as.character(by)%in%truthAll$gene_id[truthAll$statusGene==0])
	   #hlp=mutate(hlp,padj=adjustShafferDTE(pvals), falseGene=as.character(by)%in%truthAll$gene_id[truthAll$statusGene==0])	   
	   #hlp$falseTx=FALSE
	   #hlp$falseTx[hlp$padj<=alphaAdjusted] <- hlp$transcript[hlp$padj<=alphaAdjusted]%in%truthAll$transcript_id[truthAll$status==0]
	   #hlpOfdr=summarize(hlp,ofdrFalse=any(c(falseGene,falseTx)))
	   #significantTxSW <- as.character(hlp$transcript[hlp$padj<=alphaAdjusted])

	   tprQval=mean(truthAll$gene_id[truthAll$status==1]%in%significantGenesQval)
	   tprTxTx=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxTx)
	   tprTxSW=mean(truthAll$transcript_id[truthAll$status==1]%in%significantTxSW)	  

	   fdrQval=mean(!significantGenesQval%in%truthAll$gene_id[truthAll$status==1])
	   #ofdrQval=mean(hlpOfdr$ofdrFalse)
	   fdrTxTx=mean(!significantTxTx%in%truthAll$transcript_id[truthAll$status==1])
	   fdrTxSW=mean(!significantTxSW%in%truthAll$transcript_id[truthAll$status==1])

	   performanceData[i,] <- c(tprQval, tprTxTx, tprTxSW, fdrQval, fdrTxTx, fdrTxSW, ofdrQval)
}

### regular FDR-TPR curve
plot(x=performanceData$fdrQval, y=performanceData$tprQval, type="l", ylab="TPR", xlab="FDR", ylim=c(0,1), lwd=2, main="Human DTE analysis", cex.lab=1.5, cex.axis=1.5, bty="l")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)],col="black")
lines(x=performanceData$fdrTx, y=performanceData$tprTxTx, col=2, lwd=2)
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)],col=2)
lines(x=performanceData$fdrTxSW, y=performanceData$tprTxSW, col=3, lwd=2)
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)],col=3)
legend("bottomright",c("gene-level","tx-level","tx-level stage-wise"),lty=1,col=1:3, cex=1.5, bty="n")

### OFDR-TPR curve: only gene level curve changes
plot(x=performanceData$ofdrQval, y=performanceData$tprQval, type="l", ylab="TPR", xlab="OFDR", ylim=c(0,1), lwd=2)
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceData$ofdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceData$ofdrQval[c(508,516,526)],y=performanceData$tprQval[c(508,516,526)],col="black")
lines(x=performanceData$fdrTx, y=performanceData$tprTxTx, col=2, lwd=2)
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTx[c(508,516,526)],y=performanceData$tprTxTx[c(508,516,526)],col=2)
lines(x=performanceData$fdrTxSW, y=performanceData$tprTxSW, col=3, lwd=2)
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceData$fdrTxSW[c(508,516,526)],y=performanceData$tprTxSW[c(508,516,526)],col=3)
legend("bottomright",c("gene-level","tx-level","tx-level stage-wise"),lty=1,col=1:3)
## conclusion: OFDR control is somewhat worse than FDR control on the gene level, which makes sense because OFDR control is more stringent by involving the transcripts in the FDR calculation for the gene level test while a regular FDR control only looks at the screening hypothesis..












# as stated in F1000 there is indeed no loss in FDR control when aggregating tx level DTE p-values to the gene level but the FDR is massive.

















