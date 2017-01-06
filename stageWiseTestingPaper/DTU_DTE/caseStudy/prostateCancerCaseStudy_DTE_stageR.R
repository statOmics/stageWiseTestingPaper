###############################################################
#### kallisto data: RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. Downloaded from http://lair.berkeley.edu/eswaran12/?
###############################################################
#####################################################
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
dataMessy <- read.csv(file="kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
dataMessy2 <- dataMessy[,c("target_id","est_counts","sample")]
dataClean <- spread(dataMessy2,key=sample,value=est_counts)
#rm(dataMessy)
rownames(dataClean) <- dataClean[,"target_id"]
dataClean <- dataClean[,-1]
## downloaded gene to tx conversion table from biomaRt website (http://www.ensembl.org/biomart/). Database used is Homo sapiens genes (GRCh38.p5)
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
## remove all zero rows
dataClean <- dataClean[!rowSums(dataClean)==0,]

txGeneData = as.data.frame(cbind(rownames(dataClean),as.character(tx2gene$Ensembl.Transcript.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]),as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)])))
colnames(txGeneData)=c("tx","transcript","gene")
barplot(table(table(txGeneData$gene)), main="Distribution of number of tx per gene")

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

### DTE analysis
library(edgeR)
d=DGEList(dataClean)
d=calcNormFactors(d)
design=model.matrix(~patient+condition)
d=estimateGLMCommonDisp(d,design)
d=estimateGLMTrendedDisp(d,design)
d=estimateGLMTagwiseDisp(d,design)
plotBCV(d)
fit=glmFit(d,design)
lrt=glmLRT(fit,coef="conditionT")
hist(lrt$table$PValue)

### stage-wise testing
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



