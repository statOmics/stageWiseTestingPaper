###############################################################
#### kallisto data: RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. Downloaded from http://lair.berkeley.edu/eswaran12/?
###############################################################
#####################################################
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
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
dataMessy <- read.csv(file="kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
#dataMessy <- read.csv(file="kallisto_table_normalized_filtered.csv",header=TRUE)
dataMessy2 <- dataMessy[,c("target_id","est_counts","sample")]
dataClean <- spread(dataMessy2,key=sample,value=est_counts)
#rm(dataMessy)
rownames(dataClean) <- dataClean[,"target_id"]
dataClean <- dataClean[,-1]
## downloaded gene to tx conversion table from biomaRt website (http://www.ensembl.org/biomart/). Database used is Homo sapiens genes (GRCh38.p5)
tx2gene <- read.delim("ensemblGeneTxGeneName.txt",header=TRUE)
mean(is.na(match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID))) #very few tx dont have a match

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
padjTx=p.adjust(lrt$table$PValue,"BH")
sum(padjTx<.05) #nr of tx found


### stage-wise testing
source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/perGeneQValue_kvdb.R")
library(Biobase)
pvals=lrt$table$PValue
genesAll=tx2gene$Ensembl.Gene.ID[match(rownames(d$counts),tx2gene$Ensembl.Transcript.ID)]
object=list()
object$groupID=genesAll
qvals=perGeneQValue_kvdb(object=object,pvals=pvals)
alpha=.05
significantGenesQval=names(qvals)[qvals<=alpha]
length(significantGenesQval)
alphaAdjusted=alpha*length(significantGenesQval)/length(qvals)
## stage-wise follow up
genesUnique <- unique(genesAll)
pvalList=list()
for(i in 1:length(genesUnique)){
    id=which(genesAll==genesUnique[i])
    pvalList[[i]]=pvals[id]
    names(pvalList[[i]]) <- rownames(lrt)[id]
}
names(pvalList) <- genesUnique
pvalList=pvalList[significantGenesQval] #only select selected genes
padjList=lapply(pvalList, function(x) p.adjust(x,method="holm"))
padjListShaffer=lapply(pvalList, function(x) adjustShafferDTE(x))
padjListBH=lapply(pvalList, function(x) p.adjust(x,method="BH"))

#number of transcripts
significantTranscriptsTx <- rownames(lrt)[padjTx<=.05]
significantTranscriptsStageWise <- names(unlist(padjListShaffer))[unlist(padjListShaffer)<=alphaAdjusted]
significantTranscriptsStageWise <- unlist(lapply(strsplit(significantTranscriptsStageWise,split=".",fixed=TRUE), function(x) x[2]))
length(significantTranscriptsTx) ; length(significantTranscriptsStageWise)

#number of genes
genesSW=significantGenesQval
genesTx=unique(as.character(tx2gene[match(rownames(lrt)[padjTx<=.05],tx2gene$Ensembl.Transcript.ID),"Ensembl.Gene.ID"]))
length(genesSW) ; length(genesTx)

#most significant gene
which.min(qvals)
d$counts[which(genesAll=="ENSG00000012822"),] 
d$counts[which.min(pvals),]

