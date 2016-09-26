###############################################################
#### kallisto data: RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. Downloaded from http://lair.berkeley.edu/eswaran12/?
###############################################################
#####################################################
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
#dataMessy <- read.csv(file="kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
dataMessy <- read.csv(file="kallisto_table_normalized_filtered.csv",header=TRUE)
dataMessy <- dataMessy[,c("target_id","tpm","sample")]
dataClean <- spread(dataMessy,key=sample,value=tpm)
rm(dataMessy)
rownames(dataClean) <- dataClean[,"target_id"]
dataClean <- dataClean[,-1]
## downloaded gene to tx conversion table from biomaRt website (http://www.ensembl.org/biomart/). Database used is Homo sapiens genes (GRCh38.p5)
tx2gene <- read.delim("ensemblGeneTxGeneName.txt",header=TRUE)
mean(is.na(match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID))) #very few tx dont have a match
adjustShaffer <- function(p){
    p <- sort(p)
    n <- length(p)
    pAdj <- vector(length=n)
    adjustments <- c(n-2,n-2) #Shaffer. if only two tx, p-value becomes zero
    if(n>2){
	adjustments2 <- c(n-(3:n)+1) #Holm
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
dataClean <- round(dataClean)
## remove tx without gene match
dataClean <- dataClean[!is.na(match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)),]
## remove all zero rows
dataClean <- dataClean[!rowSums(dataClean)==0,]
## remove genes with only one tx
geneTable <- table(as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]))
genesWithOneTx <- names(geneTable)[geneTable==1]
txFromGenesWithOneTx <- tx2gene$Ensembl.Transcript.ID[tx2gene$Ensembl.Gene.ID%in%genesWithOneTx]
dataClean <- dataClean[!rownames(dataClean)%in%as.character(txFromGenesWithOneTx),]

txGeneData = as.data.frame(cbind(rownames(dataClean),as.character(tx2gene$Ensembl.Transcript.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]),as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)])))
colnames(txGeneData)=c("tx","transcript","gene")
barplot(table(table(txGeneData$gene)), main="Distribution of number of tx per gene")

### regular DEXSeq analysis
library(DEXSeq)
dxd <- DEXSeqDataSet(countData = dataClean,
                         sampleData = sampleData,
                         design = ~ sample + exon + patient + condition:exon,
                         featureID = rownames(dataClean),
                         groupID = as.character(txGeneData$gene))
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, fitType="local")
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon + patient)
#save(dxd,file="dxdProstate_v2.RData")
#load("dxdProstate_v2.RData")
plotDispEsts(dxd)
dxr <- DEXSeqResults(dxd)
hist(dxr$pvalue)
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
padjGeneListShaffer <- lapply(pvalGeneList, function(x) adjustShaffer(x))
padjGeneListBH <- lapply(pvalGeneList, function(x) p.adjust(x,"BH"))

## nr of tx found for both methods: more in a tx-level analysis
sum(p.adjust(dxr$pvalue,"BH")<.05,na.rm=TRUE) ; sum(unlist(padjGeneListShaffer)<alphaAdjusted) ; sum(unlist(padjGeneListBH)<alphaAdjusted)
significantTranscriptsTx <- dxr$featureID[which(p.adjust(dxr$pvalue,"BH")<.05)]
significantTranscriptsStageWise <- names(unlist(padjGeneListShaffer))[unlist(padjGeneListShaffer)<alphaAdjusted]
## nr of genes found for both methods: more in a stage-wise analysis
genesSW <- names(qvalDxr)[qvalDxr<.05]
genesTx <- unique(dxr$groupID[p.adjust(dxr$pvalue,"BH")<.05])
genesTx <- genesTx[!is.na(genesTx)]
length(genesSW) ; length(genesTx)

## genes/tx found in only one method
genesOnlyStage <- genesSW[!genesSW%in%genesTx] #genes only found in a stage-wise analysis
genesOnlyTx <- genesTx[!genesTx%in%genesSW]
names(padjGeneListShaffer)=uniqueGenesStageII
txOnlyTx <- significantTranscriptsTx[!significantTranscriptsTx%in%significantTranscriptsStageWise]
genesFromOnlyTxTx <- unique(dxr$groupID[dxr$featureID%in%txOnlyTx])
mean(genesFromOnlyTxTx%in%genesSW) #nearly all genes that contain transcripts that are only significant in a tx-level analysis are also found in a SW analysis.'


