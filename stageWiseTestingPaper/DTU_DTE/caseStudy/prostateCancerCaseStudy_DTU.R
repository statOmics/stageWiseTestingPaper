###############################################################
#### kallisto data: RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. Downloaded from http://lair.berkeley.edu/eswaran12/?
###############################################################
#####################################################
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
#dataMessy <- read.csv(file="kallisto_table_normalized_filtered.csv",header=TRUE)
dataMessy <- read.csv(file="kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
dataMessy <- dataMessy[,c("target_id","est_counts","sample")]
dataClean <- spread(dataMessy,key=sample,value=est_counts)
#rm(dataMessy)
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
dataClean <- ceiling(dataClean)
## remove tx without gene match
dataClean <- dataClean[!is.na(match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)),]
## remove all zero rows
dataClean <- dataClean[!rowSums(dataClean)==0,]
## remove genes with only one tx
geneTable <- table(as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]))
genesWithOneTx <- names(geneTable)[geneTable==1]
txFromGenesWithOneTx <- tx2gene$Ensembl.Transcript.ID[tx2gene$Ensembl.Gene.ID%in%genesWithOneTx]
dataClean <- dataClean[!rownames(dataClean)%in%as.character(txFromGenesWithOneTx),]

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

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
#save(dxd,file="dxdProstate.RData")
load("dxdProstate.RData")
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
names(padjGeneListShaffer) <- uniqueGenesStageII
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

##compare
sum(genesTx%in%genesSW)
sum(genesSW%in%genesTx) #all stage-wise genes are in tx-level analysis..

## genes/tx found in only one method
#genesOnlyStage <- genesSW[!genesSW%in%genesTx] #genes only found in a stage-wise analysis
#genesOnlyTx <- genesTx[!genesTx%in%genesSW]
#names(padjGeneListShaffer)=uniqueGenesStageII
#txOnlyTx <- significantTranscriptsTx[!significantTranscriptsTx%in%significantTranscriptsStageWise] #tx you only find in tx level analysis
#genesFromOnlyTxTx <- unique(dxr$groupID[dxr$featureID%in%txOnlyTx]) #genes that contain transcripts only found in a tx level analysis
#mean(genesFromOnlyTxTx%in%genesSW) #nearly all genes that contain transcripts that are only significant in a tx-level analysis are also found in a SW analysis.'

############################################################################
## genes with significant transcripts found only in the tx-level analysis###
############################################################################
#number of transcripts for the significant genes
tab1=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSW])))
#and for the genes with significant transcripts only found in a tx-level analysis
tab2=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesFromOnlyTxTx])))
#plot shows that genes with transcripts that are only significant in a tx-level analysis usually have more isoforms.
m=matrix(0,nrow=66,ncol=2)
rownames(m)=2:67
m[match(names(tab1),rownames(m)),1]=tab1
m[match(names(tab2),rownames(m)),2]=tab2
plus30=colSums(m[30:nrow(m),])
m2=rbind(m[1:29,],plus30)
rownames(m2)[30]="31+"

barplot(m2[,1], main="genes from transcripts that have only been found in a tx-level analysis", ylab="Frequency", xlab="Number of transcripts per gene", cex.lab=1.5, cex.axis=1.5)
barplot(m2[,2],add=TRUE,col=alpha("skyblue",.5), yaxt="n", xaxt="n")
#legend("topright",c("all genes found in stage-wise analysis","genes that contain transcripts only found in a transcript level analysis"), col=c("grey","skyblue"),lty=1,lwd=3, cex=1.3)

### combine DTU and DTE analyses
# intersection between DTE and DTU
genesSWDTU=genesSW
load("genesSWDTE.Rda")
genesSWDTE=genesSW
mean(genesSWDTU%in%genesSWDTE) #only 50%
genesSWDTUNotDTE=genesSWDTU[!genesSWDTU%in%genesSWDTE]
tab2=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSWDTUNotDTE])))
#make m from barplot
plot(x=1:nrow(m),y=m[,2]/m[,1]) #fraction of genes that are DTU but not DTE decreases for a higher number of transcripts.

## the more transcripts within a gene the higher the probability it will be found to be differentially used.
tab1=table(table(as.character(txGeneData$gene)))
tab2=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSWDTU])))
m=matrix(0,nrow=66,ncol=2)
rownames(m)=2:67
m[match(names(tab1),rownames(m)),1]=tab1
m[match(names(tab2),rownames(m)),2]=tab2
plus30=colSums(m[30:nrow(m),])
m2=rbind(m[1:29,],plus30)
rownames(m2)[30]="31+"
barplot(m2[,1], main="", ylab="Frequency", xlab="Number of transcripts per gene", cex.lab=1.5, cex.axis=1.5)
barplot(m2[,2],add=TRUE,col=alpha("skyblue",.5), yaxt="n", xaxt="n")
#isoform dominance is most studied mechanism. is this only the main DS mechanism for complex genes?
complexGeneID=names(table(as.character(txGeneData$gene)))[table(as.character(txGeneData$gene))>2] #genes with 3+ tx
table(unlist(lapply(padjGeneListShaffer[names(padjGeneListShaffer)%in%complexGeneID], function(x) sum(x<alphaAdjusted)))) #also many genes with 3+ tx differentially used
sum((table(unlist(lapply(padjGeneListShaffer[names(padjGeneListShaffer)%in%complexGeneID], function(x) sum(x<alphaAdjusted))))/sum(names(padjGeneListShaffer)%in%complexGeneID))[-c(1:3)]) #proportion of these genes with 3 or more significant transcripts


#gene with most transcripts DU (8)
select(Homo.sapiens, keys="ENSG00000187244", columns="SYMBOL", keytype="ENSEMBL")
# BCAM (also called CD239) plays a role in epithelial skin cancer
# it consists of 11 transcripts of which 8 are differentially used between normal and prostate cancer cell types

#gene with second most transcripts DU (7)
select(Homo.sapiens, keys="ENSG00000163110", columns="SYMBOL", keytype="ENSEMBL")
# PDLIM5 was recently associated with prostate cancer in a large-scale GWAS over multiple stages!

### plot the PDLIM5 gene
plotDTULog <- function(gene,transcripts,condition,nrGroups,nrPerGroup,interval,txInterval){
    #gene is gene you want to plot
    #transcripts are its transcripts
    #condition is the conditions
    #nrGroups is number of groups
    #nrPerGroup is number of replicates per group
    #interval is the total space for one transcript
    #txInterval is the space between transcript
    plotData <- dataClean[transcripts,]
    plotData <- sweep(plotData,2,FUN="/",STATS=colSums(plotData))
    plotData[is.na(plotData)] <- 0
    plotData <- plotData[,order(condition)] #grouped by factor of interest
    nTx=length(transcripts)
    #x-axis values
    x=sapply(1:nTx, function(txID) c((txID-interval)+txInterval*(0:(nrPerGroup[1]-1)), (txID+interval)-txInterval*(0:(nrPerGroup[2]-1))))
    #### plot counts of gene in main
    #plot(c(x),c(t(plotData)), col=rep(1:2,each=nrPerGroup[1]), pch=rep(1:19,each=nrPerGroup[1]), xaxt="n", ylab="Fraction used", xlab="", ylim=c(0,max(plotData)+.05), main=Reduce(paste,c(gene, "\n total counts=", as.numeric(colSums(dataClean[transcripts,])))))
    ###on log scale, dont plot counts (paper)
    plot(c(x),log(c(t(plotData+1e-3))), col=rep(1:2,each=nrPerGroup[1]), pch=rep(1:19,each=nrPerGroup[1]), xaxt="n", ylab="log(fraction used + .001)", xlab="", ylim=log(c(0,1)+1e-3), main=as.character(tx2gene$Associated.Gene.Name[match(gene,tx2gene$Ensembl.Gene.ID)]))
    abline(v=sapply(1:nTx, function(txID) c(txID-interval-txInterval, txID+interval+txInterval)), col=alpha("grey",.3), lty=2)
    axis(1,at=colMeans(x), labels=transcripts, cex.axis=.5, las=2)
    #put a S on tx found to be significant
    transcripts=paste0(gene,".",dxr$featureID[dxr$groupID==gene])
    sigTx <- which(transcripts%in%significantTranscriptsStageWise)
    points(x=sigTx,y=rep(max(log(plotData+1e-3))+.15,length(sigTx)), pch="S", cex=.9)
}

gene="ENSG00000163110"
transcripts=dxr$featureID[dxr$groupID==gene]
plotDTULog(gene=gene, transcripts=transcripts, condition=condition, nrGroups=2, nrPerGroup=c(14,14), interval=.14, txInterval=.01)









### JUNK
#### compare found DTU genes with genes found in paper
library(Homo.sapiens)
annotGenesSW <- select(Homo.sapiens, keys=genesSW, columns="SYMBOL", keytype="ENSEMBL")
geneNamesSW <- annotGenesSW$SYMBOL
geneNamesSW <- geneNamesSW[!is.na(geneNamesSW)]
library(openxlsx)
retainedIntronGenesPaper = read.xlsx(xlsxFile="/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/highlyReliableSplicingsPaper.xlsx", sheet=6, rows=2:80, cols=1)
mean(unique(retainedIntronGenesPaper[,1]) %in% geneNamesSW)
retainedIntronGenesPaper = read.xlsx(xlsxFile="/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/highlyReliableSplicingsPaper.xlsx", sheet=5, rows=2:80, cols=1)

genesSWDTU=genesSW
load("genesSWDTE.RData")
genesSWDTE=genesSW
rm(genesSW)
onlyDtuGenes = genesSWDTU[!genesSWDTU%in%genesSWDTE]
table(unlist(lapply(padjGeneListShaffer[onlyDtuGenes], function(x) sum(x<alphaAdjusted))))/sum(unlist(lapply(padjGeneListShaffer[onlyDtuGenes], function(x) sum(x<alphaAdjusted)))) #99.6% of genes that are DTU but not DTE show switch in dominance
library(Homo.sapiens) ; library(EGSEA) ; library(AnnotationHub)
entrezIDs <- select(Homo.sapiens, keys=onlyDtuGenes, columns="ENTREZID", keytype="ENSEMBL")
entrezIDs <- entrezIDs$ENTREZID
entrezIDs <- entrezIDs[!is.na(entrezIDs)]
index = buildMSigDBIdx(entrezIDs=entrezIDs, geneSets="all", species="Homo sapiens", min.size=1)
# EGSEA only works on voom...





