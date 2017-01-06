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

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

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
genesSW=rownames(getSignificantGenes(stageRObj))
significantTranscriptsStageWise=rownames(getSignificantTx(stageRObj))

## the more transcripts within a gene the higher the probability it will be found to be differentially used.
tab1=table(table(as.character(txGeneData$gene)))
tab2=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSW])))
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
padj=getAdjustedPValues(stageRObj)
#nr of complex genes with 3 or more significant transcripts
complexGeneWith3DUTx <- sum(sapply(complexGeneID, function(complexGene) sum(padj[padj[,"geneID"]%in%complexGene,"transcript"]<=0.05,na.rm=TRUE)>=3))
#proportion over all significant complex genes
complexGeneWith3DUTx/sum(complexGeneID%in%genesSW)

# genes with most transcripts DU
tail(sort(table(padj[which(padj[,"transcript"]<=0.05),"geneID"])))

library(Homo.sapiens)
#gene with most transcripts DU (8)
select(Homo.sapiens, keys="ENSG00000187244", columns="SYMBOL", keytype="ENSEMBL")
# BCAM (also called CD239) plays a role in epithelial skin cancer
# it consists of 11 transcripts of which 8 are differentially used between normal and prostate cancer cell types

#gene with second most transcripts DU (7)
select(Homo.sapiens, keys="ENSG00000163110", columns="SYMBOL", keytype="ENSEMBL")
# PDLIM5 was recently associated with prostate cancer in a large-scale GWAS over multiple stages.

#most significant genes
head(sort(qvalDxr))
select(Homo.sapiens, keys="ENSG00000106258", columns="SYMBOL", keytype="ENSEMBL")
# CYP3A5 was also previously associated with prostate cancer and is among the most significant genes

select(Homo.sapiens, keys="ENSG00000134324", columns="SYMBOL", keytype="ENSEMBL")
# LPIN1 was previously associated with prostate cancer where it was shown to be overexpressed in comparison to normal tissue: https://www.ncbi.nlm.nih.gov/pubmed/25834103

select(Homo.sapiens, keys="ENSG00000142973", columns="SYMBOL", keytype="ENSEMBL")
#CYP4B1 not associated with prostate cancer but with bladder and breast cancer

#Hence, 2 out of 5 highly significant genes are previously associated with prostate cancer but through differential expression analysis. Here we show that they are in fact differentially spliced.

### plot the genes
plotDTULog <- function(gene,transcripts,condition,nrGroups,nrPerGroup,interval,txInterval,pch, ...){
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
    #plot(c(x),log(c(t(plotData+1e-3))), col=rep(1:2,each=nrPerGroup[1]), pch=pch, xaxt="n", ylab="log(fraction used + .001)", xlab="", ylim=log(c(1e-3,1.15)), main=as.character(tx2gene$Associated.Gene.Name[match(gene,tx2gene$Ensembl.Gene.ID)]), ...)
    plot(c(x),log(c(t(plotData+1e-3))), col=rep(1:2,each=nrPerGroup[1]), pch=pch, xaxt="n", ylab="log(fraction used + .001)", xlab="", ylim=log(c(1e-3,1.15)), ...)
    
    abline(v=sapply(1:nTx, function(txID) c(txID-interval-txInterval, txID+interval+txInterval)), col=alpha("grey",.3), lty=2)
    axis(1,at=colMeans(x), labels=transcripts, cex.axis=.5, las=2)
    #put a S on tx found to be significant
    transcripts=paste0(gene,".",dxr$featureID[dxr$groupID==gene])
    sigTx <- which(transcripts%in%significantTranscriptsStageWise)
    points(x=sigTx,y=rep(max(log(plotData+1e-3))+.2,length(sigTx)), pch="S", cex=.9)
}

#gene="ENSG00000163110"
gene="ENSG00000106258"
#gene="ENSG00000134324"
transcripts=dxr$featureID[dxr$groupID==gene]
plotDTULog(gene=gene, transcripts=transcripts, condition=condition, nrGroups=2, nrPerGroup=c(14,14), interval=.14, txInterval=.01, pch=16,cex=.75)




plotDTU <- function(gene,transcripts,condition,nrGroups,nrPerGroup,interval,txInterval){
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
    plot(c(x),c(t(plotData)), col=rep(1:2,each=nrPerGroup[1]), pch=rep(1:19,each=nrPerGroup[1]), xaxt="n", ylab="Fraction used", xlab="", ylim=c(0,max(plotData)+.05), main=Reduce(paste,c(gene, "\n total counts=", as.numeric(colSums(dataClean[transcripts,])))))
    abline(v=sapply(1:nTx, function(txID) c(txID-interval-txInterval, txID+interval+txInterval)), col=alpha("grey",.3), lty=2)
    axis(1,at=colMeans(x), labels=transcripts, cex.axis=.5, las=2)
    #put a S on tx found to be significant
    sigTx <- which(transcripts%in%significantTranscriptsStageWise)
    points(x=sigTx,y=rep(max(plotData)+.03,length(sigTx)), pch="S", cex=.9)
}
## plot high abundant transcripts
gene="ENSG00000163110"
transcripts=dxr$featureID[dxr$groupID==gene]
transcriptHigh="ENST00000437932"
nrGroups=2
nrPerGroup=c(14,14)
interval=.14
txInterval=.01
plotData <- dataClean[transcripts,]
plotData <- sweep(plotData,2,FUN="/",STATS=colSums(plotData))
plotData[is.na(plotData)] <- 0
plotData <- plotData[,order(condition)] #grouped by factor of interest
plot(x=c(seq(0,0.5,length.out=14),seq(0.8,1.3,length.out=14)),y=plotData[transcriptHigh,],col=condition[order(condition)],xaxt="n",xlab="",ylab="Fraction used",pch=16,cex=.75)
axis(1,at=0.65, labels=transcriptHigh, cex.axis=.5, las=2)


plotDTU2 <- function(gene,transcripts,transcriptsToPlot,condition,nrGroups,nrPerGroup,interval,txInterval,...){
    #gene is gene you want to plot
    #transcripts are its transcripts
    #transcriptsToPlot are the transcripts you would want to plot
    #condition is the conditions
    #nrGroups is number of groups
    #nrPerGroup is number of replicates per group
    #interval is the total space for one transcript
    #txInterval is the space between transcript
    plotData <- dataClean[transcripts,]
    plotData <- sweep(plotData,2,FUN="/",STATS=colSums(plotData))
    plotData[is.na(plotData)] <- 0
    plotData <- plotData[,order(condition)] #grouped by factor of interest
    nTx=length(transcriptsToPlot)
    plotData=plotData[transcriptsToPlot,]
    #x-axis values
    x=sapply(1:nTx, function(txID) c((txID-interval)+txInterval*(0:(nrPerGroup[1]-1)), (txID+interval)-txInterval*(0:(nrPerGroup[2]-1))))
    #### plot counts of gene in main
    plot(c(x),c(t(plotData)), col=rep(1:2,each=nrPerGroup[1]), pch=16, xaxt="n", ylab="Fraction used", xlab="", ylim=c(0,max(plotData[transcriptsToPlot,])+.005),...)
    abline(v=sapply(1:nTx, function(txID) c(txID-interval-txInterval, txID+interval+txInterval)), col=alpha("grey",.3), lty=2)
    axis(1,at=colMeans(x), labels=transcriptsToPlot, cex.axis=.5, las=2)
    #put a S on tx found to be significant
    sigTx <- which(transcripts%in%significantTranscriptsStageWise)
    points(x=sigTx,y=rep(max(plotData)+.03,length(sigTx)), pch="S", cex=.9)
}

transcriptsToPlot=c("ENST00000359265","ENST00000380180","ENST00000508216","ENST00000511767","ENST00000627587")
plotDTU2(gene=gene, transcripts=transcripts, transcriptsToPlot=transcriptsToPlot, condition=condition, nrGroups=2, nrPerGroup=c(14,14), interval=.14, txInterval=.01,cex=.75)


layout(matrix(c(1,2,2,2),nrow=1))
par(mar=c(5,4.5,4,1))
gene="ENSG00000163110"
transcripts=dxr$featureID[dxr$groupID==gene]
transcripts=c("ENST00000437932",transcripts[!transcripts%in%"ENST00000437932"]) #dominant one
plot(x=c(seq(0,0.5,length.out=14),seq(0.8,1.3,length.out=14)),y=plotData[transcriptHigh,],col=condition[order(condition)],xaxt="n",xlab="",ylab="Fraction used",pch=16,cex=.75, cex.lab=1.5)
axis(1,at=0.65, labels=transcriptHigh, cex.axis=.5, las=2)
plotDTULog(gene=gene, transcripts=transcripts, condition=condition, nrGroups=2, nrPerGroup=c(14,14), interval=.14, txInterval=.01, pch=16,cex=.75, cex.lab=1.5)
rect(0.6, -3.2, 1.3, 0.1,col="darkgray",density=0)





