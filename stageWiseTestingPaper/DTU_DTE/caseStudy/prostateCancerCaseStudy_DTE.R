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

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

txGeneData = as.data.frame(cbind(rownames(dataClean),as.character(tx2gene$Ensembl.Transcript.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)]),as.character(tx2gene$Ensembl.Gene.ID[match(rownames(dataClean),tx2gene$Ensembl.Transcript.ID)])))
colnames(txGeneData)=c("tx","transcript","gene")
barplot(table(table(txGeneData$gene)), main="Distribution of number of tx per gene")

### DTE analysis
library(edgeR)
d=DGEList(dataClean)
d=calcNormFactors(d)
#keep=rowSums(cpm(d)>2)>=6
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

library(biomaRt)
topGenes=head(names(qvals)[order(qvals)],n=10)
mart<-useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
affy=c("202763_at","209310_s_at","207500_at")
select(mart, keys=topGenes, columns=c("entrezgene","hgnc_symbol"),keytype='ensembl_gene_id')


############################################################################
## genes with significant transcripts found only in the tx-level analysis###
############################################################################
## genes/tx found in only one method
genesOnlyStage <- genesSW[!genesSW%in%genesTx] #genes only found in a stage-wise analysis
genesOnlyTx <- genesTx[!genesTx%in%genesSW]
txOnlyTx <- significantTranscriptsTx[!significantTranscriptsTx%in%significantTranscriptsStageWise] #tx you only find in tx level analysis
genesFromOnlyTxTx <- unique(as.character(tx2gene$Ensembl.Gene.ID[match(txOnlyTx,tx2gene$Ensembl.Transcript.ID)])) #genes that contain transcripts only found in a tx level analysis

#number of transcripts for the significant genes
library(scales)
tab1=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSW])))
#and for the genes with significant transcripts only found in a tx-level analysis
tab2=table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesFromOnlyTxTx])))
#plot shows that genes with transcripts that are only significant in a tx-level analysis usually have more isoforms.
m=matrix(0,nrow=66,ncol=2)
rownames(m)=1:66
m[match(names(tab1),rownames(m)),1]=tab1
m[match(names(tab2),rownames(m)),2]=tab2
plus30=colSums(m[31:nrow(m),])
m2=rbind(m[1:30,],plus30)
rownames(m2)[31]="31+"

barplot(m2[,1], main="genes from transcripts that have only been found in a tx-level analysis", ylab="Frequency", xlab="Number of transcripts per gene", cex.lab=1.5, cex.axis=1.5)
barplot(m2[,2],add=TRUE,col=alpha("skyblue",.5), yaxt="n", xaxt="n")
legend("topright",c("all genes found in stage-wise analysis","genes that contain transcripts only found in a transcript level analysis"), col=c("grey","skyblue"),lty=1,lwd=3, cex=1.3)

##
tab1 = table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSW])))
genesSWNoTx = names(padjListShaffer)[unlist(lapply(padjListShaffer,function(x) sum(x<alphaAdjusted)))==0]
tab2 = table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSWNoTx])))
m=matrix(0,nrow=66,ncol=2)
rownames(m)=1:66
m[match(names(tab1),rownames(m)),1]=tab1
m[match(names(tab2),rownames(m)),2]=tab2
plus30=colSums(m[31:nrow(m),])
m2=rbind(m[1:30,],plus30)
rownames(m2)[31]="31+"
barplot(m2[,1], main="genes from stage-wise method", ylab="Frequency", xlab="Number of transcripts per gene", cex.lab=1.5, cex.axis=1.5)
barplot(m2[,2],add=TRUE,col=alpha("skyblue",.5), yaxt="n", xaxt="n")

tab1 = table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesTx])))
tab2 = table(table(as.character(txGeneData$gene[txGeneData$gene%in%genesSW])))


### DGE analysis
library(tximport)
geneCounts=rowsum(dataClean,group=txGeneData$gene)

countsMatTx <- txi$counts
lengthMatTx <- txi$length




library(dplyr)
dataGene = group_by(dataClean, by=txGeneData$gene)
geneCounts = summarise(dataGene, 
	  ERR031017=sum(ERR031017),
	  ERR031018=sum(ERR031018),
	  ERR031019=sum(ERR031019),
	  ERR299295=sum(ERR299295),
	  ERR299296=sum(ERR299296),
	  ERR031022=sum(ERR031022),
	  ERR031023=sum(ERR031023),
	  ERR031024=sum(ERR031024),
	  ERR031025=sum(ERR031025),
	  ERR031026=sum(ERR031026),
	  ERR031027=sum(ERR031027),
	  ERR031028=sum(ERR031028),
	  ERR031029=sum(ERR031029),
	  ERR031030=sum(ERR031030),
	  ERR031031=sum(ERR031031),
	  ERR031032=sum(ERR031032),
	  ERR031033=sum(ERR031033),
	  ERR299297=sum(ERR299297),
	  ERR031035=sum(ERR031035),
	  ERR299298=sum(ERR299298),
	  ERR299299=sum(ERR299299),
	  ERR031038=sum(ERR031038),
	  ERR031039=sum(ERR031039),
	  ERR031040=sum(ERR031040),
	  ERR031041=sum(ERR031041),
	  ERR031042=sum(ERR031042),
	  ERR031043=sum(ERR031043),
	  ERR031044=sum(ERR031044))
relativeAbundances <- mutate(dataGene,
			     ERR031017=ifelse(sum(ERR031017)==0,0,ERR031017/sum(ERR031017)),
			     ERR031018=ifelse(sum(ERR031018)==0,0,ERR031018/sum(ERR031018)))











#### trying closed testing procedure
tx=rownames(d$counts[which(genesAll=="ENSG00000002726"),])
combn(x=tx,m=2)
combn(x=tx,m=2,FUN=function(x) perGeneQValue_kvdb_closedTesting(pvals[rownames(lrt)%in%x]))



eval=list()
for(i in 1:30){
    print(i)
    tx=rownames(d$counts[which(genesAll==as.character(significantGenesQval[i])),])
    pvalId <- sapply(tx, function(x) which(rownames(lrt)==x))
    pvalLoop <- pvals[pvalId]
    names(pvalLoop) <- tx
    sigTx <- pvalLoop[which(pvalLoop<=alphaAdjusted)]
    if(length(sigTx)==0){
	eval[[i]]=FALSE
	next
    } else {
	eval[[i]] <- sapply(names(sigTx), function(transcript){
	combinationsAll <- sapply(1:length(tx), function(m) combn(x=tx,m=m))
	#get combinations for particular significant tx
	if(class(combinationsAll)=="list"){
	selectedCombinationsId <- lapply(combinationsAll, function(combo) apply(combo,2,function(column) any(column%in%transcript)))
	selectedCombinations <- sapply(1:length(combinationsAll), function(j) as.matrix(combinationsAll[[j]][,selectedCombinationsId[[j]]],nrow=j))
	qvalCombinations <- lapply(selectedCombinations, function(x) apply(x,2,function(combo) perGeneQValue_kvdb_closedTesting(pvalLoop[combo])))	
	} else if(class(combinationsAll)=="matrix"){
	    selectedCombinationsId <- apply(combinationsAll,2,function(column) any(column%in%transcript))
	    selectedCombinations <- combinationsAll[,selectedCombinationsId]
	    qvalCombinations <- apply(selectedCombinations,2,function(column) perGeneQValue_kvdb_closedTesting(pvalLoop[column]))
	}
	#calculate q-values for the respective combinations
	all(unlist(qvalCombinations)<=alphaAdjusted)
})

    }
}

sum(unlist(padjListShaffer[which(names(padjListShaffer)%in%significantGenesQval[1:30])])<=alphaAdjusted)
sum(unlist(eval[1:30]))





### help functions
perGeneQValue_kvdb_closedTesting <- function(pvals, method = perGeneQValueExact) {
  wTest= which(!is.na(pvals))
  geneID    = rep("gene1",length(wTest))
  geneSplit = split(seq(along=geneID), geneID)

  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)

  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  res        = q[match(pGene, theta)]
  res = pmin(1, res)
  names(res) = names(geneSplit)
  return(res)
}

perGeneQValueExact <- function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                            m = tab[notZero],
                            n = which(notZero))
  numerator    = sum(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)
}






