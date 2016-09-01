setwd("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/DGENew/")
source("http://130.60.190.4/robinson_lab/edgeR_robust/robust_simulation.R")

## pickrell simulation: original
library(tweeDEseqCountData)
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
nSamp <- 10 
nTags <- 10e3
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(15e6,20e6,length.out=nSamp)))
data <- NBsim(foldDiff = 3,  dataset = pickrell, nTags = nTags, pDiff=.05, group = grp, verbose = TRUE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE)
selectedMethods <- c("edgeR","limma_voom")
pvals <- pval(data, method=selectedMethods, count.type="counts", mc.cores=1)
getFDR(pvals) #inflated FDR of edgeR: 11%. voom stays at 5%.

source("robust_simulation_kvdb.R")
#in the code we just sourced, I made two changes to the simulation code:
#	- it allows you to input a vector of DE genes
#	- it allows you to input a seed for sampling the genes from the estimated dataset, which allows us to link genes across timepoints.
ind <- sample(1:nTags,(nTags/20),replace=FALSE) 
data <- NBsim(foldDiff = rep(3,length(ind)), ind=ind,  dataset = pickrell, nTags = 10e3, pDiff=.05, group = grp, verbose = TRUE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE, seed=1)
pvals <- pval(data, method=selectedMethods, count.type="counts", mc.cores=1)
getFDR(pvals) #same -> no changes in simulation framework
#also, cpm and dispersion is equal for same seed across simulations.


## simulating two timepoints. To avoid non-intended flipping in simulation framework, we set pUp=1.
# 500 genes only in t1
# 500 genes only in t2
# 500 genes DE in t1 and t2 with constant fold change
# 500 genes DE in t1 and t2 with reversed fold change

nTags <- 10e3
foldDiffT1 <- rep(c(3,3,3),each=500) #t1 unique, t1 cst, t1 interaction
foldDiffT2 <- rep(c(3,3,1/3),each=500) #t2 unique, t2 cst, t2 interaction

ind1 <- sample(10e3,1000) ## t1 t2 constant (first 500) and t1 t2 interaction genes (last 500)
remaining1 <- (1:10e3)[-ind1]
indT1Unique <- sample(remaining1,500)
remaining2 <- remaining1[!remaining1%in%indT1Unique]
indT2Unique <- sample(remaining2,500)
indT1 <- c(indT1Unique,ind1)
indT2 <- c(indT2Unique,ind1)

dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = TRUE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE, seed=1)
dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = TRUE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE, seed=1)
data <- cbind(dataT1$counts,dataT2$counts)

## for performance evaluations
indInteractionAll <- paste0("ids",c(ind1[501:1000],indT1Unique,indT2Unique))
indT1All <- paste0("ids",indT1)
indT2All <- paste0("ids",indT2)
indAll <- unique(c(indInteractionAll,indT1All,indT2All))
nonDeGenes <- paste0("ids",1:nTags)
nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

time <- factor(rep(c("T1","T2"),each=10))
treat <- factor(rep(rep(c("control","treatment"),each=5),2))
design <- model.matrix(~treat*time)
d <- DGEList(data)
d <- edgeR::calcNormFactors(d)
d <- estimateGLMCommonDisp(d,design)
d <- estimateGLMTrendedDisp(d,design)
d <- estimateGLMTagwiseDisp(d,design)
plotBCV(d)
plotMDS(d,labels=paste0(treat,time))
fit <- glmFit(d,design)
L <- matrix(0,nrow=4,ncol=3)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- c("t1","t2","interaction")
L[2,1]=1
L[c(2,4),2]=1
L[4,3]=1

## standard analysis
alpha=0.05
lrtList <- list()
for(i in 1:ncol(L)) lrtList[[i]] <- glmLRT(fit,contrast=L[,i])
foundT1 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
foundT2 <- rownames(d)[p.adjust(lrtList[[2]]$table$PValue,"BH")<alpha]
foundInt <- rownames(d)[p.adjust(lrtList[[3]]$table$PValue,"BH")<alpha]
fdrT1 <- mean(!foundT1%in%indT1All)
fdrT2 <- mean(!foundT2%in%indT2All)
fdrInt <- mean(!foundInt%in%indInteractionAll)
fpT1 <- foundT1[!foundT1%in%indT1All]
fpT2 <- foundT2[!foundT2%in%indT2All]
fpInt <- foundInt[!foundInt%in%indInteractionAll]
fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
overallFDR <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
powerInteraction <- mean(indInteractionAll%in%foundInt)
allGenesFound <- unique(c(foundT1,foundT2,foundInt))
nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

## stage-wise analysis
alpha=0.05
lrt1 <- glmLRT(fit,contrast=L)
significantGenesStageI <- rownames(d)[p.adjust(lrt1$table$PValue,"BH")<alpha]
fit2 <- fit[significantGenesStageI,]
alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(d)
lrtListSW <- list()
for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
foundT1SW <- rownames(fit2)[p.adjust(lrtListSW[[1]]$table$PValue,"BH")<alphaAdjusted]
foundT2SW <- rownames(fit2)[p.adjust(lrtListSW[[2]]$table$PValue,"BH")<alphaAdjusted]
foundIntSW <- rownames(fit2)[p.adjust(lrtListSW[[3]]$table$PValue,"BH")<alphaAdjusted]
fdrT1SW <- mean(!foundT1SW%in%indT1All)
fdrT2SW <- mean(!foundT2SW%in%indT2All)
fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageISW)))
powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)


#### looping it
N=10

fdrT1 <- vector(length=N)
fdrT2 <- vector(length=N)
fdrInt <- vector(length=N)
fdrAllHyp <- vector(length=N)
overallFDR <- vector(length=N)
powerInteraction <- vector(length=N)
nullGeneFDR <- vector(length=N)

fdrT1SW <- vector(length=N)
fdrT2SW <- vector(length=N)
fdrIntSW <- vector(length=N)
fdrAllHypSW <- vector(length=N)
overallFDRSW <- vector(length=N)
powerInteractionSW <- vector(length=N)
nullGeneFDRSW <- vector(length=N)

for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nSamp)))
	nTags <- 10e3
	
	foldDiffT1 <- rep(c(3,3,3),each=500) #t1 unique, t1 cst, t1 interaction
	foldDiffT2 <- rep(c(3,3,1/3),each=500) #t2 unique, t2 cst, t2 interaction

	ind1 <- sample(10e3,1000) ## t1 t2 constant (first 500) and t1 t2 interaction genes (last 500)
	remaining1 <- (1:10e3)[-ind1]
	indT1Unique <- sample(remaining1,500)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,500)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)

	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
data <- cbind(dataT1$counts,dataT2$counts)

	## for performance evaluations
	indInteractionAll <- paste0("ids",c(ind1[501:1000],indT1Unique,indT2Unique))
	indT1All <- paste0("ids",indT1)
	indT2All <- paste0("ids",indT2)
	indAll <- unique(c(indInteractionAll,indT1All,indT2All))
	nonDeGenes <- paste0("ids",1:nTags)
	nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

	time <- factor(rep(c("T1","T2"),each=10))
	treat <- factor(rep(rep(c("control","treatment"),each=5),2))
	design <- model.matrix(~treat*time)
	d <- DGEList(data)
	d <- edgeR::calcNormFactors(d)
	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)
	plotBCV(d)
	plotMDS(d,labels=paste0(treat,time))
	fit <- glmFit(d,design)
	L <- matrix(0,nrow=4,ncol=3)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1

	
	## standard analysis
	alpha=0.05
	lrtList <- list()
	for(i in 1:ncol(L)) lrtList[[i]] <- glmLRT(fit,contrast=L[,i])
	foundT1 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
	foundT2 <- rownames(d)[p.adjust(lrtList[[2]]$table$PValue,"BH")<alpha]
	foundInt <- rownames(d)[p.adjust(lrtList[[3]]$table$PValue,"BH")<alpha]
	fdrT1[iter] <- mean(!foundT1%in%indT1All)
	fdrT2[iter] <- mean(!foundT2%in%indT2All)
	fdrInt[iter] <- mean(!foundInt%in%indInteractionAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fdrAllHyp[iter] <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
	overallFDR[iter] <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
	powerInteraction[iter] <- mean(indInteractionAll%in%foundInt)
	allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR[iter] <- length(nullGenesFound)/length(allGenesFound)
	

	## stage-wise analysis
	alpha=0.05
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(d)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(d)
	lrtListSW <- list()
	for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
	foundT1SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundT2SW <- rownames(fit2)[lrtListSW[[2]]$table$PValue<alphaAdjusted]
	foundIntSW <- rownames(fit2)[lrtListSW[[3]]$table$PValue<alphaAdjusted]
	fdrT1SW[iter] <- mean(!foundT1SW%in%indT1All)
	fdrT2SW[iter] <- mean(!foundT2SW%in%indT2All)
	fdrIntSW[iter] <- mean(!foundIntSW%in%indInteractionAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fdrAllHypSW[iter] <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	overallFDRSW[iter] <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI)))
	powerInteractionSW[iter] <- mean(indInteractionAll%in%foundIntSW)
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW[iter] <- length(nullGenesFoundSW)/length(allGenesFoundSW)

}


boxplot(cbind(fdrAllHyp,fdrAllHypSW,
	      overallFDR,overallFDRSW,
	      nullGeneFDR,nullGeneFDRSW),
	boxwex=.2, at=c(.1,.3,.6,.8,1.1,1.3))
abline(h=0.05,col=2,lty=2,lwd=2)

boxplot(cbind(powerInteraction,powerInteractionSW))


### plugging in the numbers from Hammer dataset. There we found:
#	- 5073 genes DE in both t1 and t2 without interaction (5000)
#	- 1268 genes DE only in t1 (1000)
#	- 942 genes DE only in t2 (1000)
#	- 411 genes DE for all three contrasts (500)
#	As in Hammer dataset, I will simulate 13e3 genes


N=10
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
# Hammer numbers
#nCst <- 5000
#nT1 <- 1000
#nT2 <- 1000
#nInt <- 500
grp <- factor(rep(0:1,each=5))

fdrT1 <- vector(length=N)
fdrT2 <- vector(length=N)
fdrInt <- vector(length=N)
fdrAllHyp <- vector(length=N)
overallFDR <- vector(length=N)
powerInteraction <- vector(length=N)
nullGeneFDR <- vector(length=N)
powerT1 <- vector(length=N)
powerT2 <- vector(length=N)

fdrT1SW <- vector(length=N)
fdrT2SW <- vector(length=N)
fdrIntSW <- vector(length=N)
fdrAllHypSW <- vector(length=N)
overallFDRSW <- vector(length=N)
powerInteractionSW <- vector(length=N)
nullGeneFDRSW <- vector(length=N)
powerT1SW <- vector(length=N)
powerT2SW <- vector(length=N)
fdrScreenSW <- vector(length=N)

for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	#libSize = sample(round(seq(5e6,10e6,length.out=nreps*4))) #low libSize	
	nTags <- 13e3
	
	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 3),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	#foldDiffT1 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 2),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 1/2),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	foldDiffT1 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 0.8),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	foldDiffT2 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 1.3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))


	ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
	remaining1 <- (1:nTags)[-ind1]
	indT1Unique <- sample(remaining1,nT1)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,nT2)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)

	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[1:10], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[11:20], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	data <- cbind(dataT1$counts,dataT2$counts)

	## for performance evaluations
	indInteractionAll <- paste0("ids",c(ind1[(nCst+1):length(ind1)],indT1Unique,indT2Unique))
	indT1All <- paste0("ids",indT1)
	indT2All <- paste0("ids",indT2)
	indAll <- unique(c(indInteractionAll,indT1All,indT2All))
	nonDeGenes <- paste0("ids",1:nTags)
	nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

	time <- factor(rep(c("T1","T2"),each=nreps*2))
	treat <- factor(rep(rep(c("control","treatment"),each=nreps),2))
	design <- model.matrix(~treat*time)
	d <- DGEList(data)
	d <- edgeR::calcNormFactors(d)
	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)
	plotBCV(d)
	plotMDS(d,labels=paste0(treat,time))
	fit <- glmFit(d,design)
	L <- matrix(0,nrow=4,ncol=3)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1

	
	## standard analysis
	alpha=0.05
	lrtList <- list()
	for(i in 1:ncol(L)) lrtList[[i]] <- glmLRT(fit,contrast=L[,i])
	foundT1 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
	foundT2 <- rownames(d)[p.adjust(lrtList[[2]]$table$PValue,"BH")<alpha]
	foundInt <- rownames(d)[p.adjust(lrtList[[3]]$table$PValue,"BH")<alpha]
	fdrT1[iter] <- mean(!foundT1%in%indT1All)
	fdrT2[iter] <- mean(!foundT2%in%indT2All)
	fdrInt[iter] <- mean(!foundInt%in%indInteractionAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fdrAllHyp[iter] <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
	overallFDR[iter] <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR[iter] <- length(nullGenesFound)/length(allGenesFound)
	
	powerInteraction[iter] <- mean(indInteractionAll%in%foundInt)
	powerT1[iter] <- mean(indT1All%in%foundT1)
	powerT2[iter] <- mean(indT2All%in%foundT2)

	## stage-wise analysis
	alpha=0.05
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(fit)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	lrtListSW <- list()
	for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
	foundT1SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundT2SW <- rownames(fit2)[lrtListSW[[2]]$table$PValue<alphaAdjusted]
	foundIntSW <- rownames(fit2)[lrtListSW[[3]]$table$PValue<alphaAdjusted]
	fdrScreenSW[iter] <- mean(significantGenesStageI%in%nonDeGenes)	
	fdrT1SW[iter] <- mean(!foundT1SW%in%indT1All)
	fdrT2SW[iter] <- mean(!foundT2SW%in%indT2All)
	fdrIntSW[iter] <- mean(!foundIntSW%in%indInteractionAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fdrAllHypSW[iter] <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	overallFDRSW[iter] <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI)))
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW[iter] <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW[iter] <- mean(indInteractionAll%in%foundIntSW)
	powerT1SW[iter] <- mean(indT1All%in%foundT1SW)
	powerT2SW[iter] <- mean(indT2All%in%foundT2SW)


}

pdf("simulationFC2_highLibSize.pdf")
boxplot(cbind(fdrAllHyp,fdrAllHypSW,
	      overallFDR,overallFDRSW,
	      nullGeneFDR,nullGeneFDRSW),
	boxwex=.2, at=c(.1,.3,.6,.8,1.1,1.3),
	main=paste("nCst",nCst,"nT1",nT1,"nT2",nT2,"nInt",nInt,"FC=",2))
abline(h=0.05,col=2,lty=2,lwd=2)
boxplot(cbind(powerInteraction,powerInteractionSW), main=paste("nCst",nCst,"nT1",nT1,"nT2",nT2,"nInt",nInt,"FC=",2))
boxplot(cbind(fdrT1,fdrT2,fdrInt), main="Standard", ylab="FDR")
boxplot(cbind(fdrT1SW,fdrT2SW,fdrIntSW), main="Stage-wise", ylab="FDR")
dev.off()


### evaluate performance on multiple significance levels

N=30
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
# Hammer numbers
#nCst <- 5000
#nT1 <- 1000
#nT2 <- 1000
#nInt <- 500
grp <- factor(rep(0:1,each=5))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")



for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	#libSize = sample(round(seq(5e6,10e6,length.out=nreps*4))) #low libSize	
	nTags <- 13e3
	
	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 3),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	#foldDiffT1 <- unlist(c(mapply(rep,c(2,1/2, 2,1/2, 2),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt))))
	#foldDiffT2 <- unlist(c(mapply(rep,c(2,1/2, 2,1/2, 1/2),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt))))
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 0.8),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	#foldDiffT1 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 0.8),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 1.3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/1.75),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt))) #nice
	foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.75),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))





	ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
	remaining1 <- (1:nTags)[-ind1]
	indT1Unique <- sample(remaining1,nT1)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,nT2)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)

	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[1:10], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[11:20], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	data <- cbind(dataT1$counts,dataT2$counts)

	## for performance evaluations
	indInteractionAll <- paste0("ids",c(ind1[(nCst+1):length(ind1)],indT1Unique,indT2Unique))
	indT1All <- paste0("ids",indT1)
	indT2All <- paste0("ids",indT2)
	indAll <- unique(c(indInteractionAll,indT1All,indT2All))
	nonDeGenes <- paste0("ids",1:nTags)
	nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

	time <- factor(rep(c("T1","T2"),each=nreps*2))
	treat <- factor(rep(rep(c("control","treatment"),each=nreps),2))
	design <- model.matrix(~treat*time)

	##edgeR prep.
	d <- DGEList(data)
	d <- edgeR::calcNormFactors(d)
	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)
	#plotBCV(d)
	#plotMDS(d,labels=paste0(treat,time))
	fit <- glmFit(d,design)
	L <- matrix(0,nrow=4,ncol=3)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1
	
	##limma prep.
	v <- voom(d,design, plot=FALSE)
	fitLimma <- lmFit(v,design)
	colnames(design)[4] <- "treatTime2"
	contrast.matrix <- makeContrasts(t1=treattreatment, 
					 t2=treattreatment+treatTime2,
					 interaction=treatTime2, levels=design)
	fitLimma <- contrasts.fit(fitLimma,contrast.matrix)
	fitLimma <- eBayes(fitLimma)
	
	## edgeR analysis
	## standard analysis
	resultsMatRegular01[iter,] <- doRegularAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatRegular05[iter,] <- doRegularAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatRegular10[iter,] <- doRegularAnalysis(data,L,alpha=0.10,design,fit)	
	## stage-wise analysis
	resultsMatSW01[iter,] <- doStageWiseAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatSW05[iter,] <- doStageWiseAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatSW10[iter,] <- doStageWiseAnalysis(data,L,alpha=0.10,design,fit)

	## limma analysis
	## standard analysis
	resultsMatRegular01Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)	
	## stage-wise analysis
	resultsMatSW01Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)

	
}

pdf("simulationMultLevels_niceFC_limmaEdgeR.pdf")
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="False Discovery Rate", main=paste("edgeR","nCst",nCst,"nT1",nT1,"nT2",nT2,"nInt",nInt,"FC 2"))
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",.8))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2)
text(x=0.25,y=0.1,"All Hypotheses")
text(x=3,y=0.1,"Gene level")
text(x=5.6,y=0.1,"Null gene level")


boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery rate", main="", xlab="False discovery rate cut-off")
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",.8))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2)
text(x=0.25,y=0.12,"all hypotheses")
text(x=3,y=0.12,"OFDR")
text(x=5.6,y=0.12,"null gene OFDR")
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=.8, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerInteraction"],resultsMatSW01[,"powerInteractionSW"],
	      resultsMatRegular05[,"powerInteraction"],resultsMatSW05[,"powerInteractionSW"],
	      resultsMatRegular10[,"powerInteraction"],resultsMatSW10[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power Interaction effect", main="edgeR")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=.8, lwd=2)


boxplot(cbind(resultsMatRegular01Limma[,"powerInteraction"],resultsMatSW01Limma[,"powerInteractionSW"],
	      resultsMatRegular05Limma[,"powerInteraction"],resultsMatSW05Limma[,"powerInteractionSW"],
	      resultsMatRegular10Limma[,"powerInteraction"],resultsMatSW10Limma[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power Interaction effect", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=.8, lwd=2)
mean(resultsMatSW01Limma[,"powerInteractionSW"]-resultsMatRegular01Limma[,"powerInteraction"])
#[1] 0.06504444
mean(resultsMatSW05Limma[,"powerInteractionSW"]-resultsMatRegular05Limma[,"powerInteraction"])
#[1] 0.04766667
mean(resultsMatSW10Limma[,"powerInteractionSW"]-resultsMatRegular10Limma[,"powerInteraction"])
#[1] 0.04271111


## OFDR plot introduction
#par(bty="l")
boxplot(resultsMatRegular05Limma[,"overallFDR"],ylab="OFDR", yaxt="n", bty="l")
axis(2,at=c(0.05,0.055,0.06))
abline(h=0.05,col="red",lty=2,lwd=2)


boxplot(cbind(resultsMatRegular05[,"fdrT1"], resultsMatRegular05[,"fdrT2"], resultsMatRegular05[,"fdrInteraction"]), labels=c("t1","t2","interaction"), ylab="FDR")
boxplot(cbind(resultsMatSW05[,"fdrT1SW"], resultsMatSW05[,"fdrT2SW"], resultsMatSW05[,"fdrInteractionSW"]), labels=c("t1","t2","interaction"), ylab="FDR")

## main effects contrasts
#t1
boxplot(cbind(resultsMatRegular01Limma[,"powerT1"],resultsMatSW01Limma[,"powerT1SW"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT1SW"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT1SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=.8, lwd=2)

#t2
boxplot(cbind(resultsMatRegular01Limma[,"powerT2"],resultsMatSW01Limma[,"powerT2SW"],
	      resultsMatRegular05Limma[,"powerT2"],resultsMatSW05Limma[,"powerT2SW"],
	      resultsMatRegular10Limma[,"powerT2"],resultsMatSW10Limma[,"powerT2SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=.8, lwd=2)


dev.off()





##### sampling interaction FC from beta
returnshapeParamKvdb <- function(nDE, min, max, mean, p95){
    gamma <- (max-mean)/(mean-min)
    targetfunction <- function(x){
      pbeta((p95-min)/(max-min),x,gamma*x)-0.95
    }
    alpha <- uniroot(targetfunction,c(0,10^6))$root
    beta <- alpha * gamma
    x <- rbeta(nDE,alpha,beta)
    y <- min+(max-min)*x
    foldSeq <- 2^y
    return(foldSeq)
  }
N=10
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
# Hammer numbers
#nCst <- 5000
#nT1 <- 1000
#nT2 <- 1000
#nInt <- 500
grp <- factor(rep(0:1,each=5))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")



for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	#libSize = sample(round(seq(5e6,10e6,length.out=nreps*4))) #low libSize	
	nTags <- 13e3
	
	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 3),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	#foldDiffT1 <- unlist(c(mapply(rep,c(2,1/2, 2,1/2, 2),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt))))
	#foldDiffT2 <- unlist(c(mapply(rep,c(2,1/2, 2,1/2, 1/2),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt))))
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 0.8),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	#foldDiffT1 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 0.8),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(2,1/2, 2,1/2, 1.3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3),c(nT1/2,nT1/2, nCst/2,nCst/2))) #nice
	foldDiffT1 <- c(foldDiffT1,1/returnshapeParamKvdb(nDE=nInt,min=log2(.001),max=1,mean=.2,p95=.8))
	foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3),c(nT2/2,nT2/2, nCst/2,nCst/2)))
	foldDiffT2 <- c(foldDiffT2,returnshapeParamKvdb(nDE=nInt,min=log2(.001),max=1,mean=.2,p95=.8))	




	ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
	remaining1 <- (1:nTags)[-ind1]
	indT1Unique <- sample(remaining1,nT1)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,nT2)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)

	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[1:10], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[11:20], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	data <- cbind(dataT1$counts,dataT2$counts)

	## for performance evaluations
	indInteractionAll <- paste0("ids",c(ind1[(nCst+1):length(ind1)],indT1Unique,indT2Unique))
	indT1All <- paste0("ids",indT1)
	indT2All <- paste0("ids",indT2)
	indAll <- unique(c(indInteractionAll,indT1All,indT2All))
	nonDeGenes <- paste0("ids",1:nTags)
	nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

	time <- factor(rep(c("T1","T2"),each=nreps*2))
	treat <- factor(rep(rep(c("control","treatment"),each=nreps),2))
	design <- model.matrix(~treat*time)

	##edgeR prep.
	d <- DGEList(data)
	d <- edgeR::calcNormFactors(d)
	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)
	#plotBCV(d)
	#plotMDS(d,labels=paste0(treat,time))
	fit <- glmFit(d,design)
	L <- matrix(0,nrow=4,ncol=3)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1
	
	##limma prep.
	v <- voom(d,design, plot=FALSE)
	fitLimma <- lmFit(v,design)
	colnames(design)[4] <- "treatTime2"
	contrast.matrix <- makeContrasts(t1=treattreatment, 
					 t2=treattreatment+treatTime2,
					 interaction=treatTime2, levels=design)
	fitLimma <- contrasts.fit(fitLimma,contrast.matrix)
	fitLimma <- eBayes(fitLimma)
	
	## edgeR analysis
	## standard analysis
	resultsMatRegular01[iter,] <- doRegularAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatRegular05[iter,] <- doRegularAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatRegular10[iter,] <- doRegularAnalysis(data,L,alpha=0.10,design,fit)	
	## stage-wise analysis
	resultsMatSW01[iter,] <- doStageWiseAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatSW05[iter,] <- doStageWiseAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatSW10[iter,] <- doStageWiseAnalysis(data,L,alpha=0.10,design,fit)

	## limma analysis
	## standard analysis
	resultsMatRegular01Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)	
	## stage-wise analysis
	resultsMatSW01Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
	
}

# Indeed this simulation gives a higher power difference for interaction effects, but now the standard analysis also controls the OFDR..





















### simulate yourself: code below gives same output as using the entire simulation framework.
nTags=13e3
nreps=5
nlibs=4*nreps
lib.size=libSize

datasetPickrell <- getDataset(exprs(pickrell.eset), drop.extreme.dispersion=0.1, drop.low.lambda=TRUE)
AveLogCPM <- datasetPickrell$dataset.AveLogCPM
dispersion <- datasetPickrell$dataset.dispersion
set.seed(1)
id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
Lambda <- 2^(AveLogCPM[id_r])
Lambda <- Lambda/sum(Lambda)
Dispersion <- dispersion[id_r]
Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
## add fold changes
foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 3),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/3),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
remaining1 <- (1:nTags)[-ind1]
indT1Unique <- sample(remaining1,nT1)
remaining2 <- remaining1[!remaining1%in%indT1Unique]
indT2Unique <- sample(remaining2,nT2)
indT1 <- c(indT1Unique,ind1)
indT2 <- c(indT2Unique,ind1)
Lambda[indT1,6:10] <- Lambda[indT1,6:10]*foldDiffT1
Lambda[indT2,16:20] <- Lambda[indT2,16:20]*foldDiffT2
## simulate
counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs) 
rownames(counts) <- paste("ids", 1:nTags, sep = "")




##### simulate yourself and customize further:
#	- add a number of genes with main time effects
#	- sample fold changes from Hammer dataset
nTags=13e3
nreps=5
nlibs=4*nreps
lib.size=sample(round(seq(5e6,10e6,length.out=nreps*4)))
nT1 <- 1000
nT2 <- 1000
nCst <- 2000
nInt <- 1000
nTime <- 500
load("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/contrastsSim.RData")
resultsShaffer <- contrastsSim[[1]]
foldChangesHammer <- contrastsSim[[2]]
datasetPickrell <- getDataset(exprs(pickrell.eset), drop.extreme.dispersion=0.1, drop.low.lambda=TRUE)


## sample parameters
AveLogCPM <- datasetPickrell$dataset.AveLogCPM
dispersion <- datasetPickrell$dataset.dispersion
set.seed(iter)
id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
Lambda <- 2^(AveLogCPM[id_r])
Lambda <- Lambda/sum(Lambda)
Dispersion <- dispersion[id_r]
Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
## sample fold changes
fcT1 <- sample(2^foldChangesHammer[resultsShaffer[,1]!=0,1],size=nT1,replace=TRUE) 
fcT2 <- sample(2^foldChangesHammer[resultsShaffer[,2]!=0,2],size=nT2,replace=TRUE)
fcCst <- sample(2^c(foldChangesHammer[resultsShaffer[,1]!=0 & resultsShaffer[,2]!=0 & resultsShaffer[,3]==0,1:2]),size=nCst,replace=TRUE)
fcInt <- 2^foldChangesHammer[resultsShaffer[,3]!=0,1:2]
fcIntT1T2 <- fcInt[sample(nrow(fcInt),size=nInt,replace=TRUE),]
fcIntT1 <- fcIntT1T2[,1]
fcIntT2 <- fcIntT1T2[,2]
foldDiffT1 <- c(fcT1,fcCst,fcIntT1)
foldDiffT2 <- c(fcT2,fcCst,fcIntT2)
## sample DE genes
indCst <- sample(nTags,nCst)
remaining <- (1:nTags)[-indCst]
indInt <- sample(remaining,nInt)
ind1 <- c(indCst,indInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
remaining1 <- remaining[!remaining%in%ind1]
indT1Unique <- sample(remaining1,nT1)
remaining2 <- remaining1[!remaining1%in%indT1Unique]
indT2Unique <- sample(remaining2,nT2)
indT1 <- c(indT1Unique,ind1) # t1 unique, t1t2 constant, t1 interaction
indT2 <- c(indT2Unique,ind1)
## up and downregulate simultaneously
fcDirT1 <- ifelse(foldDiffT1>1,1,-1)
fcDirT2 <- ifelse(foldDiffT2>1,1,-1)
Lambda[indT1,6:10] <- Lambda[indT1,6:10]*exp(log(foldDiffT1)/2*fcDirT1)
Lambda[indT1,1:5] <- Lambda[indT1,1:5]*exp(log(foldDiffT1)/2*(-fcDirT1))
Lambda[indT2,16:20] <- Lambda[indT2,16:20]*exp(log(foldDiffT2)/2*fcDirT2)
Lambda[indT2,11:15] <- Lambda[indT2,11:15]*exp(log(foldDiffT2)/2*(-fcDirT2))
#Lambda[indT1,6:10] <- Lambda[indT1,6:10]*foldDiffT1
#Lambda[indT2,16:20] <- Lambda[indT2,16:20]*foldDiffT2
## simulate counts
counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs) 
rownames(counts) <- paste("ids", 1:nTags, sep = "")
data=counts



##### looping it
N=10
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
nTime <- 500
grp <- factor(rep(0:1,each=5))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=10)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW")




for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4)))
	#lib.size = sample(round(seq(5e6,10e6,length.out=nreps*4)))	
	nTags <- 13e3
	
## sample parameters
AveLogCPM <- datasetPickrell$dataset.AveLogCPM
dispersion <- datasetPickrell$dataset.dispersion
set.seed(iter)
id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
Lambda <- 2^(AveLogCPM[id_r])
Lambda <- Lambda/sum(Lambda)
Dispersion <- dispersion[id_r]
Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
## sample fold changes
fcT1 <- sample(2^foldChangesHammer[resultsShaffer[,1]!=0,1],size=nT1,replace=TRUE) 
fcT2 <- sample(2^foldChangesHammer[resultsShaffer[,2]!=0,2],size=nT2,replace=TRUE)
fcCst <- sample(2^c(foldChangesHammer[resultsShaffer[,1]!=0 & resultsShaffer[,2]!=0 & resultsShaffer[,3]==0,1:2]),size=nCst,replace=TRUE)
fcInt <- 2^foldChangesHammer[resultsShaffer[,3]!=0,1:2]
fcIntT1T2 <- fcInt[sample(nrow(fcInt),size=nInt,replace=TRUE),]
fcIntT1 <- fcIntT1T2[,1]
fcIntT2 <- fcIntT1T2[,2]
foldDiffT1 <- c(fcT1,fcCst,fcIntT1)
foldDiffT2 <- c(fcT2,fcCst,fcIntT2)
## sample DE genes
indCst <- sample(nTags,nCst)
remaining <- (1:nTags)[-indCst]
indInt <- sample(remaining,nInt)
ind1 <- c(indCst,indInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
remaining1 <- remaining[!remaining%in%ind1]
indT1Unique <- sample(remaining1,nT1)
remaining2 <- remaining1[!remaining1%in%indT1Unique]
indT2Unique <- sample(remaining2,nT2)
indT1 <- c(indT1Unique,ind1) # t1 unique, t1t2 constant, t1 interaction
indT2 <- c(indT2Unique,ind1)
## up and downregulate simultaneously
fcDirT1 <- ifelse(foldDiffT1>1,1,-1)
fcDirT2 <- ifelse(foldDiffT2>1,1,-1)
Lambda[indT1,6:10] <- Lambda[indT1,6:10]*exp(log(foldDiffT1)/2*fcDirT1)
Lambda[indT1,1:5] <- Lambda[indT1,1:5]*exp(log(foldDiffT1)/2*(-fcDirT1))
Lambda[indT2,16:20] <- Lambda[indT2,16:20]*exp(log(foldDiffT2)/2*fcDirT2)
Lambda[indT2,11:15] <- Lambda[indT2,11:15]*exp(log(foldDiffT2)/2*(-fcDirT2))
#Lambda[indT1,6:10] <- Lambda[indT1,6:10]*foldDiffT1
#Lambda[indT2,16:20] <- Lambda[indT2,16:20]*foldDiffT2
## simulate counts
counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs) 
rownames(counts) <- paste("ids", 1:nTags, sep = "")
data=counts

	## for performance evaluations
	indInteractionAll <- paste0("ids",c(ind1[(nCst+1):length(ind1)],indT1Unique,indT2Unique))
	indT1All <- paste0("ids",indT1)
	indT2All <- paste0("ids",indT2)
	indAll <- unique(c(indInteractionAll,indT1All,indT2All))
	nonDeGenes <- paste0("ids",1:nTags)
	nonDeGenes <- nonDeGenes[!nonDeGenes%in%indAll]

	time <- factor(rep(c("T1","T2"),each=nreps*2))
	treat <- factor(rep(rep(c("control","treatment"),each=nreps),2))
	design <- model.matrix(~treat*time)

	##edgeR prep.
	d <- DGEList(data)
	d <- edgeR::calcNormFactors(d)
	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)
	#plotBCV(d)
	#plotMDS(d,labels=paste0(treat,time))
	fit <- glmFit(d,design)
	L <- matrix(0,nrow=4,ncol=3)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1
	
	##limma prep.
	v <- voom(d,design, plot=FALSE)
	fitLimma <- lmFit(v,design)
	colnames(design)[4] <- "treatTime2"
	contrast.matrix <- makeContrasts(t1=treattreatment, 
					 t2=treattreatment+treatTime2,
					 interaction=treatTime2, levels=design)
	fitLimma <- contrasts.fit(fitLimma,contrast.matrix)
	fitLimma <- eBayes(fitLimma)
	
	## edgeR analysis
	## standard analysis
	resultsMatRegular01[iter,] <- doRegularAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatRegular05[iter,] <- doRegularAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatRegular10[iter,] <- doRegularAnalysis(data,L,alpha=0.10,design,fit)	
	## stage-wise analysis
	resultsMatSW01[iter,] <- doStageWiseAnalysis(data,L,alpha=0.01,design,fit)
	resultsMatSW05[iter,] <- doStageWiseAnalysis(data,L,alpha=0.05,design,fit)
	resultsMatSW10[iter,] <- doStageWiseAnalysis(data,L,alpha=0.10,design,fit)

	## limma analysis
	## standard analysis
	resultsMatRegular01Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)	
	## stage-wise analysis
	resultsMatSW01Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
}


boxplot(cbind(fdrAllHyp,fdrAllHypSW,
	      overallFDR,overallFDRSW,
	      nullGeneFDR,nullGeneFDRSW),
	boxwex=.2, at=c(.1,.3,.6,.8,1.1,1.3),
	main=paste("nCst",nCst,"nT1",nT1,"nT2",nT2,"nInt",nInt,"FC=",2))
abline(h=0.05,col=2,lty=2,lwd=2)
boxplot(cbind(powerInteraction,powerInteractionSW), main=paste("nCst",nCst,"nT1",nT1,"nT2",nT2,"nInt",nInt,"FC=",2))
boxplot(cbind(fdrT1,fdrT2,fdrInt), main="Standard", ylab="FDR")
boxplot(cbind(fdrT1SW,fdrT2SW,fdrIntSW), main="Stage-wise", ylab="FDR")



