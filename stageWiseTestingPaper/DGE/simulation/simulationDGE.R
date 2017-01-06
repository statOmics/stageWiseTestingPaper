setwd("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper_public/stageWiseTestingPaper/DGE/simulation")
source("http://130.60.190.4/robinson_lab/edgeR_robust/robust_simulation.R")
library(scales)
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


### evaluate performance on multiple significance levels
source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper_public/stageWiseTestingPaper/DGE/simulation/simulationDGE_helpFunctions.R")
N=30
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
grp <- factor(rep(0:1,each=5))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")

fdrExtraGenesRegular01 <- vector(length=N)
fdrExtraGenesRegular05 <- vector(length=N)
fdrExtraGenesRegular10 <- vector(length=N)

fdrExtraGenesRegular01Limma <- vector(length=N)
fdrExtraGenesRegular05Limma <- vector(length=N)
fdrExtraGenesRegular10Limma <- vector(length=N)


for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	nTags <- 13e3
	
	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/1.75),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.75),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))

	#indicators
	ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
	remaining1 <- (1:nTags)[-ind1]
	indT1Unique <- sample(remaining1,nT1)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,nT2)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)
	
	#simulate
	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[1:10], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[11:20], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	data <- cbind(dataT1$counts,dataT2$counts)

	##aggregate for performance evaluations
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

	## compare extra genes found in standard to stage-wise analysis
	fdrExtraGenesRegular01[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.01,design,fit=fit)
	fdrExtraGenesRegular05[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.05,design,fit=fit)
	fdrExtraGenesRegular10[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.1,design,fit=fit)
	
	fdrExtraGenesRegular01Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.01,design,fit=fitLimma)
	fdrExtraGenesRegular05Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.05,design,fit=fitLimma)
	fdrExtraGenesRegular10Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.1,design,fit=fitLimma)	
	
}

## OFDR and power interaction plot for limma
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
text(x=0.25,y=0.12,"all hypotheses", cex=2)
text(x=2.8,y=0.12,"OFDR", cex=2)
text(x=5.4,y=0.12,"null genes", cex=2)
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerInteraction"],resultsMatSW01Limma[,"powerInteractionSW"],
	      resultsMatRegular05Limma[,"powerInteraction"],resultsMatSW05Limma[,"powerInteractionSW"],
	      resultsMatRegular10Limma[,"powerInteraction"],resultsMatSW10Limma[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)


### OFDR and power interaction plot for edgeR
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off", ylim=c(0.005,0.13))
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
text(x=0.25,y=0.12,"all hypotheses", cex=2)
text(x=2.8,y=0.12,"OFDR", cex=2)
text(x=5.4,y=0.12,"null genes", cex=2)
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerInteraction"],resultsMatSW01[,"powerInteractionSW"],
	      resultsMatRegular05[,"powerInteraction"],resultsMatSW05[,"powerInteractionSW"],
	      resultsMatRegular10[,"powerInteraction"],resultsMatSW10[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)


### t1 and t2 power limma
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01Limma[,"powerT1"],resultsMatSW01Limma[,"powerT1SW"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT1SW"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT1SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerT2"],resultsMatSW01Limma[,"powerT2SW"],
	      resultsMatRegular05Limma[,"powerT2"],resultsMatSW05Limma[,"powerT2SW"],
	      resultsMatRegular10Limma[,"powerT2"],resultsMatSW10Limma[,"powerT2SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)



### t1 and t2 power edgeR
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01[,"powerT1"],resultsMatSW01[,"powerT1SW"],
	      resultsMatRegular05[,"powerT1"],resultsMatSW05[,"powerT1SW"],
	      resultsMatRegular10[,"powerT1"],resultsMatSW10[,"powerT1SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerT2"],resultsMatSW01[,"powerT2SW"],
	      resultsMatRegular05[,"powerT2"],resultsMatSW05[,"powerT2SW"],
	      resultsMatRegular10[,"powerT2"],resultsMatSW10[,"powerT2SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

### nr of genes found and FDR of extra genes for limma
boxplot(cbind(resultsMatRegular01Limma[,"totalGenesFound"],resultsMatSW01Limma[,"totalGenesFoundSW"],
	      resultsMatRegular05Limma[,"totalGenesFound"],resultsMatSW05Limma[,"totalGenesFoundSW"],
	      resultsMatRegular10Limma[,"totalGenesFound"],resultsMatSW10Limma[,"totalGenesFoundSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(fdrExtraGenesRegular01Limma,
	      fdrExtraGenesRegular05Limma,
	      fdrExtraGenesRegular10Limma),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))


### nr of genes found and FDR of extra genes for edgeR
boxplot(cbind(resultsMatRegular01[,"totalGenesFound"],resultsMatSW01[,"totalGenesFoundSW"],
	      resultsMatRegular05[,"totalGenesFound"],resultsMatSW05[,"totalGenesFoundSW"],
	      resultsMatRegular10[,"totalGenesFound"],resultsMatSW10[,"totalGenesFoundSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(fdrExtraGenesRegular01,
	      fdrExtraGenesRegular05,
	      fdrExtraGenesRegular10),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))



### same simulation but with three replicates.
N=30
nreps <- 3
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
grp <- factor(rep(0:1,each=3))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW")

fdrExtraGenesRegular01 <- vector(length=N)
fdrExtraGenesRegular05 <- vector(length=N)
fdrExtraGenesRegular10 <- vector(length=N)

fdrExtraGenesRegular01Limma <- vector(length=N)
fdrExtraGenesRegular05Limma <- vector(length=N)
fdrExtraGenesRegular10Limma <- vector(length=N)


for(iter in 1:N){
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	nTags <- 13e3
	
	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/1.75),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.75),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))

	#indicators
	ind1 <- sample(nTags,nCst+nInt) ## t1 t2 constant (first elements) and t1 t2 interaction genes (last elements)
	remaining1 <- (1:nTags)[-ind1]
	indT1Unique <- sample(remaining1,nT1)
	remaining2 <- remaining1[!remaining1%in%indT1Unique]
	indT2Unique <- sample(remaining2,nT2)
	indT1 <- c(indT1Unique,ind1)
	indT2 <- c(indT2Unique,ind1)
	
	#simulate
	dataT1 <- NBsim(foldDiff = foldDiffT1, ind=indT1,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[1:(nreps*2)], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	dataT2 <- NBsim(foldDiff = foldDiffT2, ind=indT2,  dataset = pickrell, nTags = nTags, pDiff=.05, pUp=1, group = grp, verbose = FALSE,  drop.extreme.dispersion = 0.1, lib.size=libSize[(nreps*2+1):(nreps*4)], drop.low.lambda=TRUE, add.outlier=FALSE, seed=iter)
	data <- cbind(dataT1$counts,dataT2$counts)

	##aggregate for performance evaluations
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

	## compare extra genes found in standard to stage-wise analysis
	fdrExtraGenesRegular01[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.01,design,fit=fit)
	fdrExtraGenesRegular05[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.05,design,fit=fit)
	fdrExtraGenesRegular10[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.1,design,fit=fit)
	
	fdrExtraGenesRegular01Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.01,design,fit=fitLimma)
	fdrExtraGenesRegular05Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.05,design,fit=fitLimma)
	fdrExtraGenesRegular10Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.1,design,fit=fitLimma)	

}

## OFDR and power interaction plot for limma
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
text(x=0.25,y=0.105,"all hypotheses", cex=2)
text(x=2.8,y=0.105,"OFDR", cex=2)
text(x=5.4,y=0.105,"null genes", cex=2)
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerInteraction"],resultsMatSW01Limma[,"powerInteractionSW"],
	      resultsMatRegular05Limma[,"powerInteraction"],resultsMatSW05Limma[,"powerInteractionSW"],
	      resultsMatRegular10Limma[,"powerInteraction"],resultsMatSW10Limma[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)

### OFDR and power interaction plot for edgeR
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off", ylim=c(0.005,0.13))
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=2, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=2, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=2, lwd=2)
text(x=0.25,y=0.12,"all hypotheses", cex=2)
text(x=2.8,y=0.12,"OFDR", cex=2)
text(x=5.4,y=0.12,"null genes", cex=2)
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerInteraction"],resultsMatSW01[,"powerInteractionSW"],
	      resultsMatRegular05[,"powerInteraction"],resultsMatSW05[,"powerInteractionSW"],
	      resultsMatRegular10[,"powerInteraction"],resultsMatSW10[,"powerInteractionSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2.5, lwd=2)


### t1 and t2 power limma
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01Limma[,"powerT1"],resultsMatSW01Limma[,"powerT1SW"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT1SW"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT1SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerT2"],resultsMatSW01Limma[,"powerT2SW"],
	      resultsMatRegular05Limma[,"powerT2"],resultsMatSW05Limma[,"powerT2SW"],
	      resultsMatRegular10Limma[,"powerT2"],resultsMatSW10Limma[,"powerT2SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)



### t1 and t2 power edgeR
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01[,"powerT1"],resultsMatSW01[,"powerT1SW"],
	      resultsMatRegular05[,"powerT1"],resultsMatSW05[,"powerT1SW"],
	      resultsMatRegular10[,"powerT1"],resultsMatSW10[,"powerT1SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerT2"],resultsMatSW01[,"powerT2SW"],
	      resultsMatRegular05[,"powerT2"],resultsMatSW05[,"powerT2SW"],
	      resultsMatRegular10[,"powerT2"],resultsMatSW10[,"powerT2SW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

### nr of genes found and FDR of extra genes for limma
boxplot(cbind(resultsMatRegular01Limma[,"totalGenesFound"],resultsMatSW01Limma[,"totalGenesFoundSW"],
	      resultsMatRegular05Limma[,"totalGenesFound"],resultsMatSW05Limma[,"totalGenesFoundSW"],
	      resultsMatRegular10Limma[,"totalGenesFound"],resultsMatSW10Limma[,"totalGenesFoundSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(fdrExtraGenesRegular01Limma,
	      fdrExtraGenesRegular05Limma,
	      fdrExtraGenesRegular10Limma),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))


### nr of genes found and FDR of extra genes for edgeR
boxplot(cbind(resultsMatRegular01[,"totalGenesFound"],resultsMatSW01[,"totalGenesFoundSW"],
	      resultsMatRegular05[,"totalGenesFound"],resultsMatSW05[,"totalGenesFoundSW"],
	      resultsMatRegular10[,"totalGenesFound"],resultsMatSW10[,"totalGenesFoundSW"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.1,.1),3), border=rep(c("black","steelblue"),3), col=alpha(rep(c("black","steelblue"),3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("bottomright",c("Standard","Stage-wise"),lty=1,col=c("black","steelblue"), bty="n", cex=2, lwd=2)

boxplot(cbind(fdrExtraGenesRegular01,
	      fdrExtraGenesRegular05,
	      fdrExtraGenesRegular10),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))




