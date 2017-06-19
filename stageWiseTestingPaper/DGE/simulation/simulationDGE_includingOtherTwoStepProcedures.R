setwd("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper_public/stageWiseTestingPaper/DGE/simulation")
source("http://130.60.190.4/robinson_lab/edgeR_robust/robust_simulation.R")
library(scales) ; library(RColorBrewer)
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

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")

resultsMatJiangDoerge01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

fdrExtraGenesRegular01 <- vector(length=N)
fdrExtraGenesRegular05 <- vector(length=N)
fdrExtraGenesRegular10 <- vector(length=N)

fdrExtraGenesRegular01Limma <- vector(length=N)
fdrExtraGenesRegular05Limma <- vector(length=N)
fdrExtraGenesRegular10Limma <- vector(length=N)


for(iter in 1:N){
	message(paste0("simulation ",iter,"\n"))
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
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fit)
	resultsMatJiangDoerge05[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fit)
	resultsMatJiangDoerge10[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fit)


	## limma analysis
	## standard analysis
	resultsMatRegular01Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
	## stage-wise analysis
	resultsMatSW01Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fitLimma)
	resultsMatJiangDoerge05Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fitLimma)
	resultsMatJiangDoerge10Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fitLimma)

	## compare extra genes found in standard to stage-wise analysis
	fdrExtraGenesRegular01[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.01,design,fit=fit)
	fdrExtraGenesRegular05[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.05,design,fit=fit)
	fdrExtraGenesRegular10[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.1,design,fit=fit)

	fdrExtraGenesRegular01Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.01,design,fit=fitLimma)
	fdrExtraGenesRegular05Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.05,design,fit=fitLimma)
	fdrExtraGenesRegular10Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.1,design,fit=fitLimma)
}

## OFDR and power interaction plot for limma
colorBrewerCols=brewer.pal(8,"Dark2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],resultsMatJiangDoerge01Limma[,"fdrAllHyp"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],resultsMatJiangDoerge05Limma[,"fdrAllHyp"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],resultsMatJiangDoerge10Limma[,"fdrAllHyp"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],resultsMatJiangDoerge01Limma[,"overallFDR"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],resultsMatJiangDoerge05Limma[,"overallFDR"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],resultsMatJiangDoerge10Limma[,"overallFDR"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge01Limma[,"nullGeneFDR"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge05Limma[,"nullGeneFDR"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge10Limma[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
text(x=0.25,y=0.12,"all hypotheses", cex=1.5)
text(x=2.8,y=0.12,"OFDR", cex=1.5)
text(x=5.4,y=0.12,"null genes", cex=1.5)
 legend(x=-0.85,y=0.115,c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
mtext(text="A",side=3,line=1,adj=0,cex=1.5)

boxplot(cbind(
	resultsMatRegular01Limma[,"powerInteraction"],resultsMatSW01Limma[,"powerInteractionSW"],resultsMatJiangDoerge01Limma[,"powerInteraction"],
	resultsMatRegular05Limma[,"powerInteraction"],resultsMatSW05Limma[,"powerInteractionSW"],resultsMatJiangDoerge05Limma[,"powerInteraction"],
	resultsMatRegular10Limma[,"powerInteraction"],resultsMatSW10Limma[,"powerInteractionSW"],resultsMatJiangDoerge10Limma[,"powerInteraction"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
mtext(text="B",side=3,line=1,adj=0,cex=1.5)


### OFDR and power interaction plot for edgeR
colorBrewerCols=brewer.pal(8,"Set2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],resultsMatJiangDoerge01[,"fdrAllHyp"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],resultsMatJiangDoerge05[,"fdrAllHyp"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],resultsMatJiangDoerge10[,"fdrAllHyp"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],resultsMatJiangDoerge01[,"overallFDR"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],resultsMatJiangDoerge05[,"overallFDR"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],resultsMatJiangDoerge10[,"overallFDR"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],resultsMatJiangDoerge01[,"nullGeneFDR"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],resultsMatJiangDoerge05[,"nullGeneFDR"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"],resultsMatJiangDoerge10[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
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
text(x=0.25,y=0.12,"all hypotheses", cex=1.5)
text(x=2.8,y=0.12,"OFDR", cex=1.5)
text(x=5.4,y=0.12,"null genes", cex=1.5)
 legend(x=-0.85,y=0.115,c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

boxplot(cbind(
	resultsMatRegular01[,"powerInteraction"],resultsMatSW01[,"powerInteractionSW"],resultsMatJiangDoerge01[,"powerInteraction"],
	resultsMatRegular05[,"powerInteraction"],resultsMatSW05[,"powerInteractionSW"],resultsMatJiangDoerge05[,"powerInteraction"],
	resultsMatRegular10[,"powerInteraction"],resultsMatSW10[,"powerInteractionSW"],resultsMatJiangDoerge10[,"powerInteraction"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

### t1 and t2 power limma
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01Limma[,"powerT1"],resultsMatSW01Limma[,"powerT1SW"],resultsMatJiangDoerge01Limma[,"powerT1"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT1SW"],resultsMatJiangDoerge05Limma[,"powerT1"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT1SW"],resultsMatJiangDoerge10Limma[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerT2"],resultsMatSW01Limma[,"powerT2SW"],resultsMatJiangDoerge01Limma[,"powerT2"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT2SW"],resultsMatJiangDoerge05Limma[,"powerT2"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT2SW"],resultsMatJiangDoerge10Limma[,"powerT2"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)


### t1 and t2 power edgeR
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01[,"powerT1"],resultsMatSW01[,"powerT1SW"], resultsMatJiangDoerge01[,"powerT1"],
	      resultsMatRegular05[,"powerT1"],resultsMatSW05[,"powerT1SW"],resultsMatJiangDoerge05[,"powerT1"],
	      resultsMatRegular10[,"powerT1"],resultsMatSW10[,"powerT1SW"],resultsMatJiangDoerge10[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerT2"],resultsMatSW01[,"powerT2SW"], resultsMatJiangDoerge01[,"powerT2"],
	      resultsMatRegular05[,"powerT2"],resultsMatSW05[,"powerT2SW"],resultsMatJiangDoerge05[,"powerT2"],
	      resultsMatRegular10[,"powerT2"],resultsMatSW10[,"powerT2SW"],resultsMatJiangDoerge10[,"powerT2"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)


### nr of genes found limma and edgeR
boxplot(cbind(resultsMatRegular01Limma[,"totalGenesFound"],resultsMatSW01Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge01Limma[,"totalGenesFound"],
	      resultsMatRegular05Limma[,"totalGenesFound"],resultsMatSW05Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge05Limma[,"totalGenesFound"],
	      resultsMatRegular10Limma[,"totalGenesFound"],resultsMatSW10Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge10Limma[,"totalGenesFound"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)

boxplot(cbind(resultsMatRegular01[,"totalGenesFound"],resultsMatSW01[,"totalGenesFoundSW"],resultsMatJiangDoerge01[,"totalGenesFound"],
	      resultsMatRegular05[,"totalGenesFound"],resultsMatSW05[,"totalGenesFoundSW"],resultsMatJiangDoerge05[,"totalGenesFound"],
	      resultsMatRegular10[,"totalGenesFound"],resultsMatSW10[,"totalGenesFoundSW"],resultsMatJiangDoerge10[,"totalGenesFound"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)



### total number of correct genes
boxplot(cbind(resultsMatRegular01Limma[,"totalGenesFound"]-resultsMatRegular01Limma[,"nrFalsePositiveNullGenes"],
	resultsMatSW01Limma[,"totalGenesFoundSW"]-resultsMatSW01Limma[,"nrFalsePositiveNullGenesSW"],
	resultsMatJiangDoerge01Limma[,"totalGenesFound"]-resultsMatJiangDoerge01Limma[,"nrFalsePositiveNullGenes"],
	resultsMatRegular05Limma[,"totalGenesFound"]-resultsMatRegular05Limma[,"nrFalsePositiveNullGenes"],
	resultsMatSW05Limma[,"totalGenesFoundSW"]-resultsMatSW05Limma[,"nrFalsePositiveNullGenesSW"],
	resultsMatJiangDoerge05Limma[,"totalGenesFound"]-resultsMatJiangDoerge05Limma[,"nrFalsePositiveNullGenes"],
	resultsMatRegular10Limma[,"totalGenesFound"]-resultsMatRegular10Limma[,"nrFalsePositiveNullGenes"],
	resultsMatSW10Limma[,"totalGenesFoundSW"]-resultsMatSW10Limma[,"nrFalsePositiveNullGenesSW"],
	resultsMatJiangDoerge10Limma[,"totalGenesFound"]-resultsMatJiangDoerge10Limma[,"nrFalsePositiveNullGenes"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Number of correct genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1, lwd=2)


### FDP of extra genes for limma and edgeR
boxplot(cbind(fdrExtraGenesRegular01Limma,
	      fdrExtraGenesRegular05Limma,
	      fdrExtraGenesRegular10Limma),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))

boxplot(cbind(fdrExtraGenesRegular01,
	      fdrExtraGenesRegular05,
	      fdrExtraGenesRegular10),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))


#### FDR of contrasts: interaction effect
par(mfrow=c(1,2))
boxplot(cbind(
	resultsMatRegular01[,"fdrInt"],resultsMatSW01[,"fdrInteractionSW"],resultsMatJiangDoerge01[,"fdrInt"],
	resultsMatRegular05[,"fdrInt"],resultsMatSW05[,"fdrInteractionSW"],resultsMatJiangDoerge05[,"fdrInt"],
	resultsMatRegular10[,"fdrInt"],resultsMatSW10[,"fdrInteractionSW"],resultsMatJiangDoerge10[,"fdrInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR interaction effect", main="edgeR, 5 vs. 5", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

boxplot(cbind(
	resultsMatRegular01Limma[,"fdrInt"],resultsMatSW01Limma[,"fdrInteractionSW"],resultsMatJiangDoerge01Limma[,"fdrInt"],
	resultsMatRegular05Limma[,"fdrInt"],resultsMatSW05Limma[,"fdrInteractionSW"],resultsMatJiangDoerge05Limma[,"fdrInt"],
	resultsMatRegular10Limma[,"fdrInt"],resultsMatSW10Limma[,"fdrInteractionSW"],resultsMatJiangDoerge10Limma[,"fdrInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR interaction effect", main="limma-voom, 5 vs. 5", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)




#### FDR of contrasts: t1
boxplot(cbind(
	resultsMatRegular01Limma[1:5,"fdrT1"],resultsMatSW01Limma[1:5,"fdrT1SW"],resultsMatJiangDoerge01Limma[1:5,"fdrT1"],
	resultsMatRegular05Limma[1:5,"fdrT1"],resultsMatSW05Limma[1:5,"fdrT1SW"],resultsMatJiangDoerge05Limma[1:5,"fdrT1"],
	resultsMatRegular10Limma[1:5,"fdrT1"],resultsMatSW10Limma[1:5,"fdrT1SW"],resultsMatJiangDoerge10Limma[1:5,"fdrT1"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR t1", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

#### FDR of contrasts: t2
boxplot(cbind(
	resultsMatRegular01Limma[1:5,"fdrT2"],resultsMatSW01Limma[1:5,"fdrT2SW"],resultsMatJiangDoerge01Limma[1:5,"fdrT2"],
	resultsMatRegular05Limma[1:5,"fdrT2"],resultsMatSW05Limma[1:5,"fdrT2SW"],resultsMatJiangDoerge05Limma[1:5,"fdrT2"],
	resultsMatRegular10Limma[1:5,"fdrT2"],resultsMatSW10Limma[1:5,"fdrT2SW"],resultsMatJiangDoerge10Limma[1:5,"fdrT2"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR t2", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

### proportion of false interaction genes without main effects
par(mfrow=c(1,2))
boxplot(cbind(
	resultsMatRegular01[,"propNullInt"],resultsMatSW01[,"propNullInt"],resultsMatJiangDoerge01[,"propNullInt"],
	resultsMatRegular05[,"propNullInt"],resultsMatSW05[,"propNullInt"],resultsMatJiangDoerge05[,"propNullInt"],
	resultsMatRegular10[,"propNullInt"],resultsMatSW10[,"propNullInt"],resultsMatJiangDoerge10[,"propNullInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Proportion false interaction genes that are null genes", main="edgeR, 5 vs. 5", bty="l", cex.lab=1.25, cex.axis=1.25)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

boxplot(cbind(
	resultsMatRegular01Limma[,"propNullInt"],resultsMatSW01Limma[,"propNullInt"],resultsMatJiangDoerge01Limma[,"propNullInt"],
	resultsMatRegular05Limma[,"propNullInt"],resultsMatSW05Limma[,"propNullInt"],resultsMatJiangDoerge05Limma[,"propNullInt"],
	resultsMatRegular10Limma[,"propNullInt"],resultsMatSW10Limma[,"propNullInt"],resultsMatJiangDoerge10Limma[,"propNullInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Proportion false interaction genes that are null genes", main="limma-voom, 5 vs. 5", bty="l", cex.lab=1.25, cex.axis=1.25)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)



####################################################
#### FDR-TPR curve on the first simulation #########
####################################################
seed=1
set.seed(seed)
nreps <- 5
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
grp <- factor(rep(0:1,each=5))
iter=seed
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


	pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
	resultsEdgeRConventional=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsEdgeRStageWise=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt
		=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsEdgeRStageWiseJiangDoerge=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaConventional=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaStageWise=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt
		=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaStageWiseJiangDoerge=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)


	## edgeR analysis
	## standard analysis
	resultsEdgeRConventional[,2:ncol(resultsEdgeRConventional)] <- t(sapply(pvalSeq, function(alpha) doRegularAnalysis(data,L,alpha=alpha,design,fit)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

	## stage-wise analysis
	resultsEdgeRStageWise[,2:ncol(resultsEdgeRStageWise)] <- t(sapply(pvalSeq, function(alpha) doStageWiseAnalysis(data,L,alpha=alpha,design,fit)[c("overallFDRSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "powerInteractionSW", "powerT1SW", "powerT2SW", "fdrAllHypSW")]))

	## Jiang & Doerge stage-wise analysis
	resultsEdgeRStageWiseJiangDoerge[,2:ncol(resultsEdgeRStageWiseJiangDoerge)] <- t(sapply(pvalSeq, function(alpha) doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*alpha,alpha2=0.2*alpha,design,fit)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))


	## limma analysis
	## standard analysis
	resultsLimmaConventional[,2:ncol(resultsLimmaConventional)] <- t(sapply(pvalSeq, function(alpha) doRegularAnalysisLimma(data,alpha=alpha,design,fitLimma)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

	## stage-wise analysis
	resultsLimmaStageWise[,2:ncol(resultsLimmaStageWise)] <- t(sapply(pvalSeq, function(alpha) doStageWiseAnalysisLimma(data,alpha=alpha,design,fitLimma)[c("overallFDRSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "powerInteractionSW", "powerT1SW", "powerT2SW", "fdrAllHypSW")]))

	## Jiang & Doerge stage-wise analysis
	resultsLimmaStageWiseJiangDoerge[,2:ncol(resultsLimmaStageWiseJiangDoerge)] <- t(sapply(pvalSeq, function(alpha) doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*alpha,alpha2=0.2*alpha,design,fitLimma)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

########### edgeR
colorBrewerCols=brewer.pal(8,"Dark2")
alphaLevels=c(508, 516, 526) #1, 5 and 10% FDR
#t1, FDP t1
library(scales)
#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
#par(mfrow=c(1,3))
plot(x=resultsEdgeRConventional$FDPt1, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 5vs5, t1, FDR hyp.", xlab="FDP on t1 contrast", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$FDPt1, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$FDPt1, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$FDPt1[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$FDPt1[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$FDPt1[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$FDPt1[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
#t1, OFDP
plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 5vs5, t1, OFDR", xlab="overall FDP", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#t1, fdr all hypotheses
plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 5vs5, t1, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()

png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_t2.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
#t2, OFDP
plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprT2,type="l", main="edgeR 5vs5, t2, OFDR", xlab="overall FDP", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#T2, fdr all hypotheses
plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprT2,type="l", main="edgeR 5vs5, t2, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()


#interaction, FDP int
#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
#par(mfrow=c(1,3))
plot(x=resultsEdgeRConventional$FDPInt, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 5vs5, interaction, FDR hyp.", xlab="FDP on interaction", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$FDPInt, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$FDPInt, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$FDPInt[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$FDPInt[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$FDPInt[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$FDPInt[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#interaction, OFDP
png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 5vs5, interaction, OFDR", xlab="Overall FDP", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#interaction, fdr all hypotheses
plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 5vs5, interaction, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()


######## limma
#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
#par(mfrow=c(1,3))
plot(x=resultsLimmaConventional$FDPt1, y=resultsLimmaConventional$tprT1,type="l", main="Limma 5vs5, t1, FDR hyp.", xlab="FDP on t1 contrast", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$FDPt1, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$FDPt1, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$FDPt1[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$FDPt1[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$FDPt1[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$FDPt1[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)

#t1, OFDP
png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprT1,type="l", main="limma-voom 5vs5, t1, OFDR", xlab="Overall FDP", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#t1, FDP over all hyp
plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprT1,type="l", main="limma-voom 5vs5, t1, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()


#t2, OFDP
png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_t2.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprT2,type="l", main="limma-voom 5vs5, t2, OFDR", xlab="Overall FDP", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#T2, FDP over all hyp
plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprT2,type="l", main="limma-voom 5vs5, t2, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()



#interaction, FDP int
#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
#par(mfrow=c(1,3))
plot(x=resultsLimmaConventional$FDPInt, y=resultsLimmaConventional$tprInt,type="l", main="Limma 5vs5, interaction, FDR hyp.", xlab="FDP on interaction", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$FDPInt, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$FDPInt, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$FDPInt[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$FDPInt[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$FDPInt[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$FDPInt[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)

png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
par(mfrow=c(1,2))
#interaction, OFDP
plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprInt,type="l", main="limma-voom 5vs5, interaction, OFDR", xlab="Overall FDP", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
#interaction, FDP all hyp
plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprInt,type="l", main="limma-voom 5vs5, interaction, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
dev.off()

###### version 2 of resultsDGE main figures
## OFDR and power interaction plot for limma
pdf("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/resultsDGE_limmaFiveReps_v2.pdf", width=11,height=7)
colorBrewerCols=brewer.pal(8,"Set2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],resultsMatJiangDoerge01Limma[,"fdrAllHyp"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],resultsMatJiangDoerge05Limma[,"fdrAllHyp"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],resultsMatJiangDoerge10Limma[,"fdrAllHyp"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],resultsMatJiangDoerge01Limma[,"overallFDR"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],resultsMatJiangDoerge05Limma[,"overallFDR"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],resultsMatJiangDoerge10Limma[,"overallFDR"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge01Limma[,"nullGeneFDR"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge05Limma[,"nullGeneFDR"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge10Limma[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
axis(2,at=c(0.01,0.05,0.1))
axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
abline(v=c(2.3,4.7),col=alpha("grey",1))
lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
text(x=0.25,y=0.12,"all hypotheses", cex=1.5)
text(x=2.8,y=0.12,"OFDR", cex=1.5)
text(x=5.4,y=0.12,"null genes", cex=1.5)
 legend(x=-0.7,y=0.1175,c("Conventional","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.5, lwd=2)
mtext(text="A",side=3,line=1,adj=0,cex=1.5)

plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprInt,type="l", main="", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=2, cex.axis=2)
abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
#working points
pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
legend("bottomright",c("Conventional","SW: Heller","SW: Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
mtext(text="B",side=3,line=1,adj=-.15,cex=1.5)
dev.off()



### OFDR and power interaction plot for edgeR
png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/resultsDGE_edgeRFiveReps_v2.png", width=11,height=7,units="in",res=300)
colorBrewerCols=brewer.pal(8,"Set2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],resultsMatJiangDoerge01[,"fdrAllHyp"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],resultsMatJiangDoerge05[,"fdrAllHyp"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],resultsMatJiangDoerge10[,"fdrAllHyp"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],resultsMatJiangDoerge01[,"overallFDR"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],resultsMatJiangDoerge05[,"overallFDR"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],resultsMatJiangDoerge10[,"overallFDR"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],resultsMatJiangDoerge01[,"nullGeneFDR"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],resultsMatJiangDoerge05[,"nullGeneFDR"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"],resultsMatJiangDoerge10[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
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
text(x=0.25,y=0.13,"all hypotheses", cex=1.5)
text(x=2.8,y=0.13,"OFDR", cex=1.5)
text(x=5.4,y=0.13,"null genes", cex=1.5)
 legend(x=-0.85,y=0.125,c("Conventional","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.5, lwd=2)

 plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprInt,type="l", main="", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=2, cex.axis=2)
 abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
 lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
 lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
 #working points
 pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
 points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
 points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
 pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
 points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
 points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
 pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
 points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
 points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
 legend("bottomright",c("Conventional","SW: Heller","SW: Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
 mtext(text="B",side=3,line=1,adj=-.15,cex=1.5)
 dev.off()






### same simulation but with three replicates.
N=30
nreps <- 3
nCst <- 2000
nT1 <- 1000
nT2 <- 1000
nInt <- 1000
grp <- factor(rep(0:1,each=3))

resultsMatRegular01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW01Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW05Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW10Limma) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")

resultsMatJiangDoerge01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge01Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge01Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge05Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge05Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge10Limma <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge10Limma) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

fdrExtraGenesRegular01 <- vector(length=N)
fdrExtraGenesRegular05 <- vector(length=N)
fdrExtraGenesRegular10 <- vector(length=N)

fdrExtraGenesRegular01Limma <- vector(length=N)
fdrExtraGenesRegular05Limma <- vector(length=N)
fdrExtraGenesRegular10Limma <- vector(length=N)


for(iter in 1:N){
	message(paste0("simulation ",iter,"\n"))
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
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fit)
	resultsMatJiangDoerge05[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fit)
	resultsMatJiangDoerge10[iter,] <- doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fit)


	## limma analysis
	## standard analysis
	resultsMatRegular01Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10Limma[iter,] <- doRegularAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
	## stage-wise analysis
	resultsMatSW01Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10Limma[iter,] <- doStageWiseAnalysisLimma(data,alpha=0.10,design,fit=fitLimma)
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fitLimma)
	resultsMatJiangDoerge05Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fitLimma)
	resultsMatJiangDoerge10Limma[iter,] <- doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fitLimma)

	## compare extra genes found in standard to stage-wise analysis
	fdrExtraGenesRegular01[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.01,design,fit=fit)
	fdrExtraGenesRegular05[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.05,design,fit=fit)
	fdrExtraGenesRegular10[iter] = compareRegularAndStageWiseAnalysisEdgeR(data,L,alpha=0.1,design,fit=fit)

	fdrExtraGenesRegular01Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.01,design,fit=fitLimma)
	fdrExtraGenesRegular05Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.05,design,fit=fitLimma)
	fdrExtraGenesRegular10Limma[iter] = compareRegularAndStageWiseAnalysisLimma(data,L,alpha=0.1,design,fit=fitLimma)

}

## OFDR and power interaction plot for limma
colorBrewerCols=brewer.pal(8,"Set2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],resultsMatJiangDoerge01Limma[,"fdrAllHyp"],
	      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],resultsMatJiangDoerge05Limma[,"fdrAllHyp"],
	      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],resultsMatJiangDoerge10Limma[,"fdrAllHyp"],
	      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],resultsMatJiangDoerge01Limma[,"overallFDR"],
	      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],resultsMatJiangDoerge05Limma[,"overallFDR"],
	      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],resultsMatJiangDoerge10Limma[,"overallFDR"],
	      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge01Limma[,"nullGeneFDR"],
	      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge05Limma[,"nullGeneFDR"],
	      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge10Limma[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
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
text(x=0.25,y=0.105,"all hypotheses", cex=1.5)
text(x=2.8,y=0.105,"OFDR", cex=1.5)
text(x=5.4,y=0.105,"null genes", cex=1.5)
 legend(x=-0.85,y=0.1,c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

boxplot(cbind(
	resultsMatRegular01Limma[,"powerInteraction"],resultsMatSW01Limma[,"powerInteractionSW"],resultsMatJiangDoerge01Limma[,"powerInteraction"],
	resultsMatRegular05Limma[,"powerInteraction"],resultsMatSW05Limma[,"powerInteractionSW"],resultsMatJiangDoerge05Limma[,"powerInteraction"],
	resultsMatRegular10Limma[,"powerInteraction"],resultsMatSW10Limma[,"powerInteractionSW"],resultsMatJiangDoerge10Limma[,"powerInteraction"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)


### OFDR and power interaction plot for edgeR
colorBrewerCols=brewer.pal(8,"Set2")
layout(matrix(c(1,1,2),nrow=1,ncol=3))
par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],resultsMatJiangDoerge01[,"fdrAllHyp"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],resultsMatJiangDoerge05[,"fdrAllHyp"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],resultsMatJiangDoerge10[,"fdrAllHyp"],
	      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],resultsMatJiangDoerge01[,"overallFDR"],
	      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],resultsMatJiangDoerge05[,"overallFDR"],
	      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],resultsMatJiangDoerge10[,"overallFDR"],
	      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],resultsMatJiangDoerge01[,"nullGeneFDR"],
	      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],resultsMatJiangDoerge05[,"nullGeneFDR"],
	      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"],resultsMatJiangDoerge10[,"nullGeneFDR"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
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
text(x=0.25,y=0.13,"all hypotheses", cex=1.5)
text(x=2.8,y=0.13,"OFDR", cex=1.5)
text(x=5.4,y=0.13,"null genes", cex=1.5)
 legend(x=-0.85,y=0.125,c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

boxplot(cbind(
	resultsMatRegular01[,"powerInteraction"],resultsMatSW01[,"powerInteractionSW"],resultsMatJiangDoerge01[,"powerInteraction"],
	resultsMatRegular05[,"powerInteraction"],resultsMatSW05[,"powerInteractionSW"],resultsMatJiangDoerge05[,"powerInteraction"],
	resultsMatRegular10[,"powerInteraction"],resultsMatSW10[,"powerInteractionSW"],resultsMatJiangDoerge10[,"powerInteraction"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power interaction effect", main="", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

### t1 and t2 power limma
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01Limma[,"powerT1"],resultsMatSW01Limma[,"powerT1SW"],resultsMatJiangDoerge01Limma[,"powerT1"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT1SW"],resultsMatJiangDoerge05Limma[,"powerT1"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT1SW"],resultsMatJiangDoerge10Limma[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)

boxplot(cbind(resultsMatRegular01Limma[,"powerT2"],resultsMatSW01Limma[,"powerT2SW"],resultsMatJiangDoerge01Limma[,"powerT2"],
	      resultsMatRegular05Limma[,"powerT1"],resultsMatSW05Limma[,"powerT2SW"],resultsMatJiangDoerge05Limma[,"powerT2"],
	      resultsMatRegular10Limma[,"powerT1"],resultsMatSW10Limma[,"powerT2SW"],resultsMatJiangDoerge10Limma[,"powerT2"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)


### t1 and t2 power edgeR
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01[,"powerT1"],resultsMatSW01[,"powerT1SW"], resultsMatJiangDoerge01[,"powerT1"],
	      resultsMatRegular05[,"powerT1"],resultsMatSW05[,"powerT1SW"],resultsMatJiangDoerge05[,"powerT1"],
	      resultsMatRegular10[,"powerT1"],resultsMatSW10[,"powerT1SW"],resultsMatJiangDoerge10[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)

boxplot(cbind(resultsMatRegular01[,"powerT2"],resultsMatSW01[,"powerT2SW"], resultsMatJiangDoerge01[,"powerT2"],
	      resultsMatRegular05[,"powerT2"],resultsMatSW05[,"powerT2SW"],resultsMatJiangDoerge05[,"powerT2"],
	      resultsMatRegular10[,"powerT2"],resultsMatSW10[,"powerT2SW"],resultsMatJiangDoerge10[,"powerT2"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t2", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)


### nr of genes found limma and edgeR
boxplot(cbind(resultsMatRegular01Limma[,"totalGenesFound"],resultsMatSW01Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge01Limma[,"totalGenesFound"],
	      resultsMatRegular05Limma[,"totalGenesFound"],resultsMatSW05Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge05Limma[,"totalGenesFound"],
	      resultsMatRegular10Limma[,"totalGenesFound"],resultsMatSW10Limma[,"totalGenesFoundSW"],resultsMatJiangDoerge10Limma[,"totalGenesFound"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)

boxplot(cbind(resultsMatRegular01[,"totalGenesFound"],resultsMatSW01[,"totalGenesFoundSW"],resultsMatJiangDoerge01[,"totalGenesFound"],
	      resultsMatRegular05[,"totalGenesFound"],resultsMatSW05[,"totalGenesFoundSW"],resultsMatJiangDoerge05[,"totalGenesFound"],
	      resultsMatRegular10[,"totalGenesFound"],resultsMatSW10[,"totalGenesFoundSW"],resultsMatJiangDoerge10[,"totalGenesFound"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Number of genes found", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)



### FDR of extra genes for limma and edgeR
boxplot(cbind(fdrExtraGenesRegular01Limma,
	      fdrExtraGenesRegular05Limma,
	      fdrExtraGenesRegular10Limma),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))


boxplot(cbind(fdrExtraGenesRegular01,
	      fdrExtraGenesRegular05,
	      fdrExtraGenesRegular10),
	      boxwex=.2, at=c(0.3,1.1,1.9), col=alpha("black",.2), xaxt="n", ylab="FDP of the extra genes", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))


#### FDR of contrasts: interaction effect
par(mfrow=c(1,2))
boxplot(cbind(
	resultsMatRegular01[,"fdrInt"],resultsMatSW01[,"fdrInteractionSW"],resultsMatJiangDoerge01[,"fdrInt"],
	resultsMatRegular05[,"fdrInt"],resultsMatSW05[,"fdrInteractionSW"],resultsMatJiangDoerge05[,"fdrInt"],
	resultsMatRegular10[,"fdrInt"],resultsMatSW10[,"fdrInteractionSW"],resultsMatJiangDoerge10[,"fdrInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR interaction effect", main="edgeR, 3 vs. 3", bty="l", ylim=c(0,0.17))
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

boxplot(cbind(
	resultsMatRegular01Limma[,"fdrInt"],resultsMatSW01Limma[,"fdrInteractionSW"],resultsMatJiangDoerge01Limma[,"fdrInt"],
	resultsMatRegular05Limma[,"fdrInt"],resultsMatSW05Limma[,"fdrInteractionSW"],resultsMatJiangDoerge05Limma[,"fdrInt"],
	resultsMatRegular10Limma[,"fdrInt"],resultsMatSW10Limma[,"fdrInteractionSW"],resultsMatJiangDoerge10Limma[,"fdrInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="FDR interaction effect", main="limma-voom, 3 vs. 3", bty="l")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)


### null genes in false positive interaction genes
par(mfrow=c(1,2))
boxplot(cbind(
	resultsMatRegular01[,"propNullInt"],resultsMatSW01[,"propNullInt"],resultsMatJiangDoerge01[,"propNullInt"],
	resultsMatRegular05[,"propNullInt"],resultsMatSW05[,"propNullInt"],resultsMatJiangDoerge05[,"propNullInt"],
	resultsMatRegular10[,"propNullInt"],resultsMatSW10[,"propNullInt"],resultsMatJiangDoerge10[,"propNullInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Proportion false interaction genes that are null genes", main="edgeR, 3 vs. 3", bty="l", cex.lab=1.25, cex.axis=1.25)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)

boxplot(cbind(
	resultsMatRegular01Limma[,"propNullInt"],resultsMatSW01Limma[,"propNullInt"],resultsMatJiangDoerge01Limma[,"propNullInt"],
	resultsMatRegular05Limma[,"propNullInt"],resultsMatSW05Limma[,"propNullInt"],resultsMatJiangDoerge05Limma[,"propNullInt"],
	resultsMatRegular10Limma[,"propNullInt"],resultsMatSW10Limma[,"propNullInt"],resultsMatJiangDoerge10Limma[,"propNullInt"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Proportion false interaction genes that are null genes", main="limma-voom, 3 vs. 3", bty="l", cex.lab=1.25, cex.axis=1.25)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.lab=1.25, cex.axis=1.25)
abline(v=c(0.7,1.5),col=alpha("grey",1))
legend("bottomright",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.25, lwd=2)




####################################################
#### FDR-TPR curve on the first simulation #########
####################################################
iter=1
set.seed(iter)
nreps=3
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


	pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
	resultsEdgeRConventional=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsEdgeRStageWise=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt
		=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsEdgeRStageWiseJiangDoerge=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaConventional=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaStageWise=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt
		=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)
	resultsLimmaStageWiseJiangDoerge=data.frame(alpha=pvalSeq, OFDP=NA, FDPInt=NA, FDPt1=NA, FDPt2=NA, tprInt=NA, tprT1=NA, tprT2=NA, fdpAll=NA)


	## edgeR analysis
	## standard analysis
	resultsEdgeRConventional[,2:ncol(resultsEdgeRConventional)] <- t(sapply(pvalSeq, function(alpha) doRegularAnalysis(data,L,alpha=alpha,design,fit)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

	## stage-wise analysis
	resultsEdgeRStageWise[,2:ncol(resultsEdgeRStageWise)] <- t(sapply(pvalSeq, function(alpha) doStageWiseAnalysis(data,L,alpha=alpha,design,fit)[c("overallFDRSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "powerInteractionSW", "powerT1SW", "powerT2SW", "fdrAllHypSW")]))

	## Jiang & Doerge stage-wise analysis
	resultsEdgeRStageWiseJiangDoerge[,2:ncol(resultsEdgeRStageWiseJiangDoerge)] <- t(sapply(pvalSeq, function(alpha) doStageWiseJiangDoergeAnalysis(data,L,alpha1=0.8*alpha,alpha2=0.2*alpha,design,fit)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))


	## limma analysis
	## standard analysis
	resultsLimmaConventional[,2:ncol(resultsLimmaConventional)] <- t(sapply(pvalSeq, function(alpha) doRegularAnalysisLimma(data,alpha=alpha,design,fitLimma)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

	## stage-wise analysis
	resultsLimmaStageWise[,2:ncol(resultsLimmaStageWise)] <- t(sapply(pvalSeq, function(alpha) doStageWiseAnalysisLimma(data,alpha=alpha,design,fitLimma)[c("overallFDRSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "powerInteractionSW", "powerT1SW", "powerT2SW", "fdrAllHypSW")]))

	## Jiang & Doerge stage-wise analysis
	resultsLimmaStageWiseJiangDoerge[,2:ncol(resultsLimmaStageWiseJiangDoerge)] <- t(sapply(pvalSeq, function(alpha) doStageWiseJiangDoergeAnalysisLimma(data,L,alpha1=0.8*alpha,alpha2=0.2*alpha,design,fitLimma)[c("overallFDR", "fdrInteraction", "fdrT1", "fdrT2", "powerInteraction", "powerT1", "powerT2", "fdrAllHyp")]))

	########### edgeR
	colorBrewerCols=brewer.pal(8,"Dark2")
	alphaLevels=c(508, 516, 526) #1, 5 and 10% FDR
	#t1, FDP t1
	library(scales)
	#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
	#par(mfrow=c(1,3))
	plot(x=resultsEdgeRConventional$FDPt1, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 3vs3, t1, FDR hyp.", xlab="FDP on t1 contrast", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$FDPt1, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$FDPt1, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$FDPt1[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$FDPt1[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$FDPt1[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$FDPt1[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$FDPt1[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_threeSamples_t1.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	#t1, OFDP
	plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 3vs3, t1, OFDR", xlab="overall FDP", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#t1, fdr all hypotheses
	plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprT1,type="l", main="edgeR 3vs3, t1, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()

	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_threeSamples_t2.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	#t2, OFDP
	plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprT2,type="l", main="edgeR 3vs3, t2, OFDR", xlab="overall FDP", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#T2, fdr all hypotheses
	plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprT2,type="l", main="edgeR 3vs3, t2, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()


	#interaction, FDP int
	#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
	#par(mfrow=c(1,3))
	plot(x=resultsEdgeRConventional$FDPInt, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 5vs5, interaction, FDR hyp.", xlab="FDP on interaction", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$FDPInt, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$FDPInt, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$FDPInt[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$FDPInt[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$FDPInt[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$FDPInt[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$FDPInt[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#interaction, OFDP
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_edgeR_threeSamples_interaction.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	plot(x=resultsEdgeRConventional$OFDP, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 3vs3, interaction, OFDR", xlab="Overall FDP", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$OFDP, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$OFDP, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$OFDP[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$OFDP[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$OFDP[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#interaction, fdr all hypotheses
	plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprInt,type="l", main="edgeR 3vs3, interaction, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()


	######## limma
	#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_t1.png", width=11,height=7,units="in",res=300)
	#par(mfrow=c(1,3))
	plot(x=resultsLimmaConventional$FDPt1, y=resultsLimmaConventional$tprT1,type="l", main="Limma 5vs5, t1, FDR hyp.", xlab="FDP on t1 contrast", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$FDPt1, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$FDPt1, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$FDPt1[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$FDPt1[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$FDPt1[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$FDPt1[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$FDPt1[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)

	#t1, OFDP
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_threeSamples_t1.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprT1,type="l", main="limma-voom 3vs3, t1, OFDR", xlab="Overall FDP", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#t1, FDP over all hyp
	plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprT1,type="l", main="limma-voom 3vs3, t1, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t1 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprT1,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprT1,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT1[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()


	#t2, OFDP
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_threeSamples_t2.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprT2,type="l", main="limma-voom 3vs3, t2, OFDR", xlab="Overall FDP", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#T2, FDP over all hyp
	plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprT2,type="l", main="limma-voom 3vs3, t2, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR t2 contrast", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprT2,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprT2,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprT2[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()



	#interaction, FDP int
	#png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_fiveSamples_interaction.png", width=11,height=7,units="in",res=300)
	#par(mfrow=c(1,3))
	plot(x=resultsLimmaConventional$FDPInt, y=resultsLimmaConventional$tprInt,type="l", main="Limma 5vs5, interaction, FDR hyp.", xlab="FDP on interaction", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$FDPInt, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$FDPInt, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$FDPInt[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$FDPInt[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$FDPInt[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$FDPInt[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$FDPInt[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)

	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/fdrTprDGE_limma_threeSamples_interaction.png", width=11,height=7,units="in",res=300)
	par(mfrow=c(1,2))
	#interaction, OFDP
	plot(x=resultsLimmaConventional$OFDP, y=resultsLimmaConventional$tprInt,type="l", main="limma-voom 3vs3, interaction, OFDR", xlab="Overall FDP", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$OFDP, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$OFDP, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$OFDP[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$OFDP[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$OFDP[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$OFDP[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$OFDP[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	#interaction, FDP all hyp
	plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprInt,type="l", main="limma-voom 3vs3, interaction, FDR all hyp.", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=1.25, cex.axis=1.25)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW Heller","SW Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	dev.off()

	###### version 2 of resultsDGE main figures
	## OFDR and power interaction plot for limma
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/resultsDGE_limmaThreeReps_v2.png", width=11,height=7,units="in",res=300)
	colorBrewerCols=brewer.pal(8,"Set2")
	layout(matrix(c(1,1,2),nrow=1,ncol=3))
	par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
	boxplot(cbind(resultsMatRegular01Limma[,"fdrAllHyp"],resultsMatSW01Limma[,"fdrAllHypSW"],resultsMatJiangDoerge01Limma[,"fdrAllHyp"],
		      resultsMatRegular05Limma[,"fdrAllHyp"],resultsMatSW05Limma[,"fdrAllHypSW"],resultsMatJiangDoerge05Limma[,"fdrAllHyp"],
		      resultsMatRegular10Limma[,"fdrAllHyp"],resultsMatSW10Limma[,"fdrAllHypSW"],resultsMatJiangDoerge10Limma[,"fdrAllHyp"],
		      resultsMatRegular01Limma[,"overallFDR"],resultsMatSW01Limma[,"overallFDRSW"],resultsMatJiangDoerge01Limma[,"overallFDR"],
		      resultsMatRegular05Limma[,"overallFDR"],resultsMatSW05Limma[,"overallFDRSW"],resultsMatJiangDoerge05Limma[,"overallFDR"],
		      resultsMatRegular10Limma[,"overallFDR"],resultsMatSW10Limma[,"overallFDRSW"],resultsMatJiangDoerge10Limma[,"overallFDR"],
		      resultsMatRegular01Limma[,"nullGeneFDR"],resultsMatSW01Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge01Limma[,"nullGeneFDR"],
		      resultsMatRegular05Limma[,"nullGeneFDR"],resultsMatSW05Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge05Limma[,"nullGeneFDR"],
		      resultsMatRegular10Limma[,"nullGeneFDR"],resultsMatSW10Limma[,"nullGeneFDRSW"],resultsMatJiangDoerge10Limma[,"nullGeneFDR"]),
		boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
	axis(2,at=c(0.01,0.05,0.1))
	axis(1,at=c(seq(0.3,6.7,by=0.8)), labels=rep(c("1%","5%","10%"),3))
	abline(v=c(2.3,4.7),col=alpha("grey",1))
	lines(x=c(0.1,0.5),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
	lines(x=c(0.9,1.3),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
	lines(x=c(1.7,2.1),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
	lines(x=c(2.5,2.9),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
	lines(x=c(3.3,3.7),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
	lines(x=c(4.1,4.5),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
	lines(x=c(4.9,5.3),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
	lines(x=c(5.7,6.1),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
	lines(x=c(6.5,6.9),y=rep(0.1,each=2),col=2,lty=3, lwd=2)
	text(x=0.25,y=0.12,"all hypotheses", cex=1.5)
	text(x=2.8,y=0.12,"OFDR", cex=1.5)
	text(x=5.4,y=0.12,"null genes", cex=1.5)
	 legend(x=-0.7,y=0.1,c("Conventional","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.5, lwd=2)
	mtext(text="A",side=3,line=1,adj=0,cex=1.5)

	plot(x=resultsLimmaConventional$fdpAll, y=resultsLimmaConventional$tprInt,type="l", main="", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=2, cex.axis=2)
	abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	lines(x=resultsLimmaStageWise$fdpAll, y=resultsLimmaStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	lines(x=resultsLimmaStageWiseJiangDoerge$fdpAll, y=resultsLimmaStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	#working points
	pchHlp = resultsLimmaConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaConventional$fdpAll[alphaLevels],y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaConventional$fdpAll[alphaLevels], y=resultsLimmaConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	pchHlp = resultsLimmaStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels],y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWise$fdpAll[alphaLevels], y=resultsLimmaStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	pchHlp = resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	points(x=resultsLimmaStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsLimmaStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	legend("bottomright",c("Conventional","SW: Heller","SW: Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	mtext(text="B",side=3,line=1,adj=-.15,cex=1.5)
	dev.off()



	### OFDR and power interaction plot for edgeR
	png("~/Dropbox/phdKoen/stagewisetesting/figures/supplementary/resultsDGE_edgeRThreeReps_v2.png", width=11,height=7,units="in",res=300)
	colorBrewerCols=brewer.pal(8,"Set2")
	layout(matrix(c(1,1,2),nrow=1,ncol=3))
	par(bty="l", cex.axis=2, cex.lab=2, mar=c(5,4.4,4,2)+0.1)
	boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],resultsMatJiangDoerge01[,"fdrAllHyp"],
		      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],resultsMatJiangDoerge05[,"fdrAllHyp"],
		      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],resultsMatJiangDoerge10[,"fdrAllHyp"],
		      resultsMatRegular01[,"overallFDR"],resultsMatSW01[,"overallFDRSW"],resultsMatJiangDoerge01[,"overallFDR"],
		      resultsMatRegular05[,"overallFDR"],resultsMatSW05[,"overallFDRSW"],resultsMatJiangDoerge05[,"overallFDR"],
		      resultsMatRegular10[,"overallFDR"],resultsMatSW10[,"overallFDRSW"],resultsMatJiangDoerge10[,"overallFDR"],
		      resultsMatRegular01[,"nullGeneFDR"],resultsMatSW01[,"nullGeneFDRSW"],resultsMatJiangDoerge01[,"nullGeneFDR"],
		      resultsMatRegular05[,"nullGeneFDR"],resultsMatSW05[,"nullGeneFDRSW"],resultsMatJiangDoerge05[,"nullGeneFDR"],
		      resultsMatRegular10[,"nullGeneFDR"],resultsMatSW10[,"nullGeneFDRSW"],resultsMatJiangDoerge10[,"nullGeneFDR"]),
		boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=3)+rep(c(-.25,0,.25),9), border=rep(colorBrewerCols[c(3,1,2)],9), col=alpha(rep(colorBrewerCols[c(3,1,2)],9),.2), xaxt="n", yaxt="n", ylab="Empirical false discovery proportion", main="", xlab="False discovery rate cut-off")
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
	text(x=0.25,y=0.13,"all hypotheses", cex=1.5)
	text(x=2.8,y=0.13,"OFDR", cex=1.5)
	text(x=5.4,y=0.13,"null genes", cex=1.5)
	 legend(x=-0.85,y=0.125,c("Conventional","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", cex=1.5, lwd=2)

	 plot(x=resultsEdgeRConventional$fdpAll, y=resultsEdgeRConventional$tprInt,type="l", main="", xlab="FDP over all hypotheses", ylab="TPR interaction", lwd=2, col=colorBrewerCols[3], cex.lab=2, cex.axis=2)
	 abline(v=c(0.01,0.05,0.1),col=alpha("grey",0.6), lwd=2, lty=2)
	 lines(x=resultsEdgeRStageWise$fdpAll, y=resultsEdgeRStageWise$tprInt,col=colorBrewerCols[1], lwd=2)
	 lines(x=resultsEdgeRStageWiseJiangDoerge$fdpAll, y=resultsEdgeRStageWiseJiangDoerge$tprInt,col=colorBrewerCols[2], lwd=2)
	 #working points
	 pchHlp = resultsEdgeRConventional$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	 points(x=resultsEdgeRConventional$fdpAll[alphaLevels],y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	 points(x=resultsEdgeRConventional$fdpAll[alphaLevels], y=resultsEdgeRConventional$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], cex=1.5, col=colorBrewerCols[3])
	 pchHlp = resultsEdgeRStageWise$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	 points(x=resultsEdgeRStageWise$fdpAll[alphaLevels],y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	 points(x=resultsEdgeRStageWise$fdpAll[alphaLevels], y=resultsEdgeRStageWise$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[1], cex=1.5)
	 pchHlp = resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels] <= c(0.01,0.05,0.1)
	 points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels],y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(19,0)[pchHlp+1],col="white")
	 points(x=resultsEdgeRStageWiseJiangDoerge$fdpAll[alphaLevels], y=resultsEdgeRStageWiseJiangDoerge$tprInt[alphaLevels], pch=c(1,19)[pchHlp+1], col=colorBrewerCols[2], cex=1.5)
	 legend("bottomright",c("Conventional","SW: Heller","SW: Jiang & Doerge"),col=colorBrewerCols[c(3,1,2)],lty=1,lwd=2, bty="n", cex=1.5)
	 mtext(text="B",side=3,line=1,adj=-.15,cex=1.5)
	 dev.off()





####################################################################
############### ADD TEST FOR AVERAGE DE ############################
####################################################################
N=3
nreps <- 5
nCst <- 5000
nT1 <- 100
nT2 <- 1000
nInt <- 100
grp <- factor(rep(0:1,each=5))

resultsMatRegular01Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular01Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatRegular05Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular05Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatRegular10Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular10Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatRegular01LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular01LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatRegular05LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular05LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatRegular10LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatRegular10LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")

resultsMatSW01Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW01Avg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")
resultsMatSW05Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW05Avg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")
resultsMatSW10Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW10Avg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")
resultsMatSW01LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW01LimmaAvg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")
resultsMatSW05LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW05LimmaAvg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")
resultsMatSW10LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatSW10LimmaAvg) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt", "fdrAvgSW", "powerAvgSW")

resultsMatJiangDoerge01Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge01Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatJiangDoerge05Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge05Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatJiangDoerge10Avg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge10Avg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatJiangDoerge01LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge01LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatJiangDoerge05LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge05LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")
resultsMatJiangDoerge10LimmaAvg <- matrix(NA,nrow=N,ncol=15)
colnames(resultsMatJiangDoerge10LimmaAvg) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt", "fdrAvg", "powerAvg")


for(iter in 1:N){
	message(paste0("simulation ",iter,"\n"))
	set.seed(iter)
	libSize = sample(round(seq(15e6,20e6,length.out=nreps*4))) #high libSize
	nTags <- 13e3

	#Fold changes: t1 unique (50%), t1 unique (50%), t1 cst (50%), t1 cst (50%), interaction
	#foldDiffT1 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1/1.75),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	#foldDiffT2 <- unlist(mapply(rep,c(3,1/3, 3,1/3, 1.75),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))
	foldDiffT1 <- unlist(mapply(rep,c(5,1/5, 2,1/2, 1/1.25),c(nT1/2,nT1/2, nCst/2,nCst/2, nInt)))
	foldDiffT2 <- unlist(mapply(rep,c(5,1/5, 2,1/2, 1.25),c(nT2/2,nT2/2, nCst/2,nCst/2, nInt)))

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
	indAvgAll <- c(indT1All,indT2All) #also unique in t1/t2 is avg.
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
	L <- matrix(0,nrow=4,ncol=4)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t1","t2","interaction", "avgDE")
	L[2,1]=1
	L[c(2,4),2]=1
	L[4,3]=1
	L[,4]=(L[,1]+L[,2])/2

	##limma prep.
	v <- voom(d,design, plot=FALSE)
	fitLimma <- lmFit(v,design)
	colnames(design)[4] <- "treatTime2"
	contrast.matrix <- makeContrasts(t1=treattreatment,
					 t2=treattreatment+treatTime2,
					 interaction=treatTime2,
					 avgDE=treattreatment+0.5*treatTime2, levels=design)
	fitLimma <- contrasts.fit(fitLimma,contrast.matrix)
	fitLimma <- eBayes(fitLimma)

	## edgeR analysis
	## standard analysis
	resultsMatRegular01Avg[iter,] <- doRegularAnalysisAvg(data,L,alpha=0.01,design,fit)
	resultsMatRegular05Avg[iter,] <- doRegularAnalysisAvg(data,L,alpha=0.05,design,fit)
	resultsMatRegular10Avg[iter,] <- doRegularAnalysisAvg(data,L,alpha=0.10,design,fit)
	## stage-wise analysis
	resultsMatSW01Avg[iter,] <- doStageWiseAnalysisAvg(data,L,alpha=0.01,design,fit)
	resultsMatSW05Avg[iter,] <- doStageWiseAnalysisAvg(data,L,alpha=0.05,design,fit)
	resultsMatSW10Avg[iter,] <- doStageWiseAnalysisAvg(data,L,alpha=0.10,design,fit)
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01Avg[iter,] <- doStageWiseJiangDoergeAnalysisAvg(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fit)
	resultsMatJiangDoerge05Avg[iter,] <- doStageWiseJiangDoergeAnalysisAvg(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fit)
	resultsMatJiangDoerge10Avg[iter,] <- doStageWiseJiangDoergeAnalysisAvg(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fit)


	## limma analysis
	## standard analysis
	resultsMatRegular01LimmaAvg[iter,] <- doRegularAnalysisLimmaAvg(data,alpha=0.01,design,fit=fitLimma)
	resultsMatRegular05LimmaAvg[iter,] <- doRegularAnalysisLimmaAvg(data,alpha=0.05,design,fit=fitLimma)
	resultsMatRegular10LimmaAvg[iter,] <- doRegularAnalysisLimmaAvg(data,alpha=0.10,design,fit=fitLimma)
	## stage-wise analysis
	resultsMatSW01LimmaAvg[iter,] <- doStageWiseAnalysisLimmaAvg(data,alpha=0.01,design,fit=fitLimma)
	resultsMatSW05LimmaAvg[iter,] <- doStageWiseAnalysisLimmaAvg(data,alpha=0.05,design,fit=fitLimma)
	resultsMatSW10LimmaAvg[iter,] <- doStageWiseAnalysisLimmaAvg(data,alpha=0.10,design,fit=fitLimma)
	## Jiang & Doerge stage-wise analysis
	resultsMatJiangDoerge01LimmaAvg[iter,] <- doStageWiseJiangDoergeAnalysisLimmaAvg(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fitLimma)
	resultsMatJiangDoerge05LimmaAvg[iter,] <- doStageWiseJiangDoergeAnalysisLimmaAvg(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fitLimma)
	resultsMatJiangDoerge10LimmaAvg[iter,] <- doStageWiseJiangDoergeAnalysisLimmaAvg(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fitLimma)

}

###### compare results with regular 5 vs. 5 simulation
### power t1
par(mfrow=c(1,2))
boxplot(cbind(resultsMatRegular01[,"powerT1"],resultsMatSW01[,"powerT1SW"], resultsMatJiangDoerge01[,"powerT1"],
	      resultsMatRegular05[,"powerT1"],resultsMatSW05[,"powerT1SW"],resultsMatJiangDoerge05[,"powerT1"],
	      resultsMatRegular10[,"powerT1"],resultsMatSW10[,"powerT1SW"],resultsMatJiangDoerge10[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)

boxplot(cbind(resultsMatRegular01Avg[,"powerT1"],resultsMatSW01Avg[,"powerT1SW"], resultsMatJiangDoerge01Avg[,"powerT1"],
	      resultsMatRegular05Avg[,"powerT1"],resultsMatSW05Avg[,"powerT1SW"],resultsMatJiangDoerge05Avg[,"powerT1"],
	      resultsMatRegular10Avg[,"powerT1"],resultsMatSW10Avg[,"powerT1SW"],resultsMatJiangDoerge10Avg[,"powerT1"]),
	      boxwex=.2, at=rep(c(0.3,1.1,1.9),each=3)+rep(c(-.25,0,.25),3), border=rep(colorBrewerCols[c(3,1,2)],3), col=alpha(rep(colorBrewerCols[c(3,1,2)],3),.2), xaxt="n", ylab="Power t1", main="")
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"))
abline(v=c(0.7,1.5),col=alpha("grey",.8))
legend("topleft",c("Standard","SW: Heller","SW: Jiang & Doerge"),lty=1,col=colorBrewerCols[c(3,1,2)], bty="n", lwd=2)


par(mfrow=c(1,3))

boxplot(cbind(resultsMatJiangDoerge01[1:3,"powerInteraction"],resultsMatJiangDoerge01Avg[1:3,"powerInteraction"],
	resultsMatJiangDoerge05[1:3,"powerInteraction"],resultsMatJiangDoerge05Avg[1:3,"powerInteraction"],
	resultsMatJiangDoerge10[1:3,"powerInteraction"],resultsMatJiangDoerge10Avg[1:3,"powerInteraction"]
	), main="Jiang Doerge, interaction")

boxplot(cbind(resultsMatJiangDoerge01[1:3,"powerT1"],resultsMatJiangDoerge01Avg[1:3,"powerT1"],
	resultsMatJiangDoerge05[1:3,"powerT1"],resultsMatJiangDoerge05Avg[1:3,"powerT1"],
	resultsMatJiangDoerge10[1:3,"powerT1"],resultsMatJiangDoerge10Avg[1:3,"powerT1"]
	))

boxplot(cbind(resultsMatJiangDoerge01[1:3,"powerT2"],resultsMatJiangDoerge01Avg[1:3,"powerT2"],
	resultsMatJiangDoerge05[1:3,"powerT2"],resultsMatJiangDoerge05Avg[1:3,"powerT2"],
	resultsMatJiangDoerge10[1:3,"powerT2"],resultsMatJiangDoerge10Avg[1:3,"powerT2"]
	))

boxplot(cbind(resultsMatSW01[1:3,"powerT1SW"],resultsMatSW01Avg[1:3,"powerT1SW"],
	resultsMatSW05[1:3,"powerT1SW"],resultsMatSW05Avg[1:3,"powerT1SW"],
	resultsMatSW10[1:3,"powerT1SW"],resultsMatSW10Avg[1:3,"powerT1SW"]
	))
