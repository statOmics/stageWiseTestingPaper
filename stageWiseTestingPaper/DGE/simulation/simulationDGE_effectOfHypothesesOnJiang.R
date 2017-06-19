source("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/githubPaper_public/stageWiseTestingPaper/DGE/simulation/simulationDGE_helpFunctions.R")
N=30
nreps <- 5
nCst <- 0
nT1 <- 0 # 0 or 3000
nT2 <- 1000
nInt <- 1000
grp <- factor(rep(0:1,each=5))

resultsMatRegular01_t2Int <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular01_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05_t2Int <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular05_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10_t2Int <- matrix(NA,nrow=N,ncol=11)
colnames(resultsMatRegular10_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
 resultsMatSW01_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatSW01_t2Int) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT2SW", "fdrAllHypSW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
 resultsMatSW05_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatSW05_t2Int) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT2SW", "fdrAllHypSW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
 resultsMatSW10_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatSW10_t2Int) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT2SW", "fdrAllHypSW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
 resultsMatJiangDoerge01_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatJiangDoerge01_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
  resultsMatJiangDoerge05_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatJiangDoerge05_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
 resultsMatJiangDoerge10_t2Int <- matrix(NA,nrow=N,ncol=11)
 colnames(resultsMatJiangDoerge10_t2Int) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT2", "fdrInt", "fdrAllHyp", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

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
	L <- matrix(0,nrow=4,ncol=2)
	rownames(L) <- colnames(fit$coefficients)
	colnames(L) <- c("t2","interaction")
	L[c(2,4),1]=1
	L[4,2]=1
	
	# ##limma prep.
	# v <- voom(d,design, plot=FALSE)
	# fitLimma <- lmFit(v,design)
	# colnames(design)[4] <- "treatTime2"
	# contrast.matrix <- makeContrasts(t1=treattreatment, 
	# 				 t2=treattreatment+treatTime2,
	# 				 interaction=treatTime2, levels=design)
	# fitLimma <- contrasts.fit(fitLimma,contrast.matrix)
	# fitLimma <- eBayes(fitLimma)
	
	## edgeR analysis
	## standard analysis
	resultsMatRegular01_t2Int[iter,] <- doRegularAnalysis_t2Int(data,L,alpha=0.01,design,fit)
	resultsMatRegular05_t2Int[iter,] <- doRegularAnalysis_t2Int(data,L,alpha=0.05,design,fit)
	resultsMatRegular10_t2Int[iter,] <- doRegularAnalysis_t2Int(data,L,alpha=0.10,design,fit)	
	## stage-wise analysis
	 resultsMatSW01_t2Int[iter,] <- doStageWiseAnalysis_t2Int(data,L,alpha=0.01,design,fit)
	 resultsMatSW05_t2Int[iter,] <- doStageWiseAnalysis_t2Int(data,L,alpha=0.05,design,fit)
	 resultsMatSW10_t2Int[iter,] <- doStageWiseAnalysis_t2Int(data,L,alpha=0.10,design,fit)
	# ## Jiang & Doerge stage-wise analysis
	 resultsMatJiangDoerge01_t2Int[iter,] <- doStageWiseJiangDoergeAnalysis_t2Int(data,L,alpha1=0.8*0.01,alpha2=0.2*0.01,design,fit)
	 resultsMatJiangDoerge05_t2Int[iter,] <- doStageWiseJiangDoergeAnalysis_t2Int(data,L,alpha1=0.8*0.05,alpha2=0.2*0.05,design,fit)
	 resultsMatJiangDoerge10_t2Int[iter,] <- doStageWiseJiangDoergeAnalysis_t2Int(data,L,alpha1=0.8*0.1,alpha2=0.2*0.1,design,fit)
}

#store if no T1 genes
resultsMatJiangDoerge01_noT1 = resultsMatJiangDoerge01_t2Int
resultsMatJiangDoerge05_noT1 = resultsMatJiangDoerge05_t2Int
resultsMatJiangDoerge10_noT1 = resultsMatJiangDoerge10_t2Int
resultsMatRegular01_noT1 = resultsMatRegular01_t2Int
resultsMatRegular05_noT1 = resultsMatRegular05_t2Int
resultsMatRegular10_noT1 = resultsMatRegular10_t2Int
resultsMatSW01_noT1 = resultsMatSW01_t2Int
resultsMatSW05_noT1 = resultsMatSW05_t2Int
resultsMatSW10_noT1 = resultsMatSW10_t2Int
#store if 3000 t1 genes
resultsMatJiangDoerge01_3000T1 = resultsMatJiangDoerge01_t2Int
resultsMatJiangDoerge05_3000T1 = resultsMatJiangDoerge05_t2Int
resultsMatJiangDoerge10_3000T1 = resultsMatJiangDoerge10_t2Int
resultsMatRegular01_3000T1 = resultsMatRegular01_t2Int
resultsMatRegular05_3000T1 = resultsMatRegular05_t2Int
resultsMatRegular10_3000T1 = resultsMatRegular10_t2Int
resultsMatSW01_3000T1 = resultsMatSW01_t2Int
resultsMatSW05_3000T1 = resultsMatSW05_t2Int
resultsMatSW10_3000T1 = resultsMatSW10_t2Int


boxplot(cbind(resultsMatJiangDoerge01[1:3,"powerInteraction"],resultsMatJiangDoerge01_t2Int[1:3,"powerInteraction"],
	resultsMatJiangDoerge05[1:3,"powerInteraction"],resultsMatJiangDoerge05_t2Int[1:3,"powerInteraction"],
	resultsMatJiangDoerge10[1:3,"powerInteraction"],resultsMatJiangDoerge10_t2Int[1:3,"powerInteraction"]
	))


boxplot(cbind(resultsMatJiangDoerge01[1:3,"powerT2"],resultsMatJiangDoerge01_t2Int[1:3,"powerT2"],
	resultsMatJiangDoerge05[1:3,"powerT2"],resultsMatJiangDoerge05_t2Int[1:3,"powerT2"],
	resultsMatJiangDoerge10[1:3,"powerT2"],resultsMatJiangDoerge10_t2Int[1:3,"powerT2"]
	))



par(mfrow=c(1,3))
boxplot(cbind(resultsMatRegular01_noT1[,"powerT2"],resultsMatRegular01_3000T1[,"powerT2"],
	resultsMatRegular05_noT1[,"powerT2"],resultsMatRegular05_3000T1[,"powerT2"],
	resultsMatRegular10_noT1[,"powerT2"],resultsMatRegular10_3000T1[,"powerT2"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.12,.12),3), border=colorBrewerCols[3], col=alpha(colorBrewerCols[3],.2), xaxt="n", ylab="Power t2", main="Conventional, edgeR, 5 vs. 5", bty="l", cex.axis=1.25, cex.lab=1.5, cex.main=1.5
	)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.axis=1.25, cex.lab=1.5)
abline(v=c(0.7,1.5),col=alpha("grey",1))
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)


boxplot(cbind(resultsMatSW01_noT1[,"powerT2SW"],resultsMatSW01_3000T1[,"powerT2SW"],
	resultsMatSW05_noT1[,"powerT2SW"],resultsMatSW05_3000T1[,"powerT2SW"],
	resultsMatSW10_noT1[,"powerT2SW"],resultsMatSW10_3000T1[,"powerT2SW"]
	),
boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.12,.12),3), border=colorBrewerCols[1], col=alpha(colorBrewerCols[1],.2), xaxt="n", ylab="Power t2", main="Heller, edgeR, 5 vs. 5", bty="l", cex.axis=1.25, cex.lab=1.5, cex.main=1.5)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.axis=1.25, cex.lab=1.5)
abline(v=c(0.7,1.5),col=alpha("grey",1))
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)

boxplot(cbind(resultsMatJiangDoerge01_noT1[,"powerT2"],resultsMatJiangDoerge01_3000T1[,"powerT2"],
	resultsMatJiangDoerge05_noT1[,"powerT2"],resultsMatJiangDoerge05_3000T1[,"powerT2"],
	resultsMatJiangDoerge10_noT1[,"powerT2"],resultsMatJiangDoerge10_3000T1[,"powerT2"]),
	boxwex=.2, at=rep(c(0.3,1.1,1.9),each=2)+rep(c(-.12,.12),3), border=colorBrewerCols[2], col=alpha(colorBrewerCols[2],.2), xaxt="n", ylab="Power t2", main="Jiang & Doerge, edgeR, 5 vs. 5", bty="l", cex.axis=1.25, cex.lab=1.5, cex.main=1.5
	)
axis(1,at=c(0.3,1.1,1.9), labels=c("1%","5%","10%"), cex.axis=1.25, cex.lab=1.5)
abline(v=c(0.7,1.5),col=alpha("grey",1))
lines(x=c(0,0.6),y=rep(0.01,each=2),col=2,lty=3, lwd=2)
lines(x=c(0.8,1.4),y=rep(0.05,each=2),col=2,lty=3, lwd=2)
lines(x=c(1.6,2.2),y=rep(0.1,each=2),col=2,lty=3, lwd=2)



###### also testing t1
resultsMatRegular01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatRegular10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

resultsMatSW01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW01) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW05) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatSW10) <- c("fdrScreenSW", "overallFDRSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW","powerT1SW", "powerT2SW", "nullGeneFDRSW", "nrFalsePositiveNullGenesSW", "totalGenesFoundSW", "propNullInt")

resultsMatJiangDoerge01 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge01) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge05 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge05) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")
resultsMatJiangDoerge10 <- matrix(NA,nrow=N,ncol=13)
colnames(resultsMatJiangDoerge10) <- c("overallFDR", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrInt", "fdrAllHyp", "powerT1", "powerT2", "nullGeneFDR", "nrFalsePositiveNullGenes", "totalGenesFound", "propNullInt")

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

}



