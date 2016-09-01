doRegularAnalysis <- function(data,L,alpha,design, fit){
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
	allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)
	
	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound)))

}

doRegularAnalysisLimma <- function(data,alpha,design, fit){
    	ttList <- list()
	for(i in 1:3) ttList[[i]] <- topTable(fit,coef=i,sort.by="none",number=Inf)
	foundT1 <- rownames(fit)[p.adjust(ttList[[1]]$P.Value,"BH")<alpha]
	foundT2 <- rownames(fit)[p.adjust(ttList[[2]]$P.Value,"BH")<alpha]
	foundInt <- rownames(fit)[p.adjust(ttList[[3]]$P.Value,"BH")<alpha]
	fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
	overallFDR <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)
	
	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)

	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound)))
}



doStageWiseAnalysis <- function(data,L,alpha,design, fit){
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(fit)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	lrtListSW <- list()
	for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
	foundT1SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundT2SW <- rownames(fit2)[lrtListSW[[2]]$table$PValue<alphaAdjusted]
	foundIntSW <- rownames(fit2)[lrtListSW[[3]]$table$PValue<alphaAdjusted]
	fdrT1SW <- mean(!foundT1SW%in%indT1All)
	fdrT2SW <- mean(!foundT2SW%in%indT2All)
	fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fdrScreenSW <- mean(significantGenesStageI%in%nonDeGenes)
	fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI)))
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
	powerT1SW <- mean(indT1All%in%foundT1SW)
	powerT2SW <- mean(indT2All%in%foundT2SW)

return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW)))
}


doStageWiseAnalysisLimma <- function(data,alpha,design, fit){
    	significantGenesStageI <- rownames(fit)[p.adjust(fit$F.p.value,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	ttListSW <- list()
	for(i in 1:3) ttListSW[[i]] <- topTable(fit2,coef=i,sort.by="none",number=Inf)
	foundT1SW <- rownames(fit2)[ttListSW[[1]]$P.Value<alphaAdjusted]
	foundT2SW <- rownames(fit2)[ttListSW[[2]]$P.Value<alphaAdjusted]
	foundIntSW <- rownames(fit2)[ttListSW[[3]]$P.Value<alphaAdjusted]
	fdrT1SW <- mean(!foundT1SW%in%indT1All)
	fdrT2SW <- mean(!foundT2SW%in%indT2All)
	fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fdrScreenSW <- mean(significantGenesStageI%in%nonDeGenes)
	fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI)))
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
	powerT1SW <- mean(indT1All%in%foundT1SW)
	powerT2SW <- mean(indT2All%in%foundT2SW)

return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW)))
}





##############
doSimulation <- function(nCst,nInt,N, nreps){
resultsMatRegular01 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatRegular01) <- c("fdrGene", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "powerGene", "nrFalsePositives", "fractionFalsePositiveNullGenes", "nrFalsePositiveNullGenes")
resultsMatRegular05 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatRegular05) <- c("fdrGene", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "powerGene", "nrFalsePositives", "fractionFalsePositiveNullGenes", "nrFalsePositiveNullGenes")
resultsMatRegular10 <- matrix(NA,nrow=N,ncol=12)
colnames(resultsMatRegular10) <- c("fdrGene", "powerInteraction", "fdrInteraction", "fdrT1", "fdrT2", "fdrAllHyp", "powerT1", "powerT2", "powerGene", "nrFalsePositives", "fractionFalsePositiveNullGenes", "nrFalsePositiveNullGenes")



resultsMatSW01 <- matrix(NA,nrow=N,ncol=14)
colnames(resultsMatSW01) <- c("fdrScreen", "fdrGeneSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW", "fracPositiveNullGenesSW", "powerT1SW", "powerT2SW", "powerGeneSW", "nrFalsePositivesSW", "fractionFalsePositiveNullGenesSW", "nrFalsePositiveNullGenesSW")
resultsMatSW05 <- matrix(NA,nrow=N,ncol=14)
colnames(resultsMatSW05) <- c("fdrScreen", "fdrGeneSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW", "fracPositiveNullGenesSW", "powerT1SW", "powerT2SW", "powerGeneSW", "nrFalsePositivesSW", "fractionFalsePositiveNullGenesSW", "nrFalsePositiveNullGenesSW")
resultsMatSW10 <- matrix(NA,nrow=N,ncol=14)
colnames(resultsMatSW10) <- c("fdrScreen", "fdrGeneSW", "powerInteractionSW", "fdrInteractionSW", "fdrT1SW", "fdrT2SW", "fdrAllHypSW", "fracPositiveNullGenesSW", "powerT1SW", "powerT2SW", "powerGeneSW", "nrFalsePositivesSW", "fractionFalsePositiveNullGenesSW", "nrFalsePositiveNullGenesSW")

    for(iter in 1:N){
set.seed(iter)
libSize <- runif(nSamp,25e6,30e6)
grp <- as.factor(rep(0:1, each = nreps))


# sample fold changes
fcT1Cst <- 2^rnorm(nCst, mean=sample(empiricalFC, size = nCst, replace = TRUE), sd=densFit$bw)
while(all(abs(log2(fcT1Cst))>1)==FALSE){
    nLow <- sum(abs(log2(fcT1Cst))<1)
    fcT1Cst[abs(log2(fcT1Cst))<1] <- 2^rnorm(nLow, mean=sample(empiricalFC, size = nLow, replace = TRUE), sd=densFit$bw)
    #print(nLow)
}
fcT2Cst <- fcT1Cst
#interaction: sample T1 FC and multiply by sampled interaction FC.
fcInt <- 2^rnorm(nInt, mean=sample(empiricalFCInt, size = nInt, replace = TRUE), sd=densInt$bw)
while(all(abs(log2(fcInt))>1)==FALSE){
    nLow <- sum(abs(log2(fcInt))<1)
    fcInt[abs(log2(fcInt))<1] <- 2^rnorm(nLow, mean=sample(empiricalFCInt, size = nLow, replace = TRUE), sd=densInt$bw)
   # print(nLow)
}
fcIntT1 <- sample(fcT1Cst,nInt,replace=TRUE)
fcIntT2 <- fcIntT1*fcInt

# build DE gene vectors
deT1Cst <- sample(1:nTags,nCst)
remaining1 <- (1:nTags)[-deT1Cst]
deT2Cst <- deT1Cst
deInt <- sample(remaining1,nInt)
remaining4 <- remaining1[!remaining1%in%deInt] #non DE genes

# aggregate fold changes for simulation
fcSim1 <- c(fcT1Cst,fcIntT1)
fcSim2 <- c(fcT2Cst,fcIntT2)


# aggregate indices of DE genes for simulation
DEind1 <- c(deT1Cst,deInt) #constant over both timepoints, unique in timepoint, interaction.
DEind2 <- c(deT2Cst,deInt)



# aggregate indices for performance calculations
deT1All <- c(deT1Cst,deInt)
deT2All <- c(deT2Cst,deInt)
deInt <- c(deInt) #unique in one timepoint is also interaction


dataTime1 <- NBsimSeed(foldDiff = fcSim1, ind=DEind1, dataset = exprs(pickrell.eset), nTags = nTags, group = grp, verbose = FALSE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE, seed=iter)
dataTime2 <- NBsimSeed(foldDiff = fcSim2, ind=DEind2, dataset = exprs(pickrell.eset), nTags = nTags, group = grp, verbose = FALSE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE, seed=iter)
head(cbind(dataTime1$AveLogCPM,dataTime2$AveLogCPM)) #genes are linked across timepoints.


## design
dataAll <- cbind(dataTime1$counts,dataTime2$counts)
treat=factor(rep(c("control","SNL","control","SNL"),each=nreps), levels=c("control","SNL"))
time = factor(rep(c("mo2","w2"),each=nreps*2),levels=c("w2","mo2"))
design <- model.matrix(~treat*time)
L=matrix(0,nrow=ncol(design),ncol=3)
rownames(L)=colnames(design)
colnames(L)=c("t1","t2","interaction")
L[2,1]=1
L[c(2,4),2]=1
L[4,3]=1

## fit models
dSim <- DGEList(dataAll)
dSim <- calcNormFactors(dSim)
dSim <- estimateDisp(dSim,design)
fit <- glmFit(dSim,design)


## regular analysis
resultsMatRegular01[iter,] <- doRegularAnalysis(data=dataAll,L=L,alpha=0.01,design=design,fit=fit)
resultsMatRegular05[iter,] <- doRegularAnalysis(data=dataAll,L=L,alpha=0.05,design=design,fit=fit)
resultsMatRegular10[iter,] <- doRegularAnalysis(data=dataAll,L=L,alpha=0.1,design=design,fit=fit)

## stage-wise analysis
resultsMatSW01[iter,] <- doStageWiseAnalysis(data=dataAll,L=L,alpha=.01,design=design,fit=fit)
resultsMatSW05[iter,] <- doStageWiseAnalysis(data=dataAll,L=L,alpha=.05,design=design,fit=fit)
resultsMatSW10[iter,] <- doStageWiseAnalysis(data=dataAll,L=L,alpha=.1,design=design,fit=fit)
}

print(
boxplot(cbind(resultsMatRegular01[,"fdrAllHyp"],resultsMatSW01[,"fdrAllHypSW"],
	      resultsMatRegular05[,"fdrAllHyp"],resultsMatSW05[,"fdrAllHypSW"],
	      resultsMatRegular10[,"fdrAllHyp"],resultsMatSW10[,"fdrAllHypSW"],
	      resultsMatRegular01[,"fdrGene"],resultsMatSW01[,"fdrGeneSW"],
	      resultsMatRegular05[,"fdrGene"],resultsMatSW05[,"fdrGeneSW"],
	      resultsMatRegular10[,"fdrGene"],resultsMatSW10[,"fdrGeneSW"],
	      resultsMatRegular01[,"fractionFalsePositiveNullGenes"],resultsMatSW01[,"fractionFalsePositiveNullGenesSW"],
	      resultsMatRegular05[,"fractionFalsePositiveNullGenes"],resultsMatSW05[,"fractionFalsePositiveNullGenesSW"],
	      resultsMatRegular10[,"fractionFalsePositiveNullGenes"],resultsMatSW10[,"fractionFalsePositiveNullGenesSW"]),
	boxwex=.2,at=rep(seq(0.3,6.7,by=0.8),each=2)+rep(c(-.1,.1),9), border=rep(c("black","steelblue"),9), col=alpha(rep(c("black","steelblue"),9),.2), xaxt="n", yaxt="n", ylab="False Discovery Rate")
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
text(x=0.2,y=0.13,"All Hypotheses")
text(x=3,y=0.13,"Gene level")
text(x=5.6,y=0.13,"Null gene level")
)

}


