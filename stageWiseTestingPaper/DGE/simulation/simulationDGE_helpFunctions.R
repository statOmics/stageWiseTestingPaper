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

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
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


	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))

	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
}



doStageWiseAnalysis <- function(data,L,alpha,design, fit){
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(fit)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	if(length(significantGenesStageI)==0) 	return(c(fdrScreenSW=0, overallFDRSW=0, powerInteractionSW=0, fdrInteractionSW=0, fdrT1SW=0, fdrT2SW=0, fdrAllHypSW=0, powerT1SW=0, powerT2SW=0, nullGeneFDRSW=0,   nrFalsePositiveNullGenesSW=0 , totalGenesFoundSW=0, propNullInt=0))

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

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundIntSW[!(foundIntSW%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt))
}


doStageWiseAnalysisLimma <- function(data,alpha,design, fit){
    significantGenesStageI <- rownames(fit)[p.adjust(fit$F.p.value,"BH")<alpha]
		if(length(significantGenesStageI)==0) 	return(c(fdrScreenSW=0, overallFDRSW=0, powerInteractionSW=0, fdrInteractionSW=0, fdrT1SW=0, fdrT2SW=0, fdrAllHypSW=0, powerT1SW=0, powerT2SW=0, nullGeneFDRSW=0,   nrFalsePositiveNullGenesSW=0 , totalGenesFoundSW=0, propNullInt=0))

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

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundIntSW[!(foundIntSW%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt))
}



compareRegularAndStageWiseAnalysisEdgeR <- function(data,L,alpha,design,fit){
    ### regular analysis
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


    ### stage-wise analysis
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

    ### compare
	uniqueGenesRegular = allGenesFound[!allGenesFound%in%allGenesFoundSW]
	fdrExtraGenesRegular = mean(uniqueGenesRegular%in%nonDeGenes)
	return(fdrExtraGenesRegular)
}

compareRegularAndStageWiseAnalysisLimma <- function(data,L,alpha,design,fit){
    ### regular analysis
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

    ### stage-wise analysis
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

	 ### compare
	uniqueGenesRegular = allGenesFound[!allGenesFound%in%allGenesFoundSW]
	fdrExtraGenesRegular = mean(uniqueGenesRegular%in%nonDeGenes)
	return(fdrExtraGenesRegular)
}

doStageWiseJiangDoergeAnalysis <- function(data,L,alpha1,alpha2,design,fit){
	lrtStep1 <- glmLRT(fit,contrast=L)
	significantGenesStepI <- rownames(fit)[p.adjust(lrtStep1$table$PValue,"BH")<=alpha1]
	if(length(significantGenesStepI)==0) 	return(c(overallFDR=0, powerInteraction=0, fdrInteraction=0, fdrT1=0, fdrT2=0, fdrInt=0, fdrAllHyp=0, powerT1=0, powerT2=0,   nullGeneFDR=0 , nrFalsePositiveNullGenes=0, totalGenesFound=0, fdrAvg=0, powerAvg=0))

	fitStep2 <- fit[significantGenesStepI,]
	lrtListStep2 <- list()
	for(i in 1:ncol(L)) lrtListStep2[[i]] <- glmLRT(fitStep2,contrast=L[,i])
	pvalMatrixStep2 <- do.call(cbind,lapply(lrtListStep2,function(x) x$table$PValue))
	# aggregate p-values and perform FDR on aggregated set from Step 2 genes.
	resultsStep2 = matrix(p.adjust(c(pvalMatrixStep2),"BH")<=alpha2,nrow=nrow(pvalMatrixStep2), ncol=ncol(pvalMatrixStep2), byrow=FALSE)
	dimnames(resultsStep2) <- list(rownames(fitStep2),colnames(L))
	foundT1 <- names(which(resultsStep2[,"t1"]==TRUE))
	foundT2 <- names(which(resultsStep2[,"t2"]==TRUE))
	foundInt <- names(which(resultsStep2[,"interaction"]==TRUE))
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

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
}

doStageWiseJiangDoergeAnalysisLimma <- function(data,L,alpha1,alpha2,design,fit){
	significantGenesStepI <- rownames(fit)[p.adjust(fit$F.p.value,"BH")<=alpha1]
	if(length(significantGenesStepI)==0) 	return(c(overallFDR=0, powerInteraction=0, fdrInteraction=0, fdrT1=0, fdrT2=0, fdrInt=0, fdrAllHyp=0, powerT1=0, powerT2=0,   nullGeneFDR=0 , nrFalsePositiveNullGenes=0, totalGenesFound=0, fdrAvg=0, powerAvg=0))

	fitStep2 <- fit[significantGenesStepI,]
	lrtListStep2 <- list()
	for(i in 1:ncol(L)) lrtListStep2[[i]] <- topTable(fitStep2,coef=i,sort.by="none",number=Inf)
	pvalMatrixStep2 <- do.call(cbind,lapply(lrtListStep2,function(x) x$P.Value))
	if(is.null(pvalMatrixStep2)) 	return(c(overallFDR=0, powerInteraction=0, fdrInteraction=0, fdrT1=0, fdrT2=0, fdrInt=0, fdrAllHyp=0, powerT1=0, powerT2=0,   nullGeneFDR=0 , nrFalsePositiveNullGenes=0, totalGenesFound=0))

	# aggregate p-values and perform FDR on aggregated set from Step 2 genes.
	resultsStep2 = matrix(p.adjust(c(pvalMatrixStep2),"BH")<=alpha2,nrow=nrow(pvalMatrixStep2), ncol=ncol(pvalMatrixStep2), byrow=FALSE)
	dimnames(resultsStep2) <- list(rownames(fitStep2),colnames(L))
	foundT1 <- names(which(resultsStep2[,"t1"]==TRUE))
	foundT2 <- names(which(resultsStep2[,"t2"]==TRUE))
	foundInt <- names(which(resultsStep2[,"interaction"]==TRUE))
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

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))

	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
}


##########################################
##### including average DE contrast ######
##########################################
doRegularAnalysisAvg <- function(data,L,alpha,design, fit){
	lrtList <- list()
	for(i in 1:ncol(L)) lrtList[[i]] <- glmLRT(fit,contrast=L[,i])
	foundT1 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
	foundT2 <- rownames(d)[p.adjust(lrtList[[2]]$table$PValue,"BH")<alpha]
	foundInt <- rownames(d)[p.adjust(lrtList[[3]]$table$PValue,"BH")<alpha]
	foundAvg <- rownames(d)[p.adjust(lrtList[[4]]$table$PValue,"BH")<alpha]
	fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	fdrAvg <- mean(!foundAvg%in%indAvgAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fpAvg <- foundAvg[!foundAvg%in%indAvgAll]
	fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt)+length(fpAvg))/(length(foundT1)+length(foundT2)+length(foundInt)+length(foundAvg))
	overallFDR <- length(unique(c(fpT1,fpT2,fpInt,fpAvg)))/length(unique(c(foundT1,foundT2,foundInt,foundAvg)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt,foundAvg))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)
	powerAvg <- mean(indAvgAll%in%foundAvg)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt, fdrAvg=fdrAvg, powerAvg=powerAvg))
}


doRegularAnalysisLimmaAvg <- function(data,alpha,design, fit){
    	ttList <- list()
	for(i in 1:4) ttList[[i]] <- topTable(fit,coef=i,sort.by="none",number=Inf)
	foundT1 <- rownames(fit)[p.adjust(ttList[[1]]$P.Value,"BH")<alpha]
	foundT2 <- rownames(fit)[p.adjust(ttList[[2]]$P.Value,"BH")<alpha]
	foundInt <- rownames(fit)[p.adjust(ttList[[3]]$P.Value,"BH")<alpha]
	foundAvg <- rownames(d)[p.adjust(lrtList[[4]]$table$PValue,"BH")<alpha]
	fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	fdrAvg <- mean(!foundAvg%in%indAvgAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fpAvg <- foundAvg[!foundAvg%in%indAvgAll]
	fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt)+length(fpAvg))/(length(foundT1)+length(foundT2)+length(foundInt)+length(foundAvg))
	overallFDR <- length(unique(c(fpT1,fpT2,fpInt,fpAvg)))/length(unique(c(foundT1,foundT2,foundInt,foundAvg)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt,foundAvg))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)
	powerAvg <- mean(indAvgAll%in%foundAvg)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt, fdrAvg=fdrAvg, powerAvg=powerAvg))
}


doStageWiseAnalysisAvg <- function(data,L,alpha,design, fit){
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(fit)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	lrtListSW <- list()
	for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
	foundT1SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundT2SW <- rownames(fit2)[lrtListSW[[2]]$table$PValue<alphaAdjusted]
	foundIntSW <- rownames(fit2)[lrtListSW[[3]]$table$PValue<alphaAdjusted]
	foundAvgSW <- rownames(fit2)[lrtListSW[[4]]$table$PValue<alphaAdjusted]
	fdrT1SW <- mean(!foundT1SW%in%indT1All)
	fdrT2SW <- mean(!foundT2SW%in%indT2All)
	fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
	fdrAvgSW <- mean(!foundAvgSW%in%indAvgAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fpAvgSW <- foundAvgSW[!foundAvgSW%in%indAvgAll]
	fdrScreenSW <- mean(significantGenesStageI%in%nonDeGenes)
	fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW)+length(fpAvgSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI)+length(foundAvgSW))
	overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW,fpAvgSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI,foundAvgSW)))
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI,foundAvgSW))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
	powerT1SW <- mean(indT1All%in%foundT1SW)
	powerT2SW <- mean(indT2All%in%foundT2SW)
	powerAvgSW <- mean(indAvgAll%in%foundAvgSW)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundIntSW[!(foundIntSW%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt, fdrAvgSW=fdrAvgSW, powerAvgSW=powerAvgSW))
}


doStageWiseAnalysisLimmaAvg <- function(data,alpha,design, fit){
    significantGenesStageI <- rownames(fit)[p.adjust(fit$F.p.value,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	ttListSW <- list()
	for(i in 1:4) ttListSW[[i]] <- topTable(fit2,coef=i,sort.by="none",number=Inf)
	foundT1SW <- rownames(fit2)[ttListSW[[1]]$P.Value<alphaAdjusted]
	foundT2SW <- rownames(fit2)[ttListSW[[2]]$P.Value<alphaAdjusted]
	foundIntSW <- rownames(fit2)[ttListSW[[3]]$P.Value<alphaAdjusted]
	foundAvgSW <- rownames(fit2)[ttListSW[[4]]$P.Value<alphaAdjusted]
	fdrT1SW <- mean(!foundT1SW%in%indT1All)
	fdrT2SW <- mean(!foundT2SW%in%indT2All)
	fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
	fdrAvgSW <- mean(!foundAvgSW%in%indAvgAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fpAvgSW <- foundAvgSW[!foundAvgSW%in%indAvgAll]
	fdrScreenSW <- mean(significantGenesStageI%in%nonDeGenes)
	fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW)+length(fpAvgSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI)+length(foundAvgSW))
	overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW,fpAvgSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI,foundAvgSW)))
	allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI,foundAvgSW))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
	powerT1SW <- mean(indT1All%in%foundT1SW)
	powerT2SW <- mean(indT2All%in%foundT2SW)
	powerAvgSW <- mean(indAvgAll%in%foundAvgSW)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundIntSW[!(foundIntSW%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt, fdrAvgSW=fdrAvgSW, powerAvgSW=powerAvgSW))
}

doStageWiseJiangDoergeAnalysisAvg <- function(data,L,alpha1,alpha2,design,fit){
	lrtStep1 <- glmLRT(fit,contrast=L)
	significantGenesStepI <- rownames(fit)[p.adjust(lrtStep1$table$PValue,"BH")<=alpha1]
	fitStep2 <- fit[significantGenesStepI,]
	lrtListStep2 <- list()
	for(i in 1:ncol(L)) lrtListStep2[[i]] <- glmLRT(fitStep2,contrast=L[,i])
	pvalMatrixStep2 <- do.call(cbind,lapply(lrtListStep2,function(x) x$table$PValue))
	# aggregate p-values and perform FDR on aggregated set from Step 2 genes.
	resultsStep2 = matrix(p.adjust(c(pvalMatrixStep2),"BH")<=alpha2,nrow=nrow(pvalMatrixStep2), ncol=ncol(pvalMatrixStep2), byrow=FALSE)
	dimnames(resultsStep2) <- list(rownames(fitStep2),colnames(L))
	foundT1 <- names(which(resultsStep2[,"t1"]==TRUE))
	foundT2 <- names(which(resultsStep2[,"t2"]==TRUE))
	foundInt <- names(which(resultsStep2[,"interaction"]==TRUE))
	foundAvg <- names(which(resultsStep2[,"avgDE"]==TRUE))
	fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	fdrAvg <- mean(!foundAvg%in%indAvgAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fpAvg <- foundAvg[!foundAvg%in%indAvgAll]
	fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt)+length(fpAvg))/(length(foundT1)+length(foundT2)+length(foundInt)+length(foundAvg))
	overallFDR <- length(unique(c(fpT1,fpT2,fpInt,fpAvg)))/length(unique(c(foundT1,foundT2,foundInt,foundAvg)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt,foundAvg))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)
	powerAvg <- mean(indAvgAll%in%foundAvg)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt, fdrAvg=fdrAvg, powerAvg=powerAvg))
}

doStageWiseJiangDoergeAnalysisLimmaAvg <- function(data,L,alpha1,alpha2,design,fit){
	significantGenesStepI <- rownames(fit)[p.adjust(fit$F.p.value,"BH")<=alpha1]
	fitStep2 <- fit[significantGenesStepI,]
	lrtListStep2 <- list()
	for(i in 1:ncol(L)) lrtListStep2[[i]] <- topTable(fitStep2,coef=i,sort.by="none",number=Inf)
	pvalMatrixStep2 <- do.call(cbind,lapply(lrtListStep2,function(x) x$P.Value))
	if(is.null(pvalMatrixStep2)) 	return(c(overallFDR=0, powerInteraction=0, fdrInteraction=fdrInt, fdrT1=0, fdrT2=0, fdrInt=0, fdrAllHyp=0, powerT1=0, powerT2=0,   nullGeneFDR=0 , nrFalsePositiveNullGenes=0, totalGenesFound=0, fdrAvg=0, powerAvg=0))

	# aggregate p-values and perform FDR on aggregated set from Step 2 genes.
	resultsStep2 = matrix(p.adjust(c(pvalMatrixStep2),"BH")<=alpha2,nrow=nrow(pvalMatrixStep2), ncol=ncol(pvalMatrixStep2), byrow=FALSE)
	dimnames(resultsStep2) <- list(rownames(fitStep2),colnames(L))
	foundT1 <- names(which(resultsStep2[,"t1"]==TRUE))
	foundT2 <- names(which(resultsStep2[,"t2"]==TRUE))
	foundInt <- names(which(resultsStep2[,"interaction"]==TRUE))
	foundAvg <- names(which(resultsStep2[,"avgDE"]==TRUE))
	fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	fdrAvg <- mean(!foundAvg%in%indAvgAll)
	fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	fpAvg <- foundAvg[!foundAvg%in%indAvgAll]
	fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt)+length(fpAvg))/(length(foundT1)+length(foundT2)+length(foundInt)+length(foundAvg))
	overallFDR <- length(unique(c(fpT1,fpT2,fpInt,fpAvg)))/length(unique(c(foundT1,foundT2,foundInt,foundAvg)))
	allGenesFound <- unique(c(foundT1,foundT2,foundInt,foundAvg))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)
	powerAvg <- mean(indAvgAll%in%foundAvg)


	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))

	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt, fdrAvg=fdrAvg, powerAvg=powerAvg))
}



####################################################################
#### EFFECT OF OTHER HYPOTHESES ON JIANG METHOD PERFORMANCE ########
####################################################################

doStageWiseAnalysis_t2Int <- function(data,L,alpha,design, fit){
	lrt1 <- glmLRT(fit,contrast=L)
	significantGenesStageI <- rownames(fit)[p.adjust(lrt1$table$PValue,"BH")<alpha]
	fit2 <- fit[significantGenesStageI,]
	alphaAdjusted <- alpha*length(significantGenesStageI)/nrow(fit)
	lrtListSW <- list()
	for(i in 1:ncol(L)) lrtListSW[[i]] <- glmLRT(fit2,contrast=L[,i])
	#foundT1SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundT2SW <- rownames(fit2)[lrtListSW[[1]]$table$PValue<alphaAdjusted]
	foundIntSW <- rownames(fit2)[lrtListSW[[2]]$table$PValue<alphaAdjusted]
	#fdrT1SW <- mean(!foundT1SW%in%indT1All)
	fdrT2SW <- mean(!foundT2SW%in%indT2All)
	fdrIntSW <- mean(!foundIntSW%in%indInteractionAll)
	fpScreenSW <- significantGenesStageI[!significantGenesStageI%in%indAll]
	#fpT1SW <- foundT1SW[!foundT1SW%in%indT1All]
	fpT2SW <- foundT2SW[!foundT2SW%in%indT2All]
	fpIntSW <- foundIntSW[!foundIntSW%in%indInteractionAll]
	fdrScreenSW <- mean(significantGenesStageI%in%nonDeGenes)
	#fdrAllHypSW <- (length(fpT1SW)+length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT1SW)+length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	fdrAllHypSW <- (length(fpT2SW)+length(fpIntSW)+length(fpScreenSW))/(length(foundT2SW)+length(foundIntSW)+length(significantGenesStageI))
	#overallFDRSW <- length(unique(c(fpT1SW,fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI)))
	overallFDRSW <- length(unique(c(fpT2SW,fpIntSW,fpScreenSW)))/length(unique(c(foundT2SW,foundIntSW,significantGenesStageI)))
	#allGenesFoundSW <- unique(c(foundT1SW,foundT2SW,foundIntSW,significantGenesStageI))
	allGenesFoundSW <- unique(c(foundT2SW,foundIntSW,significantGenesStageI))
	nullGenesFoundSW <- allGenesFoundSW[allGenesFoundSW%in%nonDeGenes]
	nullGeneFDRSW <- length(nullGenesFoundSW)/length(allGenesFoundSW)

	powerInteractionSW <- mean(indInteractionAll%in%foundIntSW)
	#powerT1SW <- mean(indT1All%in%foundT1SW)
	powerT2SW <- mean(indT2All%in%foundT2SW)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundIntSW[!(foundIntSW%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


#return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt))
return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW), propNullInt=propNullInt))

}



doRegularAnalysis_t2Int <- function(data,L,alpha,design, fit){
	lrtList <- list()
	for(i in 1:ncol(L)) lrtList[[i]] <- glmLRT(fit,contrast=L[,i])
	#foundT1 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
	foundT2 <- rownames(d)[p.adjust(lrtList[[1]]$table$PValue,"BH")<alpha]
	foundInt <- rownames(d)[p.adjust(lrtList[[2]]$table$PValue,"BH")<alpha]
	#fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	#fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	#fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
	fdrAllHyp <- (length(fpT2)+length(fpInt))/(length(foundT2)+length(foundInt))
	#overallFDR <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
	overallFDR <- length(unique(c(fpT2,fpInt)))/length(unique(c(foundT2,foundInt)))
	#allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	allGenesFound <- unique(c(foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	#powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
}




doStageWiseJiangDoergeAnalysis_t2Int <- function(data,L,alpha1,alpha2,design,fit){
	lrtStep1 <- glmLRT(fit,contrast=L)
	significantGenesStepI <- rownames(fit)[p.adjust(lrtStep1$table$PValue,"BH")<=alpha1]
	fitStep2 <- fit[significantGenesStepI,]
	lrtListStep2 <- list()
	for(i in 1:ncol(L)) lrtListStep2[[i]] <- glmLRT(fitStep2,contrast=L[,i])
	pvalMatrixStep2 <- do.call(cbind,lapply(lrtListStep2,function(x) x$table$PValue))
	# aggregate p-values and perform FDR on aggregated set from Step 2 genes.
	resultsStep2 = matrix(p.adjust(c(pvalMatrixStep2),"BH")<=alpha2,nrow=nrow(pvalMatrixStep2), ncol=ncol(pvalMatrixStep2), byrow=FALSE)
	dimnames(resultsStep2) <- list(rownames(fitStep2),colnames(L))
	#foundT1 <- names(which(resultsStep2[,"t1"]==TRUE))
	foundT2 <- names(which(resultsStep2[,"t2"]==TRUE))
	foundInt <- names(which(resultsStep2[,"interaction"]==TRUE))
	#fdrT1 <- mean(!foundT1%in%indT1All)
	fdrT2 <- mean(!foundT2%in%indT2All)
	fdrInt <- mean(!foundInt%in%indInteractionAll)
	#fpT1 <- foundT1[!foundT1%in%indT1All]
	fpT2 <- foundT2[!foundT2%in%indT2All]
	fpInt <- foundInt[!foundInt%in%indInteractionAll]
	#fdrAllHyp <- (length(fpT1)+length(fpT2)+length(fpInt))/(length(foundT1)+length(foundT2)+length(foundInt))
	fdrAllHyp <- (length(fpT2)+length(fpInt))/(length(foundT2)+length(foundInt))
	#overallFDR <- length(unique(c(fpT1,fpT2,fpInt)))/length(unique(c(foundT1,foundT2,foundInt)))
	overallFDR <- length(unique(c(fpT2,fpInt)))/length(unique(c(foundT2,foundInt)))
	#allGenesFound <- unique(c(foundT1,foundT2,foundInt))
	allGenesFound <- unique(c(foundT2,foundInt))
	nullGenesFound <- allGenesFound[allGenesFound%in%nonDeGenes]
	nullGeneFDR <- length(nullGenesFound)/length(allGenesFound)

	powerInteraction <- mean(indInteractionAll%in%foundInt)
	#powerT1 <- mean(indT1All%in%foundT1)
	powerT2 <- mean(indT2All%in%foundT2)

	## how many of the false positive interaction genes have main effects
	falseIntGenes = foundInt[!(foundInt%in%indInteractionAll)]
	propNullInt = 1-mean(falseIntGenes%in%c(indT1All,indT2All))


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT2=fdrT2, fdrInt=fdrInt, fdrAllHyp=fdrAllHyp, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound), propNullInt=propNullInt))
}
