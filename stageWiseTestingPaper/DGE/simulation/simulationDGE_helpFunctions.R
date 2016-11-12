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


	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound)))

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

	return(c(overallFDR=overallFDR, powerInteraction=powerInteraction, fdrInteraction=fdrInt, fdrT1=fdrT1, fdrT2=fdrT2, fdrAllHyp=fdrAllHyp, powerT1=powerT1, powerT2=powerT2,   nullGeneFDR=nullGeneFDR , nrFalsePositiveNullGenes=length(nullGenesFound), totalGenesFound=length(allGenesFound)))
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

return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW)))
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

return(c(fdrScreenSW=fdrScreenSW, overallFDRSW=overallFDRSW, powerInteractionSW=powerInteractionSW, fdrInteractionSW=fdrIntSW, fdrT1SW=fdrT1SW, fdrT2SW=fdrT2SW, fdrAllHypSW=fdrAllHypSW,  powerT1SW=powerT1SW, powerT2SW=powerT2SW, nullGeneFDRSW=nullGeneFDRSW, nrFalsePositiveNullGenesSW=length(nullGenesFoundSW), totalGenesFoundSW=length(allGenesFoundSW)))
}






