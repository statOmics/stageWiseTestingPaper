
stageWiseTest <- function(pScreen, pConfirmation, alpha, method=c("none","holm","dtu","user"), onlySignificantGenes=TRUE){
    # pList is a matrix of p-values: one row corresponds to one feature. The different columns correspond to different hypotheses where the first column represents the screeening hypothesis and the subsequent columns the specific hypotheses of interest.
    # alpha is the OFDR target level. Can be a vector of length 1 or more.
    # method is either a character vector "Holm" or a numeric vector with length equal to the number of columns in pList - 1, corresponding to all the specific hypotheses of interest.
    # onlySignificantGenes is a logical indiacting whether only the significant genes should be returned. If set to FALSE, the adjusted significance level will be included in the output and p-values from the confirmation stage should be compared to the highest adjusted significance level of the genes called significant.
    method <- match.arg(method,c("none","holm","dtu","user"))
    if(onlySignificantGenes){
	padjScreen <- p.adjust(pScreen,"BH")
    	genesStageI <- padjScreen<alpha
	padjScreen <- padjScreen[genesStageI]
	geneOrder <- order(padjScreen)
	pConfSignificant <- pConfirmation[genesStageI,]	
	pConfSignificant <- pConfSignificant[geneOrder,]
	alphaAdjusted <- sum(genesStageI,na.rm=TRUE)/length(genesStageI)*alpha
	## FWER
	if(method[1]=="holm"){
	    pAdjListSignificant <- t(apply(pConfSignificant,1, function(row){
		o <- order(row)
		n <- length(row)
		adjustment <- c(n-1,(n-1):1)
		rowAdjusted <- row[o]*adjustment
		rowAdjusted <- pmin(rowAdjusted,1)
		rowAdjusted <- cummax(rowAdjusted)
		rowBack[] <- NA
		rowBack[o] <- rowAdjusted
		rowBack
}))
    } else if(is.numeric(method) && length(method)==ncol(pList)-1){
	pAdjListSignificant <- t(apply(pConfSignificant,1,function(row){
		  o <- order(row)
		  rowAdjusted <- row[o]*method
		  rowAdjusted <- pmin(rowAdjusted,1)
		  # check monotone increase of adjusted p-values
		  rowAdjusted <- cummax(rowAdjusted)
		  rowBack[] <- NA
		  rowBack[o] <- rowAdjusted
		  rowBack
	}))
    } else stop("method must be either 'holm' or a numeric vector with length equal to the number of specific hypotheses")
	## times G divided by R
	pAdjListSignificant <- pAdjListSignificant*length(genesStageI)/sum(genesStageI,na.rm=TRUE)
	pAdjListSignificant[pAdjListSignificant>1] <- 1
	return(cbind(padjScreen[geneOrder],pAdjListSignificant))
    } else {
    padjScreen <- p.adjust(pScreen[,1],"BH")
    geneOrder <- order(padjScreen)
    genesStageI <- padjScreen<alpha
    if(method[1]=="holm"){
	pAdjList <- t(apply(pConfirmation,1,function(x) p.adjust(x,method="holm")))
	pAdjList[!genesStageI,] <- 1
    } else if(is.numeric(method) && length(method)==ncol(pList)-1){
	pAdjList <- t(apply(pConfirmation,1,function(row){
		  o <- order(row)
		  rowAdjusted <- row[o]*method
		  rowAdjusted <- pmin(rowAdjusted,1)
		  # check monotone increase of adjusted p-values
		  rowAdjusted <- cummax(rowAdjusted)		  
		  #rowBack[o] <- rowAdjusted[o]
		  rowBack[] <- NA
		  rowBack[o] <- rowAdjusted
		  rowBack
	}))
	pAdjList[!genesStageI,] <- 1
    } else stop("method must be either 'holm' or a numeric vector with length equal to the number of specific hypotheses")
    padjScreenOrdered <- padjScreen[geneOrder]
    alphaAdjusted <- (1:length(padjScreenOrdered))/length(padjScreenOrdered)*alpha
    pListAll <- cbind(padjScreenOrdered,alphaAdjusted,pAdjList[geneOrder,])
    return(pListAll)
    }
}




