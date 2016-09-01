
## prepare for checking information ##
ip <- installed.packages()
Rversions <- "R-base"
biocVersions <- "bioc2.13"
bioc <- "http://www.bioconductor.org/packages/2.13/bioc/"
sMethods <- c("edgeR", "edgeR_robust", "edgeR_rdev", "edgeR_rans")
sVersions <- c("3.5.15", "3.5.15", "3.5.15", "3.5.15")
sSources <- c("bioc2.14", "bioc2.14", "bioc2.14", "bioc2.14")
rMethods <-  c("samr_SAMseq", "ShrinkBayes")
rVersions <- c("2.0", "2.6")
library(edgeR)


getDataset <- function(counts, drop.extreme.dispersion = 0.1, drop.low.lambda = TRUE) {  
## this function generates NB parameters from real dataset ##
## it is low-level function of NBsim ##
	d <- DGEList(counts)
	d <- calcNormFactors(d)
	cp <- round(cpm(d,normalized.lib.sizes=TRUE),1)
	if(drop.low.lambda) d <- d[rowSums(cp>1) >= 2, ]
	d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))
	d <- estimateGLMCommonDisp(d)
	d <- estimateGLMTrendedDisp(d)
	d <- estimateGLMTagwiseDisp(d)
	dispersion <- d$tagwise.dispersion
	AveLogCPM <- d$AveLogCPM
	if(is.numeric(drop.extreme.dispersion))
	{   
		bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
		ids <- dispersion <= bad
		AveLogCPM <- AveLogCPM[ids]
		dispersion <- dispersion[ids]
	}
	dataset.AveLogCPM <- AveLogCPM
	dataset.dispersion <- dispersion
	dataset.lib.size <- d$samples$lib.size
	dataset.nTags <- nrow(d)
	list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags)
}

setClass("FoldList", representation("list"))
setIs("FoldList", "LargeDataObject")
## FoldList is similar to DGEList ##

NBsim <-
function(dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 10, pDiff=.1, pUp=.5, foldDiff=3, name = NULL, save.file = FALSE, file = NULL, only.add.outlier = FALSE, verbose=TRUE, ind=NULL, seed=NULL)

{   
## NBsim generate simulated count from the real dataset followed by the NB model ##		
	require(edgeR)
	group = as.factor(group)

	sample.fun <- function(object)
	{
		## it is low-level function of NBsim ##
        ## it samples from the real dataset ## 
		nlibs <- object$nlibs
		nTags <- object$nTags
		AveLogCPM <-object$dataset$dataset.AveLogCPM
		dispersion <- object$dataset$dataset.dispersion
		#### KVDB
		set.seed(seed)
		#### KVDB
        id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
		Lambda <- 2^(AveLogCPM[id_r])
		Lambda <- Lambda/sum(Lambda)
		Dispersion <- dispersion[id_r]
		id_0<- Lambda == 0
		Lambda <- Lambda[!id_0]
		Dispersion <- Dispersion[!id_0]
		Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
		object$Lambda <- Lambda
		if(!is.na(fix.dispersion))
		Dispersion <- expandAsMatrix(fix.dispersion, dim = c(nTags, nlibs))
		else Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
		object$Dispersion <- Dispersion
		object
		
	}
	diff.fun <- function(object)
	{ 
		
        ## it is low-level function of NBsim ##
        ## it creates diff genes according to foldDiff ## 
		group <- object$group
		pDiff <- object$pDiff
		pUp <-  object$pUp 
		foldDiff <- object$foldDiff
		Lambda <- object$Lambda
		nTags <- object$nTags
		g <- group == levels(group)[1]
		#### KVDB
		#ind <- sample(nTags, floor(pDiff*nTags))
		ind <- ind
		#### KVDB
		if(length(ind)>0 & !all(foldDiff == 1) ) {
			fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
			Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
			Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir)) 
            #Lambda <- t(t(Lambda)/colSums(Lambda))
			object$Lambda <- Lambda
			object$indDE <- ind
			object$indNonDE <- (1:nTags)[-ind]
			object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda))
			object$mask_DEup[ind[fcDir == 1], g] <- TRUE
			object$mask_DEup[ind[fcDir == -1], !g] <- TRUE
			object$mask_DEdown[ind[fcDir == -1], g] <- TRUE
			object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE
			object$mask_nonDE[-ind,] <- TRUE
			object$mask_DE <- object$mask_DEup | object$mask_DEdown}
		if(all(foldDiff == 1)| pDiff == 0)
		object$indDE <- NA
		object
	}
	sim.fun <- function(object)
	{   
        ## it is low-level function of NBsim ##
        ## it simulate counts using rnbinom ## 
		Lambda <- object$Lambda
		Dispersion <- object$Dispersion
		nTags <- object$nTags
		nlibs <- object$nlibs
		lib.size <- object$lib.size
		counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs) 
		rownames(counts) <- paste("ids", 1:nTags, sep = "")
		object$counts <- counts
		object
	}
			
	outlier.fun <- function(object, outlierMech, pOutlier, min.factor = 2, max.factor = 5)
        {   
        ## it is low-level function of NBsim ##
        ## it makes outlier ## 
	        outlierMech <- match.arg(outlierMech, c("S", "M", "R"))
	        dim <- dim(object$counts)
                outlier.factor <- function() runif(1, min.factor, max.factor)
                countAddOut <- object$counts
                LambdaAddOut <- object$Lambda
	        DispersionAddOut <- object$Dispersion	
	        switch(outlierMech, 
	               S = {
	                    mask_outlier <- expandAsMatrix(FALSE, dim = dim)
	                    id_r <- which(runif(dim[1]) < pOutlier)
	                    id_c <- sample(dim[2], length(id_r), replace = TRUE)
		            for(i in seq(id_r)) 
		                   mask_outlier[id_r[i], id_c[i]] <- TRUE
                            countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
                            }, 
	               R = {				   
                            mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
                            countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
                            },
         
	               M = {
                            mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
                            LambdaAddOut[mask_outlier] <- sapply(LambdaAddOut[mask_outlier], function(z) z*outlier.factor())
	                    countAddOut[mask_outlier] <- rnbinom(sum(mask_outlier), mu = t(t(LambdaAddOut)*object$lib.size)[mask_outlier], size = 1/DispersionAddOut[mask_outlier])
                           }
                       )
                 if(!object$foldDiff == 1 & !pDiff == 0)
	         {
	                indDEupOutlier <- which(apply(object$mask_DEup & mask_outlier, 1, any))
	                indDEdownOutlier <- which(apply(object$mask_DEdown & mask_outlier, 1, any))
	                indDEnoOutlier <- which(apply((object$mask_DE & !mask_outlier) , 1, all))
	                indNonDEOutlier <- which(apply(object$mask_nonDE & mask_outlier, 1, any))
	                indNonDEnoOutlier <- which(apply((object$mask_nonDE & !mask_outlier) , 1, all))
	                indDEbothOutlier <- NA
	                o <- indDEupOutlier %in% indDEdownOutlier
	                q <-  indDEdownOutlier %in% indDEupOutlier
	                if(any(o))
	                {
                              indDEupOutlier <- indDEupOutlier[!o]
		              indDEbothOutlier <- indDEupOutlier[o]	
	                      indDEdownOutlier <- indDEdownOutlier[!q]	
	                }	
	         }
	         else
	         {
                        indDEupOutlier <- indDEdownOutlier <- indDEnoOutlier <- indNonDEOutlier <- indNonDEnoOutlier <- indDEbothOutlier <- NA
                 }
	             out <- list(countAddOut = countAddOut, outlierMech = outlierMech, pOutlier = pOutlier, mask_outlier = mask_outlier, indDEupOutlier = indDEupOutlier, 
                                 indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, 
                                 indNonDEnoOutlier = indNonDEnoOutlier, LambdaAddOut = LambdaAddOut, DispersionAddOut = DispersionAddOut) 

        }
	
	calProb <- function(x, l) round(1 -(1 - x)^(1/l), 2) ## calculate probability to make sure all the outlierMech produce the same amount of outliers ##


        if(verbose) message("Preparing dataset.\n")	
	if(class(dataset) == "DGEList")
	{   
		dat <- dataset
		dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
	}
	else if(is.character(dataset)) 
	{
		load(dataset)
		dat <- get(gsub("(\\.)(\\w+)", "", basename(dataset)))
		dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
	}
	else if(is.matrix(dataset)) 
	{ 
	  if(is.null(name)) name <- deparse(substitute(dataset))	
	  dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
	  dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))		
	}
	else
	dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))

	if(!only.add.outlier)
	{
		if(is.null(lib.size))
    	  dat$lib.size <- runif(nlibs, min = 0.7*median(dat$dataset$dataset.lib.size), max = 1.3*median(dat$dataset$dataset.lib.size))

	    if(is.null(nTags))
	      dat$nTags <- dat$dataset$dataset.nTags 
        if(verbose) message("Sampling.\n")	
	      dat <- sample.fun(dat)
        if(verbose) message("Calculating differential expression.\n")	
	      dat <- diff.fun(dat)
        if(verbose) message("Simulating data.\n")	
	      dat <- sim.fun(dat)
	}
	if(add.outlier){
		outlierMech <- match.arg(outlierMech,  c("S", "R", "M"), several.ok = TRUE)
		if(length(pOutlier)== 1 & length(outlierMech) > 1 & any(outlierMech == "S"))
		{ 
			prob <- calProb(pOutlier, length(group))	
			pOutlier <- rep(pOutlier, length = length(outlierMech))
			pOutlier[!outlierMech == "S"] <- prob	
		}
		else if(!length(pOutlier) == length(outlierMech))
		stop("pOutlier is not equal to outlierMech")
                if(verbose) message("Adding outliers.\n")	
		dat[outlierMech] <- mapply(function(x, y) outlier.fun(dat, outlierMech = x, pOutlier = y, min.factor = min.factor, max.factor = max.factor), x = outlierMech, y = pOutlier, SIMPLIFY = FALSE)
	    dat$pOutlier <- pOutlier
	}

	if(save.file)
	{ 
		
		## save file for shiny app ##
		if(verbose) message("Saving file.\n")
		if(is.null(file)) 
		{ g <- paste0("g", sum(levels(group)[1] == group), "v", sum(levels(group)[2] == group))
			f <- paste0("f", foldDiff)
			if(add.outlier) o <- paste0("o", sprintf( "%02d",100*pOutlier[1L])) 
			else o <- paste0("o", sprintf( "%02d", 0 ))  
			file <- paste0(dat$name, "/", g, f, o, ".Rdata")
			dir.create(dat$name, showWarnings = FALSE)  
		}
		filenm <- eval(gsub("(\\.)(\\w+)", "", basename(file)))
		assign(filenm, dat)
		save(list = filenm, file = file)		
	}
	dat 	
}

NBsimFold <-
function(fold_seq = seq(1, 3, by = 0.25), dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 8, pDiff=.1, pUp=.5, name = NULL)
{  
    ## NBsimFold generates a series of simulated counts followed by the NB model ##
	o <- order(fold_seq)
	fold_seq <- fold_seq[o]
	fold_seq <- unique(round(fold_seq, 4))
	if(is.null(name)) name <- deparse(substitute(dataset))
	if(is.matrix(dataset)) dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion)
	if(all(fold_seq[1] == fold_seq) ) names(fold_seq) <- paste("sample_", seq(fold_seq), sep = "" )
	else names(fold_seq) <- paste("fold_", fold_seq, sep = "" )
	re_split <- function(object)
    {out <- list()
		for (what in names(object[[1L]])) {
			out[[what]] <- lapply(object, .subset2, what)
			if(identical(out[[what]][[1L]], out[[what]][[2L]]))
			out[[what]] <- out[[what]][[1L]]}
		out}
	
	rere_split <- function(object, x)
	{ object[[x]] <- re_split(object[[x]])
	  object}
		
	
	simFold <- lapply(fold_seq, function(x) NBsim(foldDiff = x, dataset = dataset, group = group, nTags = nTags, nlibs = nlibs, fix.dispersion = fix.dispersion, lib.size = lib.size, drop.low.lambda = drop.low.lambda, drop.extreme.dispersion = drop.extreme.dispersion,  add.outlier = add.outlier, outlierMech = outlierMech, pOutlier = pOutlier, min.factor = min.factor, max.factor = max.factor, pDiff = pDiff, pUp=pUp, name = name))
	
    simFold <- re_split(simFold)
    for(i in seq_along(outlierMech))
       simFold <- rere_split(simFold, outlierMech[i])
	simFold$fold_seq <- names(fold_seq)
    simFold <- new("FoldList", simFold)
    simFold
}


getVersion <-
function(x)
{
  ## it is low-level function of pval ##
  x1 <- gsub("(\\_)(\\w+)", "", x)	
  v <- unlist(lapply(x1, function(z) {options(warn = -1)
                                      desp<- packageDescription(z)
                                      if(length(desp) == 1)
                                        return("unknown")
                                      else desp$Version
                                      }))
  paste0(x,"_", v )	
}	

rmVersion <-
function(x)
{
 ## it is low-level function of pval ##
 x1 <- strsplit(x, "\\_")
 x1 <- lapply(x1, function(x) x[-length(x)])	
 sapply(x1, paste0, collapse = "_")	
	
}

odd <- function(x) 
{   
	## it is low-level function of pval ##
	y <- seq(x)
	idx <- y %% 2 != 0
	x[idx]
}

mainShow <-
function(count.type, count.name, group, pOutlier)
{ 
	## it is low-level function of pval ##
	pOutlier <- paste(100*pOutlier, "% ", "outliers", sep = "")
	group <- as.factor(group)
	group <- paste0(sum(group == levels(group)[1]), "vs", sum(group == levels(group)[2]))
	if(count.type == "counts")
	paste0("No outliers", "/", count.name, "/", group)
	else
	paste0(pOutlier, "/", count.type, "/",  count.name, "/", group)
	
}


resetPar <- function() {
    ## this re-set args of par for plot ##
	dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}




DESeq2.pfun <-
function(counts, group, design = NULL, mc.cores = 4)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  dse <- DESeq(dse)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj)
      out[is.na(out)] <- 1
      out
  }

DESeq2_rmNA.pfun <-
function(counts, group, design = NULL, mc.cores = 4)
{   
	## implement DESeq2 turning off cooksCutoff ##
	library(DESeq2)
	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	colData(dse)$group <- as.factor(colData(dse)$group)
	dse <- DESeq(dse)
	res <- results(dse, cooksCutoff = FALSE)
	cbind(pval = res$pvalue, padj = res$padj)
}

DESeq_cr.pfun <-
function(counts, group, design = NULL, mc.cores = 4)
{   
	## implement DESeq using cox-reid adjust likelihood to estimate dispersion ##
	library(DESeq)
	de <- newCountDataSet(counts, group)
	de <- estimateSizeFactors(de)
	de <- estimateDispersions(de, method = "pooled-CR")
	res <- nbinomTest(de, levels(group)[1], levels(group)[2])
    cbind(pval = res$pval, padj = res$padj)
}
DESeq_pool.pfun <-
function(counts, group, design = NULL, mc.cores = 4)
{   
	## implement DESeq using pooled method to estimate dispersion ##
	library(DESeq)
	de <- newCountDataSet(counts, group)
	de <- estimateSizeFactors(de)
	de <- estimateDispersions(de, method = "pooled")
	res <- nbinomTest(de, levels(group)[1], levels(group)[2])
    cbind(pval = res$pval, padj = res$padj)
}

DESeq_glm.pfun <-
function(counts, group, design = NULL, mc.cores = 4)
{   
    ## implement DESeq via the NB GLMs ##
	library(DESeq)
	de <- newCountDataSet(counts, group)
	de <- estimateSizeFactors(de)
	de <- estimateDispersions(de)
	fit1 = fitNbinomGLMs(de, count ~ group ) 
	fit0 = fitNbinomGLMs(de, count ~ 1 ) 
    pval <- nbinomGLMTest( fit1, fit0 )
	padj <- p.adjust(pval, method="BH" )
	cbind(pval = pval, padj = padj)
}

edgeR.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- calcNormFactors(d)
	d <- estimateGLMCommonDisp(d, design = design)
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	f <- glmFit(d, design = design)
	lr <- glmLRT(f, coef=2)
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}

edgeR_robust.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{   
    ## edgeR-robsut pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- calcNormFactors(d)
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6)
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}

edgeR_rdev.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{   
    ## edgeR-robust by deviance residual ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- calcNormFactors(d)
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6, residual.type = "deviance")
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}

edgeR_rans.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{   
	## edgeR-robust by anscombe residual ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- calcNormFactors(d)
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6, residual.type = "anscombe")
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}




limma_voom.pfun <-
function(counts, group, design = NULL, mc.cores = 4) 
{   
	## limma voom pipeline ##
	library(limma)
	nf <- calcNormFactors(counts)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf)
	fit <- lmFit(y, design)
	fit <- eBayes(fit)
	pval <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")$P.Value
	padj <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")$adj.P.Val
	cbind(pval = pval, padj = padj)
}


baySeq.pfun <-
function(counts, group, design = NULL, p.out = "pvalue", mc.cores = 4)
{   
    ## baySeq pipeline ## 
	library(baySeq)
	library(snow)
	cl <- snow:::makeCluster(mc.cores, "SOCK")
	group <- as.numeric(as.character(group))
	cd <- new("countData", data = counts, replicates = group, groups = list(NDE = rep(1, length(group)), DE = group))
	cd@libsizes <- getLibsizes(cd)
	cd <- getPriors.NB(cd, equalDispersions = TRUE, estimation = "QL", cl = cl)
	cd <- getLikelihoods.NB(cd,pET = "BIC", cl = cl)
	try(snow:::stopCluster(cl), silent = TRUE)
	cd.table <- topCounts(cd, group = "DE", number = nrow(counts))
	id <- match( row.names(counts), row.names(cd.table))
	padj <- cd.table[id,]$FDR
	cbind(pval = pval, padj = padj)	#pval, padj are identical
}


samr_SAMseq.pfun <- 
function(counts, group, design = NULL, mc.cores = 4)
{
  ## SAMseq pipeline ##
  library(samr)
  f <- SAMseq(counts, group, resp.type = "Two class unpaired", fdr.output = 1)
  f.table = rbind(f$siggenes.table$genes.up, f$siggenes.table$genes.lo)
	fdr = rep(NA, nrow(counts)) #contains NA value
  fdr[as.numeric(f.table[, "Gene Name"])] <- as.numeric(f.table[, "q-value(%)"])/100  
  padj <- fdr
  cbind(pval = pval, padj = padj) #pval, padj are identical
}


#ShrinkBayes.pfun <- 
#function(counts, group, design, p.out = "pvalue", mc.cores)
#{   library(ShrinkBayes)
#	try(sfStop(), silent = TRUE)
#	try(sfRemoveAll(), silent = TRUE)
#	try(rm(g, envir = .GlobalEnv), silent = TRUE)
#	g <- as.factor(group)
#    assign("g", g, .GlobalEnv)
#	form <- y ~ 1 + g
#	form0 <- y ~ 1
#	shrinksimul <- ShrinkSeq(form = form, dat = counts, shrinkfixed = "g", fams = "zinb", ncpus = mc.cores)
#	fitall <- FitAllShrink(form, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
#	fitall0 <- FitAllShrink(form0, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
#	ptm <- proc.time()
#	npprior <-MixtureUpdatePrior(fitall = fitall, fitall0 = fitall0, shrinkpara="g", ncpus = mc.cores)
#	print(proc.time() - ptm)
#	nppostshr <- MixtureUpdatePosterior(fitall, npprior, fitall0, ncpus = mc.cores)
#	lfdr <- SummaryWrap(nppostshr)
#	fdr <- BFDR(lfdr)
#    rm(g, envir = .GlobalEnv)
#    pval <- as.vector(fdr)
#}


ShrinkBayes.pscript <-  expression(
   {
    ## ShrinkBayes pipeline implemented on special envir ##  
	ptm <- proc.time()
    library(ShrinkBayes)
	g <- as.factor(group)
	form <- y ~ 1 + g
	form0 <- y ~ 1
	shrinksimul <- ShrinkSeq(form = form, dat = counts, shrinkfixed = "g", fams = "zinb", ncpus = mc.cores)
	fitall <- FitAllShrink(form, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
	fitall0 <- FitAllShrink(form0, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)  
	npprior <-MixtureUpdatePrior(fitall = fitall, fitall0 = fitall0, shrinkpara="g", ncpus = mc.cores)
	nppostshr <- MixtureUpdatePosterior(fitall, npprior, fitall0, ncpus = mc.cores)
	lfdr <- SummaryWrap(nppostshr)
	print(proc.time() - ptm)
	try(sfStop(), silent = TRUE)
	try(sfRemoveAll(), silent = TRUE)
	fdr <- BFDR(lfdr)
	pGlobal <- as.vector(fdr)
	rm(list = c("ptm", "g", "form", "form0", "shrinksimul", "fitall", "fitall0", "npprior" ,"nppostshr", "lfdr", "fdr"))
	padj <- pGlobal
	cbind(pval = pval, padj = padj) #pval, padj are identical
	})

EBSeq.pfun <- 
function(counts, group, design = NULL, mc.cores = 4)
{
  ## EBSeq pipeline ##
  library(EBSeq)
  group <- as.factor(group)
  sizes = MedianNorm(counts)
  f <- EBTest(Data = counts, Conditions = group, sizeFactors = sizes, maxround = 5)
  pp = GetPPMat(f)
  padj <- fdr <- 1 - pp[, "PPDE"]
	if(!length(pval) == nrow(counts)) #check rm 0 counts
	{
	 i <- match(names(padj), rownames(counts))
	 Padj <- rep(NA, nrow(counts))
	 Padj[i] <- padj
	 padj <- Padj	
	}	
  cbind(pval = pval, padj = padj) #pval, padj are identical
}	





pval <-
function(y, ...) ## evaluate DE methods ##
UseMethod("pval")
pval.default <-
function(y, group, design = NULL, method = "edgeR", mc.cores = 4, globalEnvir = FALSE, ...)
{   
	## evaluate DE methods ##
    ## return to a list of pvalue and runing time ##
    ## pvalue contains pvalue and p-adjust value ##
	gc(FALSE)
	time <- proc.time()
	group <- as.factor(group)
	if(globalEnvir) method <- paste0(method, ".pscript")
	else method <- paste0(method, ".pfun")
	p <- get(method)
	if(globalEnvir)
	{
		L <- list(counts = y, group = group, design = design, mc.cores = mc.cores, p = p)
		e <- list2env(L, envir = .GlobalEnv)
		pvalue <- with(e, eval(p))
		try(rm(list = names(L), envir = e),  silent = TRUE)
		try(rm(pGlobal, envir = e),  silent = TRUE)
	}	
	else pvalue <- p(y, group, design, mc.cores, ...)
	pvalue
	new.time <- proc.time()
	output <- list(pvalue = pvalue, time = new.time - time)
}



pval.DGEList <-
function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE, name = deparse(substitute(y)))
{   
	## evaluate DE methods ##
    ## return to a DGEList including pvalue and other necessary indexs for re-analysis and plot ## 
	library(parallel)	   
	count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
	if(count.type == "counts") 
	{counts = y$counts
	 pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA}
	else 
	{counts = y[[count.type]]$countAddOut
	 mask_outlier = y[[count.type]]$mask_outlier
	 pOutlier = y[[count.type]]$pOutlier	
	 indDEupOutlier = y[[count.type]]$indDEupOutlier
	 indDEdownOutlier = y[[count.type]]$indDEdownOutlier
	 indDEbothOutlier = y[[count.type]]$indDEbothOutlier	
	 indDEnoOutlier = y[[count.type]]$indDEnoOutlier
	 indNonDEnoOutlier = y[[count.type]]$indNonDEnoOutlier
	 indNonDEOutlier = y[[count.type]]$indNonDEOutlier}
	names(method) <- method 
	group <- y$group
	design <- y$design
	is.parallel <- method %in% parallel.method
	is.globalEnvir <- method %in% globalEnvir.method
	id.re <- !(is.parallel|is.globalEnvir)
	reduced.method <- method[id.re]
	if(any(id.re)) output <-  parallel:::mclapply(reduced.method, function(x) pval(y = counts, group = group, design = design, method = x), mc.cores = mc.cores, mc.preschedule = FALSE)
	else output <- list()
	if(any(is.parallel))
	{
	  for( i in names(method[is.parallel]))
         output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores)
	}
	if(any(is.globalEnvir))
	{  		
		for( i in names(method[is.globalEnvir]))
		output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE)
	}
	output <- output[method]
	padj <- lapply(output, function(x) x[["pvalue"]][, "padj"])
	pval <- lapply(output, function(x) x[["pvalue"]][, "pval"])
    time <- lapply(output, function(x) x[["time"]])
	output <- new("DGEList", list(pval = pval, padj = padj,  counts = counts, count.type = count.type, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time)) 
	output$main <- mainShow(count.type = count.type, count.name = y$name, group = group, pOutlier = pOutlier)
	output$methodVersion <- getVersion(method)
	if(save.file)
	{   Type <- switch(output$count.type, counts = "b", S = "s", M = "m", R = "r")
	  	mnm <- output$methodVersion
		foldernm <- rmVersion(mnm)
		dir_mnm <- dirname(name)
		filenm <- gsub("(\\.)(\\w+)", "", basename(name))
		foldernm <- paste0(dir_mnm,"/", foldernm)
		lapply(foldernm, dir.create, showWarnings = FALSE)
		mapply(function(u, v, w) {assign(u, v)
			   save(list = u, file = paste0(w, "/pval_", Type, "_", filenm, "_", u, ".Rdata"))}, u = mnm, v = output$pval, w = foldernm)
		mapply(function(u, v, w) {assign(u, v)
			   save(list = u, file = paste0(w, "/padj_", Type, "_", filenm, "_", u, ".Rdata"))}, u = mnm, v = output$padj, w = foldernm)
	}	
	output
}

pval.character <-
function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE)
{   
    ## for shiny app ##
	fnm <- y
	load(y)
	name <- gsub("(\\.)(\\w+)", "", basename(y))
	y <- get(name)
	pval.DGEList(y = y, method = method, count.type = count.type, mc.cores = mc.cores, parallel.method = parallel.method, globalEnvir.method = globalEnvir.method, save.file = save.file, name = fnm)
}
pval.FoldList <-
function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), cut.computing = TRUE)
{   
    ## evaluate DE methods for FoldList ##
	library(parallel)
	count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
	if(count.type == "counts") 
	{counts = y$counts
		pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA}
	else 
	{counts = y[[count.type]]$countAddOut
		pOutlier = y[[count.type]]$pOutlier
		mask_outlier = y[[count.type]]$mask_outlier
		indDEupOutlier = y[[count.type]]$indDEupOutlier
		indDEdownOutlier = y[[count.type]]$indDEdownOutlier
		indDEbothOutlier = y[[count.type]]$indDEbothOutlier	
		indDEnoOutlier = y[[count.type]]$indDEnoOutlier
		indNonDEnoOutlier = y[[count.type]]$indNonDEnoOutlier
		indNonDEOutlier = y[[count.type]]$indNonDEOutlier}
	
	names(method) <- method 
	group <- y$group
	design <- y$design
	is.parallel <- method %in% parallel.method
	is.globalEnvir <- method %in% globalEnvir.method
	id.re <- !(is.parallel|is.globalEnvir)
	reduced.method <- method[id.re] 
	if(any(id.re)) output <- lapply(reduced.method, function(x) parallel:::mclapply(counts, function(w) pval(y = w, method = x, group = group, design = design), mc.cores = mc.cores))
	else output <- list()
	fold_seq <- fold_seq.keep <- y$fold_seq
	if(cut.computing) fold_seq.keep <- odd(fold_seq)
	if(any(is.parallel))
	{
		for(j in fold_seq)
		{ 
		  is.keep <- j %in% fold_seq.keep
		  for( i in names(method[is.parallel]))
			{
				if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores)
				else 
				{
					output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
					output[[i]][[j]][["time"]] <- NA
				}
			}	
		}
	}
	if(any(is.globalEnvir))
	{
		for(j in fold_seq)
		{ 
			is.keep <- j %in% fold_seq.keep
			for( i in names(method[is.globalEnvir]))
			{
				if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE)
				else 
				{
					output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
					output[[i]][[j]][["time"]] <- NA
				}
			}	
		}
	}
	output <- output[method]
	padj <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "padj"]), silent = TRUE)
	pval <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "pval"]), silent = TRUE)
	time <- try(lapply(output, lapply, function(x) x[["time"]]), silent = TRUE)
	output <- new("FoldList", list(fold_seq = y$fold_seq, pval = pval, padj = padj, counts = counts, count.type = count.type, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier,indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time))
	output$main <- mainShow(count.type = count.type, count.name = y$name, group = group, pOutlier = pOutlier)
	output$methodVersion <- getVersion(method)
	output
}
getPvalVersion <- function(methodVersion, pout = "pval", count.type = "counts", datanm)
{ 
  ## for shiny app ##
  Type <- switch(count.type, counts = "b", S = "s", M = "m", R = "r")
  filenm <- paste0(pout, "_", Type, "_",  basename(datanm), "_", methodVersion, ".Rdata")
  load(paste0(dirname(datanm),"/", rmVersion(methodVersion),"/", filenm))
  get(methodVersion)	   
}	
getPval <- function(dataset,methodVersion, count.type = c("counts", "S", "R", "M"))
{
    ## for shiny app ##
	load(dataset)
	datanm <- gsub("(\\.)(\\w+)", "", dataset)
	dat <- get(basename(datanm))
	count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
	if(count.type == "counts") Dat <- new("DGEList", dat)
	else 
	{
		Dat <- new("DGEList", dat[[count.type]])
		Dat[["counts"]] <- Dat[["countAddOut"]]
	}
	Dat$method <- Dat$methodVersion <- methodVersion
	Dat$group = dat$group
	Dat$indDE = dat$indDE
	Dat$name = dat$name
	Dat$main <- mainShow(count.type = count.type, count.name = Dat$name, group = Dat$group, pOutlier = Dat$pOutlier)
	index <- c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier")
	names(index) <- index
	indDiff <- lapply(index, function(x) Dat[[x]])
	indDiff <- indDiff[!sapply(indDiff, is.null)]
	indDiff <- indDiff[!is.na(indDiff)]
	Dat$index <- names(indDiff)
	names(methodVersion) <- methodVersion
	Dat[["padj"]] <- lapply(methodVersion, getPvalVersion, pout = "padj", count.type = count.type, datanm = datanm)
	Dat[["pval"]] <- lapply(methodVersion, getPvalVersion, pout = "pval", count.type = count.type, datanm = datanm)
	Dat 
	
}	
roPlot <-
function(y, ...) 
UseMethod("roPlot")
## plot ROC curve ##

roPlot.default <-
function(y, indDiff, plot.max.fpr = 0.4, add = FALSE, cex.axis = 2, threshold = 0.05, col = 1, cex.threshold = 3, ...)
{   
    ## plot ROC curve ##
	old.par <- par(c("mar", "mgp", "cex.axis"))
	par(mar=c(4,5,3,2))
    par(mgp = c(2.6, 1, 0))
	par(cex.axis = cex.axis)
	on.exit(par(old.par))
	library(ROCR)
	if(any(is.na(y)))
	{
      y[is.na(y)] <- 1
	}	
	y = 1 - y
	label <- as.factor(rep("nonDE", length(y)))
	levels(label) <- c("nonDE", "DE")
	label[indDiff] <- "DE"
	pred <- prediction(y, label, label.ordering = c("nonDE", "DE"))
	perf <- performance(pred, "tpr", "fpr")
	if(is.null(plot.max.fpr))
	plot.max.fpr <- 1
	plot(perf, xlim = c(0, plot.max.fpr), col = col, add = add, ...)
	if(!is.null(threshold))
	{
		fpr <- approx(y = perf@x.values[[1]], x = perf@alpha.values[[1]], xout = 1- threshold)$y
		tpr <- approx(y = perf@y.values[[1]], x = perf@x.values[[1]], xout = fpr)$y
		points(x = fpr, y = tpr, pch = 4, col = col, cex = cex.threshold, ...)
	}	
}
roPlot.DGEList <-
function(y, plot.max.fpr = 0.4, pout = "padj", threshold = 0.05, selected.method = NULL, show.legend = TRUE, short.main = FALSE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.6, cex.threshold = 8)
{   
	## plot ROC curve for DGEList ##
	library(ROCR)
	main <- y$main
	if(short.main) 
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")

	}
	if(is.null(selected.method)) 
	{
		method <- y$method
		methodVersion <- y$methodVersion
	}
	else 
	{
		method <- match.arg(selected.method, y$method, several.ok = TRUE)
		methodVersion <- y$methodVersion[match(method, y$method)]
	}
	pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
	pvalue <- y[[pout]][method]
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]		
	roPlot(y = pvalue[[1L]], indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold)
	if(is.list(pvalue)) mapply(function(u, v, w) roPlot(y = u, indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold, add = TRUE), u = pvalue[-1L], v = col[-1L])
	if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)	
	if(!is.null(threshold)&show.legend) legend("topleft", paste0("P_adj_value=", threshold), col = "black", pch = 4, lty = NA, lwd = lwd, box.lwd = box.lwd, cex = cex, pt.cex = 1.5*cex)	
}



fdPlot <-
function(y, ...) 
UseMethod("fdPlot")
## plot False discovery number curve ##

fdPlot.default <-
function(y, indDiff, add=FALSE, xlab="Number of genes selected", 
ylab="Number of false discoveries", lwd=4, type="l", ... ) 
{
	
	## plot False discovery number curve ##	
	old.par <- par(c("mar", "mgp"))
	par(mar=c(4,5,3,2))
    par(mgp = c(2.6, 1, 0))
	on.exit(par(old.par))	
	x <- 1:length(indDiff)
	o <- order(y)
	w <- !o[x] %in% indDiff
	y1 <- cumsum(w)
	matplot(x, y1, xlab=xlab, ylab=ylab, lwd=lwd, type=type, add=add, ...)
}
fdPlot.DGEList <-
function(y, pout = "padj", selected.method = NULL, short.main = FALSE, show.legend = TRUE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.3, xlim = NULL)
{   
	## plot False discovery number curve for DGEList ## 
	main <- y$main
	if(short.main) 
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
	}
	if(is.null(selected.method)) 
	{
		method <- y$method
		methodVersion <- y$methodVersion
	}
	else 
	{
		method <- match.arg(selected.method, y$method, several.ok = TRUE)
		methodVersion <- y$methodVersion[match(method, y$method)]
	}
	pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
	pvalue <- y[[pout]][method]

	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]	
	fdPlot(y = pvalue[[1L]], indDiff = y$indDE,  log="y", col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, xlim = xlim)
	if(is.list(pvalue)) mapply(function(u, v, w) fdPlot(y = u, indDiff = y$indDE, log="y", col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, add = TRUE), u = pvalue[-1L], v = col[-1L])
	if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)
}

getPower <-
function(y, ...) 
UseMethod("getPower")
## plot Power curve ##
getPower.default <-
function(y, indDiff, threshold)
{ 
	## plot Power curve ##
	if(all(is.na(y)))
    NA
	else if(all(is.na(indDiff)))
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "fpr")
	else
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "power")
}
getPower.DGEList <-
function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
{
	## plot Power curve for DGEList ##
	index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
	main <- y$main
	if(short.main)
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
    
    if(is.null(selected.method))
    {
        method <- y$method
        methodVersion <- y$methodVersion
    }
    else
    {
        method <- match.arg(selected.method, y$method, several.ok = TRUE)
        methodVersion <- y$methodVersion[match(method, y$method)]
    }
	pout <- match.arg(pout, c("pval", "padj"))
	pvalue <- y[[pout]][method]
	indDiff <- y[[index]]
	
	if(byAveLogCPM)
    {
        threshold <- rep(threshold, cutCPM)
        diff <- rep(FALSE, nrow(y$counts))
        diff[indDiff] <- TRUE
        d <- DGEList(counts = y$counts, group = y$group)
        d <- calcNormFactors(d)
        AveLogCPM <- aveLogCPM(d)
        o <- order(AveLogCPM)
        l <- length(o)/cutCPM
        oo <- split(o, ceiling(seq_along(o)/l))
        AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
        nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
        nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))
        
        power <- list()
        power[["all"]] <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold[1]))
        for( i in 1:length(oo))
        {
            diff_idx <- which(diff[oo[[i]]])
            pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
            power[[nm[i]]] <- unlist(lapply(pvalue_idx, getPower, indDiff = diff_idx, threshold = threshold[i]))
        }
        power <- do.call("rbind", power)
	}
    else
    {
        power <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold))
	}
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$power <- power
	out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
	out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
    if(byAveLogCPM) matPlot(out, output = "power", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else powBarPlot(out, col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
	out
}



getFDR <-
function(y, ...)
UseMethod("getFDR")
## plot Power curve ##
getFDR.default <-
function(y, indDiff, threshold)
{
	## plot Power curve ##
	if(all(is.na(y)))
    NA
	else if(all(is.na(indDiff)))
    NA
	else
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "fdr")
}
getFDR.DGEList <-
function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
{
	## plot Power curve for DGEList ##
	index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
	main <- y$main
	if(short.main)
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
    
    if(is.null(selected.method))
    {
        method <- y$method
        methodVersion <- y$methodVersion
    }
    else
    {
        method <- match.arg(selected.method, y$method, several.ok = TRUE)
        methodVersion <- y$methodVersion[match(method, y$method)]
    }
	pout <- match.arg(pout, c("pval", "padj"))
	pvalue <- y[[pout]][method]
	indDiff <- y[[index]]
	
	if(byAveLogCPM)
	{
        threshold <- rep(threshold, cutCPM)
        diff <- rep(FALSE, nrow(y$counts))
        diff[indDiff] <- TRUE
        d <- DGEList(counts = y$counts, group = y$group)
        d <- calcNormFactors(d)
        AveLogCPM <- aveLogCPM(d)
        o <- order(AveLogCPM)
        l <- length(o)/cutCPM
        oo <- split(o, ceiling(seq_along(o)/l))
        AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
        nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
        nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))
        
        fdr <- list()
        fdr[["all"]] <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold[1]))
        for( i in 1:length(oo))
        {
            diff_idx <- which(diff[oo[[i]]])
            pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
            fdr[[nm[i]]] <- unlist(lapply(pvalue_idx, getFDR, indDiff = diff_idx, threshold = threshold[i]))
        }
        fdr <- do.call("rbind", fdr)
	}
    else
    {
        fdr <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold))
	}
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$fdr <- fdr
	out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
	out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
    if(byAveLogCPM) matPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else barPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
	out
}


getTable <-
function(y, ...)
UseMethod("getTable")
## plot Power curve ##
getTable.default <-
function(y, indDiff, threshold)
{
	## plot Power curve ##
	if(all(is.na(y)))
    NA
	else if(all(is.na(indDiff)))
    NA
	else
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "table")
}
getTable.DGEList <-
function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE)
{
	## plot Power curve for DGEList ##
	index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
	main <- y$main
	if(short.main)
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
    
    if(is.null(selected.method))
    {
        method <- y$method
        methodVersion <- y$methodVersion
    }
    else
    {
        method <- match.arg(selected.method, y$method, several.ok = TRUE)
        methodVersion <- y$methodVersion[match(method, y$method)]
    }
	pout <- match.arg(pout, c("pval", "padj"))
	pvalue <- y[[pout]][method]
	indDiff <- y[[index]]
	
	if(byAveLogCPM)
	{
        diff <- rep(FALSE, nrow(y$counts))
        diff[indDiff] <- TRUE
        d <- DGEList(counts = y$counts, group = y$group)
        d <- calcNormFactors(d)
        AveLogCPM <- aveLogCPM(d)
		AveLogCPM <- cut(AveLogCPM, breaks=quantile(AveLogCPM,p=(0:cutCPM)/cutCPM), include.lowest=TRUE, dig.lab = 1)
        tab <- list()
        tab[["all"]] <- lapply(pvalue, getTable, indDiff = indDiff, threshold = threshold)
        for( i in levels(AveLogCPM))
        {
            idx <- AveLogCPM == i
            diff_idx <- which(diff[idx])
            pvalue_idx <- lapply(pvalue, function(x) x[idx])
            tab[[i]] <- lapply(pvalue_idx, getTable, indDiff = diff_idx, threshold = threshold)
        }
        tab1 <- list()
        for(i in method)
        tab1[[i]] <- lapply(tab, .subset2, i)
        tab <- tab1
	}
    else
    {
        tab <- lapply(pvalue, getTable, indDiff = indDiff, threshold = threshold)
	}
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$tab <- tab
	out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
	out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    out
}




contable <-
function(score, indDiff, threshold, output = "power")
{ 
    ## this is the low-level function of getPower to build a contengency table calculating fpr, tpr(power) ##
	output <- match.arg(output, c("table", "power", "fdr", "fpr"))
	if(any(is.na(score)))
	{
		score[is.na(score)] <- 1
	}
    if(all(score >= threshold))
        return(output <- switch(output, table = NA, power = 0, fpr = 0, fdr = NA))
    pred <- cut(score, breaks = c(0, threshold, 1))
	levels(pred) <- c("s","n")
	label <- labelS <- as.factor(rep("nonDE", length(score)))
	levels(label) <- c("nonDE", "DE")
	label[indDiff] <- "DE"
	tab <- table(pred, label)
	FP <- tab["s", "nonDE"]
	TN <- tab["n", "nonDE"]
	TP <- tab["s", "DE"]
	FN <- tab["n", "DE"]
	power <- TP/(TP + FN)
	fpr <- FP/(FP + TN)
    fdr <- FP/(FP + TP)
	output <- switch(output, table = tab, power = power, fpr = fpr, fdr = fdr)
	
}
	

barPlot <-
function(object, output = c("fdr", "power"), col = 1, cex.main = 2.1, cex.lab = 1, cex.axis = 1, cex.sub = 1, ylim = NULL,...)
{      output <- match.arg(output, c("power", "fdr"))
       ## this is the low-level function of getPower to barplot ##
	   old.par <- par(c("mar", "mgp"))
	   par(mar=c(12,5,3,2))
	   par(mgp = c(2.6, 1, 0))
	   on.exit(par(old.par))
	   score <- object[[output]]
	   if(is.null(ylim)) ylim <- c(0, min(round(1.2*max(score), 1), 1))
	   x <- barplot(object[[output]], col = col, border="orange", main = object$main, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, xaxt = "n", ylab = output, ylim = ylim, ...)
	   text(x, par("usr")[3]-0.01*(par("usr")[4]-par("usr")[3]), labels = object$methodVersion, pos = 2, xpd = TRUE, cex = cex.lab-0.5, srt=45)
	   text(0.5*(par("usr")[2] + par("usr")[1]), par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]), object$index, cex = 0.8*cex.main)	   	   
	   
}

matPlot <-
function(object, output = c("fdr", "power"), show.legend = TRUE, col = 1, cex.main = 2.1, cex.lab = 1, cex.axis = 1, cex.sub = 1, ylim = NULL, lwd = 8, cex = 2.5, pch = 20, box.lwd = 1.5, cex.legend = 1, ...)
{      	
	output <- match.arg(output, c("power", "fdr"))
    ## this is the low-level function of getPower to matplot ##
	old.par <- par(c("mar", "mgp"))
	par(mar=c(5, 5, 3, 2))
	par(mgp = c(3.2, 0.7, 0))
	on.exit(par(old.par))
	
    score <- object[[output]]
	if(is.null(ylim)) ylim <- c(0, min(round(1.2*max(score, na.rm = TRUE), 1), 1))
	if(show.legend) xlim <- c(1, nrow(score)+2)
	else xlim <- c(1, nrow(score)+0.5)
	
	matplot(score, col = col, main = object$main, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, xaxt = "n", ylab = output, xlab = "AveLogCPM", ylim = ylim, xlim = xlim, lty = 3, lwd = lwd, type = "p", cex = cex, pch = pch, ...)
	score["all",] <- NA
	matplot(score, add = TRUE, col = col, lty = 3, lwd = lwd, type = "o", cex = cex, pch = pch, ...)
	axis(1, at = seq(nrow(score)), labels = NA, cex.axis = cex.axis)
	text(seq(nrow(score)), par("usr")[3]-0.025*(par("usr")[4]-par("usr")[3]), labels = rownames(score), pos = 1, xpd = TRUE, cex = 0.7*cex.lab, srt=10)
	text(0.5*(par("usr")[2] + par("usr")[1]), par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]), object$index, cex = 0.8*cex.main)
	if(show.legend) legend("bottomright", object$methodVersion, col = col, lty = NA, pch = 19, lwd = lwd, box.lwd = box.lwd, text.font=2, cex = cex.legend) 
	
}



summPlot <- function(pval_b, pval_o, ylim = c(0, 1), byAveLogCPM = FALSE, ...)
{   
    ## this function provides a summary plots including roPLot, fdPlot, getPower separated by counts and counts with outliers ##
	old.par <- par(c("mar", "mfrow", "oma", "xpd"))
	on.exit(par(old.par))	
	if(all(is.na(pval_o$indDEbothOutlier)))
	par(mfrow = c(3, 3))
	else par(mfrow = c(4, 3))
	par(mar=c(5,5,4,2))
	par(oma=c(0,0,3, 0))
	
	fdPlot(pval_b, short.main = TRUE, show.legend = FALSE, ...)
	mtext("(a)", side=3, adj=-.1, padj=-.7, cex=2.5)
	fdPlot(pval_o, short.main = TRUE, show.legend = FALSE, ...)
	mtext("(d)", side=3, adj=-.1, padj=-.7, cex=2.5)
	getPower(pval_o, index = "indDEnoOutlier", plot = TRUE, byAveLogCPM =  byAveLogCPM, short.main = TRUE, ylim = ylim, show.legend = FALSE, ...)
	mtext("(g)", side=3, adj=-.1, padj=-.7, cex=2.5)	
	
	roPlot(pval_b, short.main = FALSE, show.legend = FALSE, ...)
	mtext("(b)", side=3, adj=-.1, padj=-.7, cex=2.5)
	roPlot(pval_o, short.main = TRUE, show.legend = TRUE, ...)
	mtext("(e)", side=3, adj=-.1, padj=-.7, cex=2.5)
	getPower(pval_o, index = "indDEupOutlier", plot = TRUE, byAveLogCPM =  byAveLogCPM, short.main = TRUE, ylim = ylim, show.legend = FALSE, ...)
	mtext("(h)", side=3, adj=-.1, padj=-.7, cex=2.5)
	
	getPower(pval_b, plot = TRUE, short.main = TRUE, byAveLogCPM =  byAveLogCPM, ylim = ylim, show.legend = FALSE, ...)
	mtext("(c)", side=3, adj=-.1, padj=-.7, cex=2.5)
	getPower(pval_o, plot = TRUE, short.main = TRUE, byAveLogCPM =  byAveLogCPM, ylim = ylim, show.legend = FALSE, ...)
	mtext("(f)", side=3, adj=-.1, padj=-.7, cex=2.5)
	getPower(pval_o, index = "indDEdownOutlier", plot = TRUE, byAveLogCPM =  byAveLogCPM, short.main = TRUE, ylim = ylim, show.legend = FALSE, ...)
	mtext("(i)", side=3, adj=-.1, padj=-.7, cex=2.5)
	
	if(!all(is.na(pval_o$indDEbothOutlier)))
	{
		plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
		plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
		getPower(pval_o, index = "indDEbothOutlier", plot = TRUE,  byAveLogCPM =  byAveLogCPM, short.main = TRUE, ylim = ylim, ...)
		mtext("(j)", side=3, adj=-.1, padj=-.7, cex=2.5)	
	}	
	
	
     #mtext(outer = TRUE, cex = 1.8, "Simulation Pickrell with high outliers, 5vs5, foldDifference=3, pDifferential=0.1, pUp = 0.5", side = 3, line = 1)
par(xpd=NA)
	rect( grconvertX(0.001, from='ndc'), grconvertY(0.005, from='ndc'), grconvertX(0.327, from='ndc'), grconvertY(0.995, from='ndc'), lwd = 2.5)
	rect( grconvertX(0.333, from='ndc'), grconvertY(0.005, from='ndc'), grconvertX(0.660, from='ndc'), grconvertY(0.995, from='ndc'), lwd = 2.5)
	rect( grconvertX(0.667, from='ndc'), grconvertY(0.005, from='ndc'), grconvertX(0.995, from='ndc'), grconvertY(0.995, from='ndc'), lwd = 2.5)
	
}




fdFdTp <- function(y, indDiff, cutoff.p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1), mc.cores = 6)
{   
	## this is the low-level function of summFold ##
	re_split2 <- function(object, rm.identical = TRUE)
        {       
              out <- list()
	    for (what in names(object[[1L]][[1L]])) {
		out[[what]] <- lapply(object, lapply,.subset2, what)
		if(rm.identical)
		{if(all(out[[what]][[1L]][[1L]] == out[[what]])) out[[what]] <- out[[what]][[1L]]}
	     }
	out
         }
		
	
	fdTp <- function(y, indDiff, cutoff.p = 0.05)
	{ 
		if(all(is.na(y)))
		list(tpr = NA, fdr = NA)
		else
		fdTp.fun(y = y, indDiff = indDiff, cutoff.p = cutoff.p) 
	}
	
		
	fdTp.fun <- function(y, indDiff, cutoff.p = 0.05)
	{   
		library(ROCR)
		if(any(is.na(y)))
		{
        #y[is.na(y)] <- runif(sum(is.na(y)))
			y[is.na(y)] <- 1
		}	
		pred <- cut(y, breaks = c(0, cutoff.p, 1))
		levels(pred) <- c("s","n")
		label <- as.factor(rep("nonDE", length(y)))
		levels(label) <- c("nonDE", "DE")
		label[indDiff] <- "DE"
		tab <- table(pred, label)
		FP <- tab["s", "nonDE"]
		TN <- tab["n", "nonDE"]
		TP <- tab["s", "DE"]
		FN <- tab["n", "DE"]
		tpr <- TP/(TP + FN)
		if((FP + TP) == 0) fdr <- 0
		else fdr <- FP/(FP + TP)
		list(tpr = tpr, fdr = fdr)
	}
	
	names(cutoff.p) <- cutoff.p
	library(parallel)
	out <- parallel:::mclapply(y, function(z) mapply(function(u, v) lapply(cutoff.p, function(w) fdTp(y = u, indDiff = v, cutoff.p = w)), u = z, v = indDiff, SIMPLIFY = FALSE), mc.cores = mc.cores)
	out <- lapply(out, re_split2, rm.identical = FALSE)
	out <- lapply(out, lapply, lapply, unlist)
	out <- lapply(out, lapply, do.call, what = "cbind")
}

MatPlot <- function(object, fdr.max = 1,  col, xlab = "", ylab = "",  main = "", add = FALSE,  lty = 1, lwd = 2, yaxt="n", cex.main = 3, cex.axis=1.5, cex.lab = 1.5, xlim = NULL, ylim = c(0, 1), type = "o", pch = 1, las = 1, cex = 1)
{
	
	
	## this is the low-level function of summFold ##
	mPlot <- function(x, y, add = FALSE, ylab = ylab, xlab = xlab, lty = 1, lwd = 2,col = 1, main = "",  cex.main = 3, cex.axis=1.5, cex.lab = 1.5, xlim = c(0, 1), ylim = c(0, 1), type = "o", pch = 1, cex = 1, las = 1)
	{
		if(!add) plot(x[,1], y[,1], lty = lty, col = col, ylab = ylab, xlab = xlab, type = "n", lwd = lwd, ylim = ylim, xlim = xlim, yaxt = yaxt, main = main, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, pch = pch, cex = cex, las = las)
		for(i in seq(ncol(y)))
		lines(x[, i], y[, i], lty = lty, col = col, type = type, lwd = lwd, pch = pch, cex = cex)
		
	}
	
	if(is.null(xlim)) 
	{  
		xmax <- max(pmin(unlist(object)[grep("fdr", names(unlist(object)))], fdr.max), na.rm = TRUE)
	    xlim = c(0, xmax)
	}	
	mPlot(pmin(object[[1]]$fdr, x.max = fdr.max), object[[1]]$tpr, xlim = xlim, add = add, lty = lty, col = col[1], ylab = ylab, xlab = xlab, type = type,lwd = lwd, ylim = ylim, main = main, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, pch = pch, cex = cex, las = las)
	if(length(object) > 1)
	mapply(function(u, v) mPlot(pmin(u$fdr, x.max = fdr.max), u$tpr, lty = 1, col = v, add = TRUE, las = las,type = type, lwd = lwd, pch = pch, cex = cex), u = object[-1], v = col[-1])
	
	
}


summFold <- function(y, pout = "padj", selected.method = NULL, cutoff.p = c(0.02, 0.05, 0.1), fdr.max = 0.2, short.main = FALSE, Plot = TRUE, col = NULL, ylim = NULL, xlim = NULL, mc.cores = 6)

{
    ## this provides a summary plot based on power verse true fdr of a series of simulated counts ##
	indDE <- y$indDE
	main <- y$main
	if(short.main) 
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
	if(is.null(selected.method)) 
	{
		method <- y$method
		methodVersion <- y$methodVersion
	}
	else 
	{
		method <- match.arg(selected.method, y$method, several.ok = TRUE)
		methodVersion <- y$methodVersion[match(method, y$method)]
	}
	pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
	pvalue <- y[[pout]][method]
	
	l <- cut(as.numeric(gsub("fold_", "", y$fold_seq)), breaks = c(0.999, 1, 1.499, 2.5, 2.999, 5, 5.999, 8))
	levels(l)[7] <- "easy"
	levels(l)[5] <- "medium"
	levels(l)[3] <- "hard"
	levels(l)[1] <- "no"
	names(cutoff.p) <- cutoff.p
	
	out_hard <- fdFdTp(lapply(pvalue, function(z) z[l=="hard"]), indDE[l=="hard"], cutoff.p = cutoff.p, mc.cores = mc.cores)  
	out_medium <- fdFdTp(lapply(pvalue, function(z) z[l=="medium"]), indDE[l=="medium"], cutoff.p = cutoff.p, mc.cores = mc.cores)
	out_easy <- fdFdTp(lapply(pvalue, function(z) z[l=="easy"]), indDE[l=="easy"], cutoff.p = cutoff.p, mc.cores = mc.cores) 
	out <- list(hard = out_hard, medium = out_medium, easy = out_easy)
	if(Plot)
	{
		require(grid)
		require(gridBase)
		if(is.null(ylim)) ylim = c(0, 1)	
		pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
		if(is.null(col)) col <- pre.col[seq(method)]	
		if(is.null(fdr.max)) fdr.max <- 3*max(cutoff.p)
		pch <- seq(cutoff.p) 
		opar=par("plt", "mar", "usr")
		plot(0, xlim=c(0, 1), ylim=c(0, 1), type="n", xaxt="n", cex.main = 3, cex.axis=1.5, cex.lab = 1.5, main = main, xlab = "FDR", ylab = "TPR")		
		on.exit(par(opar))
		vps <- baseViewports()
		pushViewport(vps$inner, vps$figure, vps$plot)
 
		labels <- c("hard", "medium", "easy")
		x <- c(0.16666, 0.5, 0.83333)
		y <- rep(0.5, 3)
		for (i in 1:3) {
			pushViewport(viewport(x=unit(x[i], "npc"), y=unit(y[i], "npc"), width= unit(0.33333, "npc"), height = unit(1, "npc")))
			grid.rect()
			par(plt=gridPLT(), new=TRUE)
			MatPlot(out[[i]], fdr.max = fdr.max , xlim = xlim, ylim = ylim, col = col, pch = pch, cex.axis=1, cex.lab = 1.5, cex = 1, lwd = 3, las = 2) 
			
			text(x = 0.5*(par("usr")[2] - par("usr")[1]), y = par("usr")[4] - 0.02,  labels[i], pos = 1, xpd = TRUE, cex = 2)
			if(i == 1) legend(x =  par("usr")[1],y = par("usr")[4]- 0.15, pch = pch, paste0("eFDR=", cutoff.p), text.font=2, pt.lwd = 2)
			upViewport()
			
		}
		legend("bottomright",  methodVersion, col = col, pch = 19, text.font=2)
		
		popViewport(3)
 
		
	}
	out
	
}




checkMethods <- function( )
{   
	## check information of DE methods ##
	methods.pfun <- grep(".pfun|.pscript",ls(envir = .GlobalEnv),value=TRUE)
	methods <- gsub(".pfun|.pscript", "", methods.pfun)
	ip <- installed.packages()
	bioc <- available.packages(contrib.url(get("bioc")))
	sMethods <- get("sMethods")
	sVersions <-get("sVersions")
	rMethods <- get("rMethods")
	rVersions <-get("rVersions")
	biocVersions <- get("biocVersions")
	Rversions <- get("Rversions")
	sSources <-  get("sSources")
	
	out <- matrix(NA , length(methods), 4)
	rownames(out) <- methods
	colnames(out) <- c("pkg", "(current)version", "(recommended)v", "source") 
	
	if( any(ids_bioc <-  (!methods %in% sMethods) &(!methods %in% rMethods)))
	{
		mp <- gsub("(\\_)(\\w+)", "", methods[ids_bioc])
	    cp <- match(mp, rownames(bioc))
	    bioc_mp <- bioc[cp] 
	    names(bioc_mp) <- mp
	    bioc_sugg <- sapply(bioc_mp, function(x) {if(is.na(x)) NA else bioc[x, "Version"]})
		vmp <- sapply(mp, function(x) {if(x %in% ip[,"Package"]) ip[x, "Version"]
					  else NA})
	  out[, "pkg"][which(ids_bioc)] <- mp
	  out[, "(current)version"][which(ids_bioc)] <- vmp
	  out[, "(recommended)v"][which(ids_bioc)] <-  bioc_sugg
	  out[, "source"][which(ids_bioc)] <- biocVersions
	  	
	}
	if( any(ids_r <-  methods %in% rMethods))
	{
		rmp <- gsub("(\\_)(\\w+)", "", methods[ids_r])
		rpkg <- gsub("(\\_)(\\w+)", "", rMethods)
		r_sugg <- rVersions[match(rmp, rpkg)]
		rvmp <- sapply(rmp, function(x) {if(x %in% ip[,"Package"]) ip[x, "Version"]
					  else NA})
		out[, "pkg"][which(ids_r)] <- rmp
		out[, "(current)version"][which(ids_r)] <- rvmp
		out[, "(recommended)v"][which(ids_r)] <-  r_sugg
		out[, "source"][which(ids_r)] <- Rversions
	  	
	}
	if( any(ids_s <-  methods %in% sMethods))
	{
		smp <- gsub("(\\_)(\\w+)", "", methods[ids_s])
		spkg <- gsub("(\\_)(\\w+)", "", sMethods)
		s_sugg <- sVersions[match(smp, spkg)]
		svmp <- sapply(smp, function(x) {if(x %in% ip[,"Package"]) ip[x, "Version"]
					   else NA})
		out[, "pkg"][which(ids_s)] <- smp
		out[, "(current)version"][which(ids_s)] <- svmp
		out[, "(recommended)v"][which(ids_s)] <-  s_sugg
		out[, "source"][which(ids_s)] <- sSources
	  	
	}
	print(out)	    
}

# give versioned intro message
system("clear") # only does something on Mac/Linux
message("----------\n")
## this is version1 of the robust_simulation code
message("Thank you for using the RNA-seq count (and outlier) simulator v0.1.0\n")
message("If you have questions or feedback, please direct them to:\n")
message("  Xiaobei Zhou (xiaobei.zhou@uzh.ch) or \n  Mark Robinson (mark.robinson@imls.uzh.ch)\n")
message("----------\n")

###count data source###
if(!"dataset.Rdata" %in% dir())
{
	message("download dataset into your local dir \n")
    message("download pickrell \n")
    library(tweeDEseqCountData)
    library(edgeR)
    data(pickrell)
    pickrell <- exprs(pickrell.eset)
    pickrell <- as.matrix(pickrell)

    message("download cheung \n")
    cheung <- read.delim("http://bowtie-bio.sourceforge.net/recount/countTables/cheung_count_table.txt", row.names=1, stringsAsFactors=FALSE)
    cheung <- as.matrix(cheung)

    cat("download gilad \n")
    gilad <- read.delim("http://bowtie-bio.sourceforge.net/recount/countTables/gilad_count_table.txt", row.names=1, stringsAsFactors=FALSE)
    gilad <- as.matrix(gilad)
    gilad <- gilad[rowSums(gilad) > 1, ]

    message("download witten \n")
    f <- "http://www.biomedcentral.com/content/supplementary/1741-7007-8-58-s2.xls"
	library(gdata)
    witten <- read.xls(f, sheet = 1, stringsAsFactors=FALSE, skip = 2, row.names = 1, na.strings = "", verbose = FALSE)
    witten <- witten[, -c(ncol(witten) - 0:2)]
    witten <- witten[-c(nrow(witten) - 0:5), ]
    witten <- as.matrix(witten)
    mode(witten) <- "numeric"
    witten[is.na(witten)] <- 0
    save(pickrell, cheung, gilad, witten, file = "dataset.Rdata")
	message("----------\n")
}	
load("dataset.Rdata")

# see what wrappers exist
message("Based on this version of simulation, the following wrappers are already defined:\n")
methods <- grep(".pfun|.pscript",ls(),value=TRUE)
print(methods)
message("----------\n")

methods <- gsub(".pfun|.pscript","",methods)
message("The following methods (as defined by a vector of strings) can be\nused in the 'method' argument of the pval() function below:\n")
print(strwrap(methods))
checkMethods()
message("----------\n")
