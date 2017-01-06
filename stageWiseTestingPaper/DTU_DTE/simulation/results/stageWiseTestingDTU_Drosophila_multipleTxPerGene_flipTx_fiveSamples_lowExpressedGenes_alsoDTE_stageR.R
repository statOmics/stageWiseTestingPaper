# Genome Biology paper of Soneson et al. (2016): Simulation for Drosophila. 
### simulation Drosophila from Genome Biology paper Soneson 2016
baseDir <- "/Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb/drosophila/diffexpression/non_null_simulation_dge/non_null_simulation/"
files=list.files(paste0(baseDir,"1_reads/reads/"),recursive=TRUE)
fastaFiles <- files[grep(x=files,pattern=".fq")]
fastaFiles <- fastaFiles[order(as.numeric(unlist(lapply(strsplit(fastaFiles,split="_"),function(x) x[2]))))]
names(fastaFiles) <- paste0("Dm_sample_",rep(1:10,each=2),"_",rep(1:2,10))
kallistoIndex="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/KallistoIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel"
txConversionFile="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/KallistoIndex/TranscriptID_conversion.txt"
kallistoDir=paste0(baseDir,"quantifications/kallisto/")
truthFile=paste0(baseDir,"3_truth/simulation_details.txt")
library(DEXSeq)
library(iCOBRA)
library(scales)
library(dplyr)

### kallisto quantification
sample <- unlist(lapply(strsplit(names(fastaFiles),split="_"),function(x) x[3]))
fileNames <- unlist(lapply(strsplit(names(fastaFiles),split=".",fixed=TRUE),function(x) x[1]))
fileNames <- sapply(fileNames,function(x) substr(x=x,1,nchar(x)-2))[seq(1,20,2)]
#for(i in 1:10){
#    pairedFasta <- fastaFiles[sample==i]
#    cmd <- paste0("kallisto quant -i ",kallistoIndex,
#	      " -o ",kallistoDir,fileNames[i],
#	      " -b 30 ",
#	      paste0(baseDir,"1_reads/reads/"),pairedFasta[1]," ",paste0(baseDir,"1_reads/reads/"),pairedFasta[2])
#    message(cmd)
#    system(cmd)
#}

### get kallisto results
files2 <- list.files(kallistoDir)
sampleDirs <- files2[grep(x=files2,pattern="Dm_sample_[123456789]")]
#sort them numerically
sampleDirs <- sampleDirs[order(as.numeric(unlist(lapply(strsplit(sampleDirs,split="_"),function(x) x[3]))))]
dir=sampleDirs[1]
hlp=read.table(paste0(kallistoDir,"/",dir,"/abundance.tsv"), header=TRUE) #for rownames
data <- as.data.frame(sapply(sampleDirs,function(dir) read.table(paste0(kallistoDir,"/",dir,"/abundance.tsv"), header=TRUE)[,"est_counts"]), row.names=hlp$target_id)
kal2tx=read.table(txConversionFile)
colnames(kal2tx) <- c("kallisto","transcript")
rownames(data) <- kal2tx$transcript[match(kal2tx$kallisto,rownames(data))]

truth=read.table(truthFile,header=TRUE)
truth_gene <- truth[,c("gene_id","gene_ds_status")]
truth_gene <- truth_gene[!duplicated(truth_gene$gene_id),]
rownames(truth_gene) <- truth_gene$gene_id

truth_tx <- truth[,c("gene_id","transcript_id","gene_ds_status","transcript_ds_status","diff_IsoPct")]
truth_tx$transcript_id <- as.character(truth_tx$transcript_id)
rownames(truth_tx) <- truth_tx[,"transcript_id"]

tx2gene <- truth[,c("gene_id","transcript_id")]
tx2gene$gene_id <- as.character(tx2gene$gene_id)
tx2gene$transcript_id <- as.character(tx2gene$transcript_id)


### DEXSeq analysis
txCount <- ceiling(data)
# remove transcripts with all zero counts,  genes with only one transcript
txCount <- txCount[!rowSums(txCount)==0,]
geneForEachTx <- tx2gene$gene_id[match(rownames(txCount),tx2gene$transcript_id)]
genesWithOneTx <- names(which(table(tx2gene$gene_id[match(rownames(txCount),tx2gene$transcript_id)])==1))
txCount <- txCount[!geneForEachTx %in% genesWithOneTx,]

geneTx <- tx2gene$gene_id[match(rownames(txCount),tx2gene$transcript_id)]
sampleData <- data.frame(condition=factor(rep(c("A","B"),each=5)))
dxd <- DEXSeqDataSet(countData = txCount, 
                         sampleData = sampleData, 
                         design = ~ sample + exon + condition:exon,
                         featureID = rownames(txCount),
                         groupID = geneTx)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxr <- DEXSeqResults(dxd)
hist(dxr$pvalue, main="", xlab="p-value")
qvalDxr <- perGeneQValue(dxr)

## stage-wise DEXSeq analysis using stageR
devtools::install_github("statOmics/stageR",auth_token="cb649b65157aa8cd235a992d99cbe9384fd0eeb2")
library(stageR)
pScreen=qvalDxr
pScreen[is.na(pScreen)]=1 #filtered p-values
pConfirmation <- matrix(dxr$pvalue,ncol=1,dimnames=list(dxr$featureID,"transcript"))
tx2gene <- as.data.frame(cbind(dxr$featureID,dxr$groupID))
#DEXSeq performs independent filtering by default. We will filter the genes that have been filtered in the DEXSeq analysis.
rowsNotFiltered=tx2gene[,2]%in%names(qvalDxr)
pConfirmation=matrix(pConfirmation[rowsNotFiltered,],ncol=1,dimnames=list(dxr$featureID[rowsNotFiltered],"transcript"))
tx2gene <- tx2gene[rowsNotFiltered,]
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
dim(getSignificantGenes(stageRObj)) #this was 816 in original code
dim(getSignificantTx(stageRObj)) #this was 1641 in original code



### characterize the genes found
trueGenes=as.character(truth_gene$gene_id[truth_gene$gene_ds_status==1])
genesSW=names(qvalDxr)[qvalDxr<.05]
genesTx <- unique(dxr$groupID[p.adjust(dxr$pvalue,"BH")<.05])
mean(trueGenes%in%genesSW)
mean(trueGenes%in%genesTx)
mean(genesSW%in%genesTx)
length(genesTx)-length(genesSW) #difference in genes found
sum(trueGenes%in%genesTx)-sum(trueGenes%in%genesSW) #difference in true genes found
(length(genesTx)-length(genesSW)-(sum(trueGenes%in%genesTx)-sum(trueGenes%in%genesSW)))/(length(genesTx)-length(genesSW)) #FDR of additional genes

## gene-level FDR
1-mean(genesSW%in%trueGenes)

### characterize the transcripts found
txSW <- rownames(getSignificantTx(stageRObj))
txTx <- dxr$featureID[p.adjust(dxr$pvalue,method="BH")<.05]
txTx <- txTx[!is.na(txTx)]
mean(txSW%in%txTx)
trueTx <- truth_tx$transcript_id[truth_tx$transcript_ds_status==1]
sum(txSW%in%trueTx)
sum(txTx%in%trueTx)
length(txTx)-length(txSW) #difference in tx found
sum(trueTx%in%txTx)-sum(trueTx%in%txSW) #difference in true tx genes found
(length(txTx)-length(txSW)-(sum(trueTx%in%txTx)-sum(trueTx%in%txSW)))/(length(txTx)-length(txSW)) #FDR of additional tx

### explore simulations
sum(truth_gene$gene_ds_status==1) #nr of genes
sum(truth_tx$transcript_ds_status==1) #nr of tx
table(table(as.character(truth_tx$gene_id[truth_tx$transcript_ds_status==1]))) #nr of tx per gene

mean(genesSW%in%genesTx)

#proportion of truly spliced transcripts found
significantTx=dxr$featureID[p.adjust(dxr$pvalue,method="BH")<.05]
significantTx=significantTx[!is.na(significantTx)]
mean(truth_tx$transcript_id[truth_tx$transcript_ds_status==1]%in%significantTx)

#which differential proportions do we find as significant?
hist(truth$diff_IsoPct[truth$transcript_ds_status==1],breaks=seq(-90,90,by=10))
hist(truth_tx[truth_tx$transcript_id[truth_tx$transcript_ds_status==1][truth_tx$transcript_id[truth_tx$transcript_ds_status==1]%in%significantTx],"diff_IsoPct"],add=TRUE, col=alpha("green",.4),breaks=seq(-90,90,by=10))

#distribution of DS genes
plot(log(truth$expected_count+1),pch=".")
points(x=(1:nrow(truth))[truth$transcript_ds_status==1],y=log(truth$expected_count[truth$transcript_ds_status==1]+1),col=2,pch="o")

### gene-level analysis

cobra <- COBRAData(padj = data.frame(kallisto_dexseq = qvalDxr,
                                     row.names = names(qvalDxr),
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth_gene, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "gene_ds_status")
cobraplot <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot)


### transcript-level analysis

padjTxDexSeq <- p.adjust(dxr$pvalue,"BH")
cobra <- COBRAData(padj = data.frame(kallisto_dexseq = padjTxDexSeq,
                                     row.names = dxr$featureID,
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth_tx, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "transcript_ds_status")
cobraplot <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   colorscheme = c("blue", "red", "green"),
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot)


## Stage-wise vs transcript-level analysis on the transcript level

### ROC stage-wise testing on transcript level: loop
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
pvalDexSeq <- dxr$pvalue
names(pvalDexSeq) <- dxr$featureID
pvalDexSeq[is.na(pvalDexSeq)] <- 1
padjDexSeq <- p.adjust(pvalDexSeq,method="BH")
tprSW <- vector(length=length(pvalSeq))
fprSW <- vector(length=length(pvalSeq))
fdrSW <- vector(length=length(pvalSeq))
stageRObjSim <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)

for(id in 1:length(pvalSeq)){
    print(id)
    ### DEXSeq stage-wise
    significantGenes <- names(qvalDxr)[which(qvalDxr<=pvalSeq[id])]
    if(length(significantGenes)==0){
	tprSW[id]=0
	fprSW[id]=0
	fdrSW[id]=0
	next
    } else {
      stageRRes <- stageWiseAdjustment(stageRObjSim, method="dtu", alpha=pvalSeq[id])
      positives <- rownames(getSignificantTx(stageRRes))

    #fast manual way
    #alphaAdjusted <- pvalSeq[id]*length(significantGenes)/length(qvalDxr)
    #genesStageII <- dxr$groupID[dxr$groupID%in%significantGenes]
    #uniqueGenesStageII <- unique(genesStageII)
    #txStageII <- dxr$featureID[dxr$groupID%in%significantGenes]
    #pvalStageII <- dxr$pvalue[dxr$featureID%in%txStageII]
    #pvalGeneList <- list()
    #for(i in 1:length(uniqueGenesStageII)){
    #    idHlp <- which(genesStageII==uniqueGenesStageII[i])
    #    pvalHlp <- pvalStageII[idHlp]
    #    names(pvalHlp) <- txStageII[idHlp]
    #    pvalGeneList[[i]] <- pvalHlp
    #}
    #padjGeneList <- lapply(pvalGeneList,function(x) adjustShaffer(x))
    #padjGeneList <- lapply(pvalGeneList,function(x) p.adjust(x,method="holm"))    
    #padj <- unlist(padjGeneList)
    #positives <- names(padj[padj<=alphaAdjusted])

    tprSW[id] <- mean(truth_tx$transcript_id[truth_tx$transcript_ds_status==1]%in%positives)  
    fprSW[id] <- mean(truth_tx$transcript_id[truth_tx$transcript_ds_status==0]%in%positives)
    truePosId <- truth_tx$transcript_ds_status[match(positives,truth_tx$transcript_id)]
    fdrSW[id] <- 1-mean(truePosId)
    }
}

### DEXSeq regular
evalDexSeqRegular <- t(sapply(pvalSeq,function(alpha){
	   positives <- names(padjDexSeq[which(padjDexSeq<alpha)])
	   tprRegular <- mean(truth_tx$transcript_id[truth_tx$transcript_ds_status==1]%in%positives)	   
	   fprRegular <- mean(truth_tx$transcript_id[truth_tx$transcript_ds_status==0]%in%positives)	   
	   truePosId <- truth_tx$transcript_ds_status[match(positives,truth_tx$transcript_id)]	   
	   fdrRegular <- 1-mean(truePosId)
	   return(cbind(tprRegular=tprRegular, fprRegular=fprRegular, fdrRegular=fdrRegular))
				   }))
colnames(evalDexSeqRegular) <- c("tpr","fpr","fdr")


## ROC curve
plot(x=fprSW,y=tprSW,type="n", xlim=c(0,0.1), bty="l",xlab="False Positive Rate", ylab="True Positive Rate")
lines(x=fprSW,y=tprSW,col=2,lwd=2)
points(x=fprSW[516],y=tprSW[516],pch=3,col=2)
lines(x=evalDexSeqRegular[,"fpr"],y=evalDexSeqRegular[,"tpr"],col=3,lwd=2)
points(x=evalDexSeqRegular[516,"fpr"],y=evalDexSeqRegular[516,"tpr"],pch=3,col=3)
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3, bty="n")

## FDR-TPR curve
plot(x=fdrSW,y=tprSW,type="n", xlab="False Discovery Rate", ylab="True Positive Rate", bty="l")
lines(x=fdrSW,y=tprSW,col=2,lwd=2)
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],pch="o",col=2)
lines(x=evalDexSeqRegular[,"fdr"],y=evalDexSeqRegular[,"tpr"],col=3,lwd=2)
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],col=3,pch="o")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3, bty="n")



# Stage-wise vs transcript level analysis on the gene level

### stage-wise testing 
pvalSeq = c(1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
positiveGenes <- as.character(truth_gene$gene_id[truth_gene$gene_ds_status==1])
negativeGenes <- as.character(truth_gene$gene_id[truth_gene$gene_ds_status==0])
evalGeneLevelSW <- t(sapply(pvalSeq, function(alpha){
  significantGenes <- names(qvalDxr)[qvalDxr<=alpha]
  tpr <- mean(positiveGenes%in%significantGenes)
  fpr <- mean(negativeGenes%in%significantGenes)
  fdr <- mean(significantGenes%in%negativeGenes)
  return(cbind(tpr,fpr,fdr))
}))
colnames(evalGeneLevelSW) <- c("tpr","fpr","fdr")


evalGeneLevelTx <- t(sapply(pvalSeq, function(alpha){
  significantGenes <- unique(dxr$groupID[p.adjust(dxr$pvalue,"BH")<=alpha])
  tpr <- mean(positiveGenes%in%significantGenes)
  fpr <- mean(negativeGenes%in%significantGenes)
  fdr <- mean(significantGenes%in%negativeGenes)
  return(cbind(tpr,fpr,fdr))
}))
colnames(evalGeneLevelTx) <- c("tpr","fpr","fdr")


## ROC curve
plot(x=evalGeneLevelSW[,"fpr"],y=evalGeneLevelSW[,"tpr"],type="n", xlim=c(0,0.1), bty="l",xlab="False Positive Rate", ylab="True Positive Rate")
lines(x=evalGeneLevelSW[,"fpr"],y=evalGeneLevelSW[,"tpr"],col=2,lwd=2)
points(x=evalGeneLevelSW[516,"fpr"],y=evalGeneLevelSW[516,"tpr"],pch=3,col=2)
lines(x=evalGeneLevelTx[,"fpr"],y=evalGeneLevelTx[,"tpr"],col=3,lwd=2)
points(x=evalGeneLevelTx[516,"fpr"],y=evalGeneLevelTx[516,"tpr"],pch=3,col=3)
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3, bty="n")

## FDR-TPR curve
plot(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"],type="n", xlab="False Discovery Rate", ylab="True Positive Rate", bty="l")
lines(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"],col=2,lwd=2)
points(x=evalGeneLevelSW[c(508,516,526),"fdr"],y=evalGeneLevelSW[c(508,516,526),"tpr"],pch="o",col=2)
lines(x=evalGeneLevelTx[,"fdr"],y=evalGeneLevelTx[,"tpr"],col=3,lwd=2)
points(x=evalGeneLevelTx[c(508,516,526),"fdr"],y=evalGeneLevelTx[c(508,516,526),"tpr"],col=3,pch="o")
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
legend("topleft",c("Stage-wise","Regular"),lty=1,col=2:3, bty="n")


### combine gene level, tx level and stage-wise analysis in one FDR-TPR plot.
#gene level
par(mar=c(5,4,4,1)+0.2)
plot(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"], type="n", xlab="False Discovery Proportion", ylab="True Positive Rate", ylim=c(0.2,1), bty="l", main="Drosophila", cex.lab=1.5)
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.8),lty=2)
lines(x=evalGeneLevelSW[,"fdr"],y=evalGeneLevelSW[,"tpr"], col="black", lwd=2)
points(x=evalGeneLevelSW[c(508,516,526),"fdr"],y=evalGeneLevelSW[c(508,516,526),"tpr"],pch=19,col="white")
points(x=evalGeneLevelSW[c(508,516,526),"fdr"],y=evalGeneLevelSW[c(508,516,526),"tpr"],col="black")
#transcript level
lines(x=evalDexSeqRegular[,"fdr"],y=evalDexSeqRegular[,"tpr"],col="red",lwd=2)
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],pch=19,col="white")
points(x=evalDexSeqRegular[c(508,516,526),"fdr"],y=evalDexSeqRegular[c(508,516,526),"tpr"],col="red")
#stage-wise transcript level
lines(x=fdrSW,y=tprSW, lwd=2, col="green3")
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],pch=19,col="white")
points(x=fdrSW[c(508,516,526)],y=tprSW[c(508,516,526)],col="green3")
legend("bottomright",c("gene level","transcript level","transcript level stage-wise"),lty=1,col=c("black","red","green3"), bty="n", cex=1.25)

