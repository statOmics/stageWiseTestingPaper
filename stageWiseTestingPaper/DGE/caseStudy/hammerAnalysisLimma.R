setwd("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/")
load("/Users/koenvandenberge/Dropbox/PhD/seminars/presentations/Janssen/hammer_eset.RData")
eset = hammer.eset ; rm(hammer.eset)
eset@phenoData@data$Time #typo. Will do it ourself
time = factor(rep(c("mo2","w2"),each=4),levels=c("w2","mo2"))
eset@phenoData@data$protocol
treat = factor(c("control","control","SNL","SNL","control","control","SNL","SNL"),levels=c("control","SNL"))
design = model.matrix(~time*treat)
rownames(design) = paste0(time,treat,"_",techRep)
colnames(design)[4] = "timeMo2xTreatSNL"

library(edgeR) ; library(Biobase) ; library(limma)
cpmOffset=2
keep = rowSums(cpm(exprs(eset))>cpmOffset)>=2 #2cpm in 2 samples
d = DGEList(exprs(eset)[keep,])
colnames(d) = rownames(design)
d = calcNormFactors(d)


## QC
library(RColorBrewer)
palette(brewer.pal(8,"Dark2"))
plotMDS(d,col=as.numeric(treat),labels=colnames(d))
plot(density(cpm(d$counts[,1],log=TRUE)),type='n')
for(i in 1:ncol(d$counts)) lines(density(cpm(d$counts[,i],log=TRUE)),col=i)
plot(density(cpm(d$counts[,1],log=TRUE)),type='n')
for(i in 1:ncol(d$counts)) lines(density(cpm(d$counts[,i],log=TRUE)),col=as.numeric(treat)[i])

## DE analysis
L=matrix(0,ncol=3,nrow=ncol(design))
colnames(L)=c("SNL-C_w2","SNL-C_mo2","SNL-C_mo2-w2")
rownames(L)=colnames(design)
L["treatSNL",1]=1
L[grep(rownames(L),pattern="SNL"),2]=1
L[,3]=L[,2]-L[,1]


## regular analysis
v = voom(d,design,plot=TRUE)
fit=lmFit(v,design)
contrast.matrix <- makeContrasts(treatSNL, treatSNL+timeMo2xTreatSNL, timeMo2xTreatSNL, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res=decideTests(fit2)
summary.TestResults(res)
uniqueGenesRegular=which(res[,1]!=0 | res[,2]!=0 | res[,3]!=0)
length(uniqueGenesRegular) #nr of unique genes regular analysis


## stage-wise testing Shaffer
nGenes=nrow(d)
alpha=0.05
tableF = topTableF(fit2, number=nGenes, sort.by="none", p.value=0.05)
genesSI=rownames(tableF) #genes significant stage I
alphaAdjusted=alpha*length(genesSI)/nGenes

pValuesStageIIShaffer=sapply(1:3,function(i) topTable(fit2, coef=i, number=nGenes, sort.by="none")$P.Value)
colnames(pValuesStageIIShaffer)=colnames(fit2$coefficients)
rownames(pValuesStageIIShaffer)=rownames(d)
# set p-values of genes not passing Stage I to 1
pValuesStageIIShaffer[!(rownames(pValuesStageIIShaffer)%in%genesSI),]=1

contrastsAll=sapply(1:3,function(i) topTable(fit2, coef=i, number=nGenes, sort.by="none")$logFC)
rownames(contrastsAll)=rownames(fit2)
resultsStageIIShaffer=(pValuesStageIIShaffer<alphaAdjusted)*sign(contrastsAll)
summary.TestResults(resultsStageIIShaffer)
uniqueGenesSW=which(resultsStageIIShaffer[,1]!=0 | resultsStageIIShaffer[,2]!=0 | resultsStageIIShaffer[,3]!=0)
length(uniqueGenesSW)

genesNotFoundStageII <- genesSI[genesSI %in% rownames(resultsStageIIShaffer)[rowSums(resultsStageIIShaffer==0)==3]]
length(genesNotFoundStageII) #stage I only genes

significantFCInteractionLimma=contrastsAll[resultsStageIIShaffer[,3]!=0,3]
hist(significantFCInteractionLimma,50)
save(significantFCInteractionLimma,file="significantFCInteractionLimma.rda")

## none of the genes that were only found in the screening hypothesis are found in a regular analysis.
mean(rowSums(res[genesNotFoundStageII,]==0)==3)

## plots
#### expression patterns between groups of genes
par(mfrow=c(1,4))

index=which(abs(resultsStageIIShaffer[,1])==1)
labels=rownames(fit2)[index]
log2fc1 = contrastsAll[index,1]
log2fc2 = contrastsAll[index,2]
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="DE at t1",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=which(abs(resultsStageIIShaffer[,2])==1)
labels=rownames(fit2)[index]
log2fc1 = contrastsAll[index,1]
log2fc2 = contrastsAll[index,2]
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="DE at t2",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=which(abs(resultsStageIIShaffer[,3])==1)
labels=rownames(fit2)[index]
log2fc1 = contrastsAll[index,1]
log2fc2 = contrastsAll[index,2]
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="Significant interaction",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=genesNotFoundStageII
labels=index
log2fc1 = contrastsAll[index,1]
log2fc2 = contrastsAll[index,2]
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="Stage I only",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.4)



