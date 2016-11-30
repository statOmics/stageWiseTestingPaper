setwd("/Users/koenvandenberge/Dropbox/PhD/Research/stageWiseTesting/")
load("/Users/koenvandenberge/Dropbox/PhD/seminars/presentations/Janssen/hammer_eset.RData")
eset = hammer.eset ; rm(hammer.eset)
eset@phenoData@data$Time #typo. Will do it ourself
time = factor(rep(c("mo2","w2"),each=4),levels=c("w2","mo2"))
eset@phenoData@data$protocol
treat = factor(c("control","control","SNL","SNL","control","control","SNL","SNL"),levels=c("control","SNL"))
techRep = eset@phenoData@data$num.tech.reps #this tells you how many technical replicates were pooled to obtain the counts
design = model.matrix(~time*treat)
rownames(design) = paste0(time,treat,"_",techRep)

library(edgeR) ; library(Biobase)
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
d = estimateDisp(d,design)
plotBCV(d)
fit = glmFit(d,design)
gof(fit,plot=TRUE)

L=matrix(0,ncol=3,nrow=ncol(design))
colnames(L)=c("SNL-C_w2","SNL-C_mo2","SNL-C_mo2-w2")
rownames(L)=colnames(design)
L["treatSNL",1]=1
L[grep(rownames(L),pattern="SNL"),2]=1
L[,3]=L[,2]-L[,1]

###regular analysis
alpha=0.05
lrt=list()
for (i in 1:ncol(L)) lrt[[i]]=glmLRT(fit,contrast=L[,i]) #each contrast seperately
lrt[[i+1]]=glmLRT(fit,contrast=L) #F-test over all contrasts
results=sapply(lrt[1:ncol(L)],function(x) decideTestsDGE(x,p.value=alpha)) #alpha
rownames(results) <- rownames(lrt[[1]])
names(results)=colnames(L)
summary.TestResults(results)
# unique genes
uniqueGenesRegular=which(results[,1]!=0 | results[,2]!=0 | results[,3]!=0)
length(uniqueGenesRegular) #nr of unique genes regular analysis

### Stagewise testing: Shaffer
nGenes=nrow(d)
alpha=0.05
tableF=topTags(lrt[[4]],n=nGenes,p.value=alpha)
genesSI=rownames(tableF) #genes significant stage I
alphaAdjusted=alpha*length(genesSI)/nGenes

pValuesStageIIShaffer=sapply(lrt[1:3],function(x) x$table$PValue)
colnames(pValuesStageIIShaffer)=colnames(L)
rownames(pValuesStageIIShaffer)=rownames(d)
# set p-values of genes not passing Stage I to 1
pValuesStageIIShaffer[!(rownames(pValuesStageIIShaffer)%in%genesSI),]=1

contrastsAll=sapply(lrt[1:3],function(x) x$table$logFC )
resultsStageIIShaffer=(pValuesStageIIShaffer<alphaAdjusted)*sign(contrastsAll)
summary.TestResults(resultsStageIIShaffer)
uniqueGenesSW=which(resultsStageIIShaffer[,1]!=0 | resultsStageIIShaffer[,2]!=0 | resultsStageIIShaffer[,3]!=0)
length(uniqueGenesSW)

genesNotFoundStageII <- genesSI[genesSI %in% rownames(resultsStageIIShaffer)[rowSums(resultsStageIIShaffer==0)==3]]
length(genesNotFoundStageII) #stage I only genes
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]!=0 & resultsStageIIShaffer[,3]==0) #genes DE in both without interaction
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]==0) #genes only in T1
sum(resultsStageIIShaffer[,1]==0 & resultsStageIIShaffer[,2]!=0) #genes only in T2
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]!=0 & resultsStageIIShaffer[,3]!=0) #genes DE in all three contrasts.
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]==0 & resultsStageIIShaffer[,3]!=0) #
sum(resultsStageIIShaffer[,1]==0 & resultsStageIIShaffer[,2]!=0 & resultsStageIIShaffer[,3]!=0) #
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]==0 & resultsStageIIShaffer[,3]==0) #


## stagewise testing: Holm correction in stage II
holmAdjP = t(apply(pValuesStageIIShaffer,1,function(x) p.adjust(x,method="holm")))
resultsStageIIHolm = (holmAdjP<alphaAdjusted)*sign(contrastsAll)
summary.TestResults(resultsStageIIHolm)


### only stage I genes and average FC
Lavg = cbind(L,0)
colnames(Lavg)[4]='averageFC'
Lavg[3,4]=1
Lavg[4,4]=.5
lrtavg = lrt[1:3]
lrtavg[[4]]=glmLRT(fit,contrast=Lavg[,4])
lrtavg[[5]]=glmLRT(fit,contrast=Lavg)

pValuesStageIIShafferAvg=sapply(lrtavg[1:4],function(x) x$table$PValue)
colnames(pValuesStageIIShafferAvg)=colnames(Lavg)
rownames(pValuesStageIIShafferAvg)=rownames(d)
#set p-values of genes not passing screening hypothesis to 1
pValuesStageIIShafferAvg[!(rownames(pValuesStageIIShafferAvg)%in%genesSI),]=1

contrastsAllAvg=sapply(lrtavg[1:4],function(x) x$table$logFC)
resultsStageIIShafferAvg=(pValuesStageIIShafferAvg<alphaAdjusted)*sign(contrastsAllAvg)
summary.TestResults(resultsStageIIShafferAvg)
#Assessing the genes that were not significant in the previous stage II test without the average test.
colMeans(abs(resultsStageIIShafferAvg[genesNotFoundStageII,]))
#all genes only found in screening stage are significant when testing for an average FC.

## none of the genes that were only found in the screening hypothesis are found in a regular analysis.
mean(rowSums(results[genesNotFoundStageII,]==0)==3)


### plots
#### expression patterns between groups of genes
par(mfrow=c(1,4))

index=which(abs(resultsStageIIShaffer[,1])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="DE at t1",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=which(abs(resultsStageIIShaffer[,2])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="DE at t2",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=which(abs(resultsStageIIShaffer[,3])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="Significant interaction",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

index=genesNotFoundStageII
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="Stage I only",ylim=c(-5,8),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)




