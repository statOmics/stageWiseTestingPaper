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
#keep = rowSums(cpm(exprs(eset))>cpmOffset)>=4 #2cpm in 4 samples
keep = rowSums(cpm(exprs(eset))>cpmOffset)>=2 #2cpm in 2 samples
d = DGEList(exprs(eset)[keep,])
colnames(d) = rownames(design)
d = calcNormFactors(d)
boxplot(colSums(exprs(eset))~time) #may explain high number of DE genes within timepoint
boxplot(colSums(exprs(eset))~treat) #may explain low number of genes for interaction effect, because more genes are DE within timepoint due to the big library size differences, but these get cancelled out when testing the interaction effect.

## QC
library(RColorBrewer)
palette(brewer.pal(8,"Dark2"))
plotMDS(d,col=as.numeric(treat),labels=colnames(d))
plot(density(cpm(d$counts[,1],log=TRUE)),type='n')
for(i in 1:ncol(d$counts)) lines(density(cpm(d$counts[,i],log=TRUE)),col=i)
plot(density(cpm(d$counts[,1],log=TRUE)),type='n')
for(i in 1:ncol(d$counts)) lines(density(cpm(d$counts[,i],log=TRUE)),col=as.numeric(treat)[i]) #timepoint 2 has higher expression for quite some genes.

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
uniqueGenesRegular=which(rowSums(results==0)==3)
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
uniqueGenesSW=which(rowSums(resultsStageIIShaffer==0)==3)

genesNotFoundStageII <- genesSI[genesSI %in% rownames(resultsStageIIShaffer)[rowSums(resultsStageIIShaffer==0)==3]]
length(genesNotFoundStageII) #stage I only genes
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]!=0 & resultsStageIIShaffer[,3]==0) #5073 genes DE in both without interaction
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]==0) #1268 genes only in T1
sum(resultsStageIIShaffer[,1]==0 & resultsStageIIShaffer[,2]!=0) #942 genes only in T2
sum(resultsStageIIShaffer[,1]!=0 & resultsStageIIShaffer[,2]!=0 & resultsStageIIShaffer[,3]!=0) #411 genes DE in all three contrasts.
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


### fold changes for simulation
Lsim <- matrix(0,nrow=ncol(design),ncol=3)
rownames(Lsim) <- colnames(design)
colnames(Lsim) <- c("deT1","deT2","time") # these three already imply interaction effect. Also can be seen because rank of matrix does not change when you add interaction contrast.
Lsim[3,1]=1
Lsim[3:4,2]=1
Lsim[c(2,4),3]=c(1,1/2)
lrtSim=list()
for(i in 1:ncol(Lsim)) lrtSim[[i]] <- glmLRT(fit,contrast=Lsim[,i])
lrtSim[[i+1]]=glmLRT(fit,contrast=Lsim)
contrastsSim=list()
contrastsSim[[1]] <- sapply(1:ncol(Lsim), function(i) decideTestsDGE(lrt[[i]],p.value=alpha))
contrastsSim[[2]] <- sapply(1:ncol(Lsim), function(i) lrt[[i]]$table$logFC)
save(contrastsSim,file="contrastsSim.RData")



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








library(gplots)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
significance=c("non-sign.","sign.")
sig1=1
sig2=1
inter=1
index = which(abs(resultsStageIIShaffer[,1])==sig1 & abs(resultsStageIIShaffer[,2])==sig2 & abs(resultsStageIIShaffer[,3])==inter)
labels=fit$genes[index,1]
heatmap.2(cpm(exprs(eset))[index,]-rowMeans(cpm(exprs(eset))[index,]), labRow=labels, labCol=rownames(design), cexRow=.75, cexCol=.75, col=hmcol,main=paste(significance[sig1+1],"FC10,",significance[sig2+1],"FC48\n",significance[inter+1],"interaction\n"),Colv=FALSE,trace="none",dendrogram="row")


## interaction genes classic
index=which(abs(results[,3])==1)
labels=fit$genes[index,1]
pdf("./presentation/interactionClassic.pdf",width=7,height=7)
heatmap.2(cpm(exprs(eset))[index,]-rowMeans(cpm(exprs(eset))[index,]), labRow=labels, labCol=rownames(design), cexRow=1.2, cexCol=1.2, col=hmcol,main="      interaction genes FDR",Colv=FALSE,trace="none",dendrogram="row", margins=c(8,8))
dev.off()

## interaction genes Shaffer
index=which(abs(resultsStageIIShaffer[,3])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="interactionGenes",ylim=range(c(log2fc1,log2fc2)),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)

pdf("./presentation/interactionShaffer.pdf",width=7,height=7)
heatmap.2(cpm(exprs(eset))[index,]-rowMeans(cpm(exprs(eset))[index,]), labCol=rownames(design), cexCol=1.2, col=hmcol,main="         interaction genes stage-wise Shaffer",Colv=FALSE,trace="none",dendrogram="row", margins=c(8,8))
dev.off()


## genes DE at t1 and interaction and not t2
par(mar=c(5.1,4.5,4.1,2.1))
index=which(abs(resultsStageIIShaffer[,1])==1 & abs(resultsStageIIShaffer[,3])==1 &  abs(resultsStageIIShaffer[,2])==0)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="",ylim=range(contrastsAll[,1:2]),xaxt="n",cex.axis=2,cex.lab=2)
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=2)
abline(h=0,lty=2,col="red",lwd=2)


## genes DE at t2 and interaction and not t1
index=which(abs(resultsStageIIShaffer[,1])==0 & abs(resultsStageIIShaffer[,3])==1 &  abs(resultsStageIIShaffer[,2])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="",ylim=range(contrastsAll[,1:2]),xaxt="n",cex.axis=2,cex.lab=2)
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=2)
abline(h=0,lty=2,lwd=2,col="red")


## genes DE at t2 and t1 and interaction
index=which(abs(resultsStageIIShaffer[,1])==1 & abs(resultsStageIIShaffer[,3])==1 &  abs(resultsStageIIShaffer[,2])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="",ylim=range(contrastsAll[,1:2]),xaxt="n",cex.axis=2,cex.lab=2)
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=2)
abline(h=0,lty=2,lwd=2,col="red")


### genes only stage I significant
matplot(t(contrastsAll[genesNotFoundStageII,1:2]),type="l",lty=1,col=1,ylab="log2(FC)",main="",ylim=range(contrastsAll[,1:2]),xaxt="n",cex.axis=2,cex.lab=2)
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=2)
abline(h=0,lty=2,lwd=2,col="red")


### plot genes only found by FDR method in timepoint1: should be genes with time-specific effects of t1 ==> you do not find them with stage-wise procedure because of dilution effect
## are also non-significant at t2. Else they would have been picked up by screening hypothesis..
index <- which(abs(results[,1])==1 & abs(resultsStageIIShaffer[,1])==0)
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="",ylim=range(contrastsAll[,1:2]),xaxt="n",cex.axis=2,cex.lab=2)
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=2)
abline(h=0,lty=2,lwd=2,col="red")







lrt[[3]]$table[index,"logFC"]

index=which(abs(results[,3])==1)
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="interactionGenes",ylim=range(c(log2fc1,log2fc2)),xaxt="n")








## average FC genes Shaffer
index=genesNotFoundStageII
labels=fit$genes[index,1]
logfc1 = fit$coef[index,3]
logfc2 = fit$coef[index,3] + fit$coef[index,4]
log2fc1 = logfc1*log2(exp(1))
log2fc2 = logfc2*log2(exp(1))
pdf("./presentation/genesNotFoundStage2.pdf")
matplot(t(cbind(log2fc1,log2fc2)),type="l",lty=1,col=1,ylab="log2(FC)",main="genes not found stage II",ylim=range(c(log2fc1,log2fc2)),xaxt="n")
axis(1,at=c(1,2),labels=c("timePoint1","timePoint2"),cex.axis=1.5)
dev.off()

## p-values?

## for simulation
##FC time effect
LavgTime = matrix(0,nrow=ncol(design),ncol=1)
rownames(LavgTime) = colnames(design)
LavgTime[c("timemo2","treatSNL","timemo2:treatSNL"),] = 1/2
hlp = glmLRT(fit,contrast=LavgTime)
signP = p.adjust(hlp$table$PValue,method="BH")<.05
fcOfSignGenesTimeEffect = hlp$table$logFC[signP] #this is on log2 scale
fcOfSignGenesTimeEffect = exp(fcOfSignGenesTimeEffect*log(2)) #FC scale
pvalsOfSignGenesTimeEffect = hlp$table$PValue[signP]
plot(x=log2(fcOfSignGenesTimeEffect), y=-log10(pvalsOfSignGenesTimeEffect)) #vulcano plot
save(fcOfSignGenesTimeEffect,file="fcOfSignGenesTimeEffect.RData")

## FC t1
signP = p.adjust(lrt[[1]]$table$PValue,method="BH")<.05
fcOfSignGenesT1 = lrt[[1]]$table$logFC[signP] #this is on log2 scale
fcOfSignGenesT1 = exp(fcOfSignGenesT1*log(2)) #FC scale
pvalsOfSignGenesT1 = lrt[[1]]$table$PValue[signP]
plot(x=log2(fcOfSignGenesT1), y=-log10(pvalsOfSignGenesT1)) #vulcano plot
save(fcOfSignGenesT1,file="fcOfSignGenesT1.RData")

## FC t2
signP = p.adjust(lrt[[2]]$table$PValue,method="BH")<.05
fcOfSignGenesT2 = lrt[[2]]$table$logFC[signP] #this is on log2 scale
fcOfSignGenesT2 = exp(fcOfSignGenesT2*log(2)) #FC scale
pvalsOfSignGenesT2 = lrt[[2]]$table$PValue[signP]
plot(x=log2(fcOfSignGenesT2), y=-log10(pvalsOfSignGenesT2)) #vulcano plot
save(fcOfSignGenesT2,file="fcOfSignGenesT2.RData")



