load("~/performanceDataDrosDTE.rda") #from DTE analysis script
load("~/performanceDataHuman.rda") #from DTE analysis script

par(mfrow=c(1,2))
par(bty="l", mar=c(5,4.5,4,1))
plot(x=performanceDataDrosophila$fdrQval, y=performanceDataDrosophila$tprQval, type="l", ylab="True Positive Rate", xlab="False Discovery Proportion", ylim=c(0,1), lwd=2, bty="l", main="", cex.lab=1.5)
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceDataDrosophila$fdrQval[c(508,516,526)],y=performanceDataDrosophila$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceDataDrosophila$fdrQval[c(508,516,526)],y=performanceDataDrosophila$tprQval[c(508,516,526)],col="black")
lines(x=performanceDataDrosophila$fdrTx, y=performanceDataDrosophila$tprTxTx, col=2, lwd=2)
points(x=performanceDataDrosophila$fdrTx[c(508,516,526)],y=performanceDataDrosophila$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceDataDrosophila$fdrTx[c(508,516,526)],y=performanceDataDrosophila$tprTxTx[c(508,516,526)],col=2)
lines(x=performanceDataDrosophila$fdrTxSW, y=performanceDataDrosophila$tprTxSW, col=3, lwd=2)
points(x=performanceDataDrosophila$fdrTxSW[c(508,516,526)],y=performanceDataDrosophila$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceDataDrosophila$fdrTxSW[c(508,516,526)],y=performanceDataDrosophila$tprTxSW[c(508,516,526)],col=3)
legend("bottomright",c("gene-level","tx-level","tx-level stage-wise"),lty=1,col=1:3, bty="n", cex=1.25)

plot(x=performanceDataHuman$fdrQval, y=performanceDataHuman$tprQval, type="l", ylab="True Positive Rate", xlab="False Discovery Proportion", ylim=c(0,1), lwd=2, main="", cex.lab=1.5)
abline(v=c(.01,.05,seq(.1,.9,.1)),col=alpha("grey",.5),lty=2)
points(x=performanceDataHuman$fdrQval[c(508,516,526)],y=performanceDataHuman$tprQval[c(508,516,526)], pch=19,col="white")
points(x=performanceDataHuman$fdrQval[c(508,516,526)],y=performanceDataHuman$tprQval[c(508,516,526)],col="black")
lines(x=performanceDataHuman$fdrTx, y=performanceDataHuman$tprTxTx, col=2, lwd=2)
points(x=performanceDataHuman$fdrTx[c(508,516,526)],y=performanceDataHuman$tprTxTx[c(508,516,526)], pch=19,col="white")
points(x=performanceDataHuman$fdrTx[c(508,516,526)],y=performanceDataHuman$tprTxTx[c(508,516,526)],col=2)
lines(x=performanceDataHuman$fdrTxSW, y=performanceDataHuman$tprTxSW, col=3, lwd=2)
points(x=performanceDataHuman$fdrTxSW[c(508,516,526)],y=performanceDataHuman$tprTxSW[c(508,516,526)], pch=19,col="white")
points(x=performanceDataHuman$fdrTxSW[c(508,516,526)],y=performanceDataHuman$tprTxSW[c(508,516,526)],col=3)
legend("bottomright",c("gene-level","tx-level","tx-level stage-wise"),lty=1,col=1:3, cex=1.25, bty="n")

