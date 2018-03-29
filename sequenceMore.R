#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

bname=args[1]
print(bname)
myPath=args[2] #inside aln/ directory
print(myPath)

pdf(file=paste0(myPath,"/../plots/moreSeqProj_",bname,".pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))


# c curve is used to compute the expected complexity curve of a mapped read file with a hypergeometric formula
cCurve<-read.table(paste0(myPath,"/c_curve_output_",bname,".txt"),stringsAsFactors=F,header=T)
plot(cCurve,main="Library complexity",col="dark green",type='l',lwd=2)
abline(a=0,b=1,lty=2)
legend("bottomright", legend=c("real data", "ideal data"), lty=c(1,2),col=c("dark green","black"),lwd=c(2,1))

#lc extrap is used to generate the expected yield for theoretical larger experiments and bounds on the
#number of distinct reads in the library and the associated confidence intervals
lcExtrap<-read.table(paste0(myPath,"/lc_extrap_output_",bname,".txt"),stringsAsFactors=F,header=T)

with(lcExtrap, plot(TOTAL_READS, EXPECTED_DISTINCT,main="Extrapolated yield",type='l',lwd=2,xlim=c(0,1.5e+09)))
lines(lcExtrap$TOTAL_READS, lcExtrap$LOWER_0.95CI,type='l',lty=2,col="grey")
lines(lcExtrap$TOTAL_READS, lcExtrap$UPPER_0.95CI,type='l',lty=2,col="grey")
polygon(c(lcExtrap$TOTAL_READS,rev(lcExtrap$TOTAL_READS)),c(lcExtrap$LOWER_0.95CI,rev(lcExtrap$UPPER_0.95CI)),
        col="light grey",lty=1,border="grey")
lines(lcExtrap$TOTAL_READS, lcExtrap$EXPECTED_DISTINCT,type='l',lwd=2)
abline(v=max(cCurve$total_reads),col="red")
abline(v=max(cCurve$total_reads)*2,col="red",lty=2)
abline(a=0,b=1,lty=2)


#bound pop is a method for estimating species richness, the total number of species or classes in the sampled population
boundPop<-read.table(paste0(myPath,"/bound_pop_output_",bname,".txt"),stringsAsFactors=F,header=T)
#abline(h=boundPop$log_mean_estimated_unobs,col="blue")
readCounts<-read.table(paste0(myPath,"/readCounts_",bname,".txt"),header=F,stringsAsFactors=F)
title(sub=paste("Total reads:",formatC(readCounts$V1,big.mark=",",format="d"),
                " Estimated Population size:",formatC(boundPop$log_mean_estimated_unobs,big.mark=",",format="d")))
legend("bottomright", legend=c("current reads", "double reads", "ideal data", "real data +95%CI"), lty=c(1,2,2,1),
       col=c("red","red","black","black"),lwd=c(1,1,1,2))


dev.off()
