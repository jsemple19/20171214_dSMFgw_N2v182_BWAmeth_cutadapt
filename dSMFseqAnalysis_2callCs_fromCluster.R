#######
# script calls methylation, and calculates fraction methylated at each C
# (correcting for issues such as duplicate counts in GCG context)
# updated: 2017-06-28
# author: Jennnifer Semple
# Scripts adapted from Drosophila scripts kindly provided by Arnaud Krebs

library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")
# collect citations for packages used
packageBib<-toBibtex(c(citation("QuasR"),
                       citation("BSgenome.Celegans.UCSC.ce11")))

setwd("~/Documents/MeisterLab/sequencingData/20171214_dSMFgw_N2v182_BWAmeth_cutadapt/")

source('./R/callAllCs.r') #to call all Cs
source('./R/useful_functionsV1.r') #load the ranges


#####    reload project using list of bam files

path='./fromCluster'
my.alignmentsDir=paste(path,'/aln/',sep='')

#sp.list=read.delim( "./QuasR_Aligned.txt",sep='\t')  # delete? redundantly creares another file??
#write.table(sp.list,'./tmp/sample_BAM.tmp',sep='\t',row.names=FALSE)

#change paths if QuasR_Aligned.txt was created on cluster
aligned<-read.table(paste0(path,'/QuasR_Aligned.txt'),stringsAsFactors=F,header=T)
aligned$FileName<-gsub("/scratch/cluster/monthly/jsemple/20171214_dSMFgw_N2v182_BWAmeth_cutadapt","~/Documents/MeisterLab/sequencingData/20171214_dSMFgw_N2v182_BWAmeth_cutadapt/fromCluster",aligned$FileName)
write.table(aligned,paste0(path,"/QuasR_Aligned.txt"),quote=F,col.names=T,row.names=F,sep='\t',append=F)

cluObj=makeCluster(3)

NOMEproj=qAlign(sampleFile=paste0(path,'/QuasR_Aligned.txt'),
              genome="BSgenome.Celegans.UCSC.ce11",
              paired="fr",
              bisulfite="dir",
              projectName="dSMF_gw_N2vF2",
              clObj=cluObj)


NOMEaln=as.data.frame(alignments(NOMEproj)$genome)
samples=NOMEaln$SampleName



#### call methylating of Cs in data
#meth_gr=qMeth(NOMEproj, mode='allC',clObj=cluObj)
# 35541247
#todayDate<-format(Sys.time(), "%Y%m%d")

if (!dir.exists(paste0(path,"/methylation_calls"))){
  dir.create(paste0(path,"/methylation_calls"))
}
if (!dir.exists(paste0(path,"/docs"))){
  dir.create(paste0(path,"/docs"))
}

## save as rds for future access
#saveRDS(meth_gr,paste0(path,'/methylation_calls/NOME_allCs.rds'))

meth_gr<-readRDS(paste0(path,'/methylation_calls/NOME_allC.rds'))

# and make some histograms
pdf(paste0(path,"/plots/hist_C_coverage.pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal],breaks=10000,xlim=c(1,200),
       main=paste(s, ": total coverage"),xlab="read counts")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth],breaks=10000,xlim=c(1,200),
       main=paste(s, ": counts of methylated Cs"),xlab="read counts")
}
dev.off()

#table of coverage. Cs0freq is number of cytosines with no coverage. rest of data is based on Cs with coverage
Ccoverage<-data.frame(sampleNames=samples,Cs0freq=0,meanCoverage=0,medianCoverage=0,stdevCoverage=0,
                      meanMethCount=0,medianMethCount=0,stdevMethCount=0)
#excluding 0s
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  Ccoverage[Ccoverage$sampleNames==s,"Cs0freq"]<-sum(mcols(meth_gr)[,columnTotal]==0)/length(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanCoverage"]<-mean(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"medianCoverage"]<-median(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"stdevCoverage"]<-sd(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanMethCount"]<-mean(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"medianMethCount"]<-median(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"stdevMethCount"]<-sd(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
}

# sampleNames    Cs0freq meanCoverage medianCoverage stdevCoverage meanMethCount medianMethCount stdevMethCount
# 1 N2_DE_gwV006 0.13088379     20.17185              8      81.23431      2.874012               0       32.53755
# 2 N2_DE_gwV007 0.02866742     21.56688             15      71.41639      3.052619               0       28.69126
# 3 F2_DE_gwV008 0.04534261     11.95661              8      54.29670      1.884294               0       22.12124
# 4 F2_DE_gwV009 0.02661173     23.28572             17      78.62135      3.163520               0       30.95672

write.csv(Ccoverage, file=paste0(path,"/docs/CytosineCoverage_excl0.csv"))


# not excluding 0s
Ccoverage<-data.frame(sampleNames=samples,Cs0freq=0,meanCoverage=0,medianCoverage=0,stdevCoverage=0,
                      meanMethCount=0,medianMethCount=0,stdevMethCount=0)

for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  Ccoverage[Ccoverage$sampleNames==s,"Cs0freq"]<-sum(mcols(meth_gr)[,columnTotal]==0)/length(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanCoverage"]<-mean(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"medianCoverage"]<-median(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"stdevCoverage"]<-sd(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanMethCount"]<-mean(mcols(meth_gr)[,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"medianMethCount"]<-median(mcols(meth_gr)[,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"stdevMethCount"]<-sd(mcols(meth_gr)[,columnMeth])
}


# sampleNames    Cs0freq meanCoverage medianCoverage stdevCoverage meanMethCount medianMethCount stdevMethCount
# 3 N2_DE_gwV006 0.13088379     17.53168              6      76.03681      2.497851               0       30.34907
# 4 N2_DE_gwV007 0.02866742     20.94862             14      70.47724      2.965108               0       28.28160
# 1 F2_DE_gwV008 0.04534261     11.41447              7      53.10974      1.798855               0       21.61746
# 2 F2_DE_gwV009 0.02661173     22.66605             16      77.65865      3.079333               0       30.54628

write.csv(Ccoverage, file=paste0(path,"/docs/CytosineCoverage.csv"))



# find sequence context of Cs using function from callAllCs.r file
cO=10 # minimal read coverage for a C (low coverage discarded)
#methFreq_grl=call_context_methylation(meth_gr,cO,genome=Celegans)
methFreq_grl<-readRDS(paste0(path,"/methylation_calls/NOME_CG-GC.rds"))


# call_context_methylation returns list of two matrices, "CG" and "GC" in which
# V1 column with fraction methylation and type column with C context
pdf(paste0(path,"/plots/hist_CG-GC_freq.pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples){
  columnMeth<-paste0(s,"_M")
  hist(unlist(mcols(methFreq_grl[["CG"]])[,columnMeth]),breaks=100,
       main=paste(s,"CmG frequency"),xlab="fraction methylated")
  hist(unlist(mcols(methFreq_grl[["GC"]])[,columnMeth]),breaks=100,
       main=paste(s,"GCm frequency"),xlab="fraction methylated")
}
dev.off()

#saveRDS(methFreq_grl,paste0(path,"/methylation_calls/NOME_CG-GC.rds"))

################## conmbine CG and GC lists ########################
# concatenate both lists
allCs<-c(methFreq_grl$CG,methFreq_grl$GC)

mergeCGGClists<-function(grl) {
  #reduce CG gr to width 1 bp
  end(grl$CG)<-start(grl$CG)
  #reduce GC gr to width 1 bp
  start(grl$GC)<-end(grl$GC)

  #combine lists, leaving out duplicate GCGs
  allCs<-c(grl$CG,grl$GC[mcols(grl$GC)$type=="GCH"])
  #sanity checks - should get same as length of allCs
  #ol<-findOverlaps(grl$CG,grl$GC)
  #estimate1<-table(mcols(grl$CG)$type)[1]+table(mcols(grl$GC)$type)[2]+table(mcols(grl$GC)$type)[1]
  #estimate2<-length(grl$CG)+length(grl$GC)-length(ol)
  allCs<-sort(allCs)
  return(allCs)
}

allCs<-mergeCGGClists(methFreq_grl)

saveRDS(allCs,paste0(path,"/methylation_calls/NOME_freq_allCs.rds"))

# #######################################
# #######################################
# # now lets call qMeth only on the amplicons
#
# #meth_gr<-readRDS(paste0(path,"/methylation_calls/NomeAmpliconCEmeth.rds"))
# #3131525 & 3340523
#
# # import data about TSS positions
# TSS=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/ampliconTSS.bed')
# # single bp GR for genes with amplicon
#
# # resize TSS to +-300
# TSS_600=resize(TSS,600,fix='center')
# strand(TSS_600)='+'
#
# meth_gr <- qMeth(NOMEproj,mode="allC",TSS_600)
# #21532
# saveRDS(meth_gr,paste0(path,'rds/Amplicon_raw_methylation.rds'))
#
# hist(meth_gr$N2_DE_ampV001_T,breaks=50)
# hist(meth_gr$N2_DE_ampV001_M,breaks=50)
#
#
# # get % methylation and context in separate CG and GC matrices
# cO=10 #minimum coverage of site
# PolII_in=call_context_methylation(meth_gr,cO,Celegans)
# #1994 & 2001
#
#
# # merge the CG and GC matrices
# PolII_inm=unlist(GRangesList(PolII_in))
# #3995
#
# hist(PolII_inm$V1,breaks=100)
# # hardly any that are completely unmethylated. some are completely methylated
#
# saveRDS(PolII_inm,paste0(path,'rds/PolII_inhib_av10.rds'))
#
# #filter out NAs
# naix=apply(as.matrix(elementMetadata(PolII_inm[,1])),1,function(x){sum(is.na(x))==1})
#
# PolII_inm<-PolII_inm[!naix]
# #2636
#
# saveRDS(PolII_inm,paste0(path,'rds/PolII_inhib_av10_naix.rds'))
#
# ####################
# ####################
# # check overlap with amplicons. plot with package functions.
# PolII_inm=readRDS(paste0(path,'rds/PolII_inhib_av10.rds'))
#
# #filter out NAs
# naix=apply(as.matrix(elementMetadata(PolII_inm[,1])),1,function(x){sum(is.na(x))==1})
#
# #look hom many amplicons are covered
#
# amplicons=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/amplicons_stranded.bed')
# #amplicons with the strand of the bisulfite DNA amplified in the amplicon
#
# table(countOverlaps(amplicons,PolII_inm[!naix,],ignore.strand=T)>0)
# #FALSE  TRUE
# #16    80
#
# #AmpOL<-findOverlaps(amplicons,PolII_inm[!naix,],ignore.strand=T)
#
#
# pdf("./plots/AmpliconViews_AmpliconBiSeq.pdf",paper="a4",height=11,width=8)
# for (i in 1:length(amplicons)) {
#    b<-AmpliconViews(NOMEproj,amplicons[i],sampleNames="N2_DE_ampV001")
#    try(plotAmpliconView(b),silent=TRUE)
# }
# dev.off()
#
# ##### dir and UD are quite different. Not in % meth but in the output reads.
#
#
# # to get individual methylation matrices for an amplicon use getCMethMatrix function
# # from callAllCs.r file
# #a<-getCMethMatrix(NOMEproj,amplicons[1],"N2")
# #aUD<-getCMethMatrix(NOMEprojUD,amplicons[1],"N2")
#
#
#
