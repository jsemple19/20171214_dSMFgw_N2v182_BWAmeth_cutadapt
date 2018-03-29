#######
# script calls methylation, and calculates fraction methylated at each C
# (correcting for issues such as duplicate counts in GCG context)
# updated: 2018-02-16
# author: Jennnifer Semple
# Scripts adapted from Drosophila scripts kindly provided by Arnaud Krebs

library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")
#modify the UCSC object to have same chr numbers as Wormbase genome fasta file
seqnames(Celegans)<-gsub("chr", "", seqnames(Celegans))
seqnames(Celegans)<-gsub("M", "MtDNA", seqnames(Celegans))
# also read is location of genome fasta file for QuasR objects
genomeFile<-"~/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"

# collect citations for packages used
packageBib<-toBibtex(c(citation("QuasR"),
                       citation("BSgenome.Celegans.UCSC.ce11")))

setwd("~/Documents/MeisterLab/sequencingData/20171214_dSMFgw_N2v182_BWAmeth")
source('./R/callAllCs.r') #to call all Cs
source('./R/useful_functionsV1.r') #load the ranges


#####    reload project using list of bam files

path='.'
my.alignmentsDir=paste(path,'/aln/',sep='')

#sp.list=read.delim( "./QuasR_Aligned.txt",sep='\t')  # delete? redundantly creares another file??
#write.table(sp.list,'./tmp/sample_BAM.tmp',sep='\t',row.names=FALSE)

cluObj=makeCluster(3)

#need to index .bam files. e.g samtools index -b  180126_SNK268_A_L001_JIB-4.filt.bam

NOMEproj=qAlign(sampleFile='./QuasR_Aligned.txt',
              genome=genomeFile,
              paired="fr",
              bisulfite="dir",
              projectName="dSMF_N2",
              clObj=cluObj)


NOMEaln=as.data.frame(alignments(NOMEproj)$genome)
samples=NOMEaln$SampleName

#### call methylating of Cs in data
meth_gr=qMeth(NOMEproj, mode='allC',clObj=cluObj)
# 35541247
#todayDate<-format(Sys.time(), "%Y%m%d")

if (!dir.exists(paste0(path,"./methylation_calls"))){
  dir.create(paste0(path,"./methylation_calls"))
}
if (!dir.exists(paste0(path,"./plots"))){
  dir.create(paste0(path,"./plots"))
}

## save as rds for future access
saveRDS(meth_gr,paste0(path,'/methylation_calls/NOME_allCs.rds'))

# and make some histograms
pdf("./plots/hist_C_coverage.pdf",width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal],breaks=100,xlim=c(1,5000),
       main=paste(s, ": total coverage"),xlab="read counts")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth],breaks=100,xlim=c(1,5000),
       main=paste(s, ": counts of methylated Cs"),xlab="read counts")
}
dev.off()


# find sequence context of Cs using function from callAllCs.r file
cO=10 # minimal read coverage for a C (low coverage discarded)
methFreq_grl=call_context_methylation(meth_gr,cO,genome=Celegans)

# call_context_methylation returns list of two matrices, "CG" and "GC" in which
# V1 column with fraction methylation and type column with C context
pdf("./plots/hist_CG-GC_freq.pdf",width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples){
  columnMeth<-paste0(s,"_M")
  hist(unlist(mcols(methFreq_grl[["CG"]])[,columnMeth]),breaks=100,
       main=paste(s,"CmG frequency"),xlab="fraction methylated")
  hist(unlist(mcols(methFreq_grl[["GC"]])[,columnMeth]),breaks=100,
       main=paste(s,"GCm frequency"),xlab="fraction methylated")
}
dev.off()

saveRDS(methFreq_grl,paste0(path,"/methylation_calls/NOME_CG-GC.rds"))

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



