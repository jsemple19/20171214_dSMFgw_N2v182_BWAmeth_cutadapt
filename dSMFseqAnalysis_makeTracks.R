# plot methFreq as bigWig
library(rtracklayer)
library(Biostrings)
library(zoo)
setwd("/Users/semple/Documents/MeisterLab/sequencingData/20171214_dSMFgw_N2v182_BWAmeth_cutadapt/")
genomeFile<-"../../GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"
path="./fromCluster"

if (!dir.exists(paste0(path,"/bed"))){
  dir.create(paste0(path,"/bed"))
}

# reading in individual CG and GC tracks
meth_grl<-readRDS(paste0(path,"/methylation_calls/NOME_CG-GC.rds"))

cgMeth<-meth_grl[["CG"]]
#3131525
gcMeth<-meth_grl[["GC"]]
#3340523

# reading in combined GC GC list
meth_all<-readRDS(paste0(path,"/methylation_calls/NOME_freq_allCs.rds"))
# 5918563

sampleNames<-names(mcols(meth_all))[grep("_M", names(mcols(meth_all)))]

gr2bed<-function(grObj,sampleNames,bedPrefix="") {
  for (n in sampleNames) {
    df.bed<-data.frame( seqnames=seqnames(grObj),
                        starts=start(grObj)-1,ends=end(grObj)-1,
                        names=c(rep(".",length(grObj))),
                        scores=mcols(grObj)[n],
                        strands=strand(grObj))
    write.table(df.bed,file=paste0(bedPrefix,n,".bed"),quote=F,sep="\t",row.names=F,col.names=F)
  }
}

gr2bed(cgMeth,sampleNames,paste0(path,"/bed/CG_"))
gr2bed(gcMeth,sampleNames,paste0(path,"/bed/GC_"))
gr2bed(meth_all,sampleNames,paste0(path,"/bed/allCs_"))

gr2bedGraph<-function(grObj,sampleNames,bedPrefix="",oneMinusScore=F) {
  for (n in sampleNames) {
    idx=which(is.na(mcols(grObj)[,n]))
    if (oneMinusScore==T) {  # to convert to %^SMF from %methylation
      scores=1-mcols(grObj[-idx])[,n]
    } else {
      scores=mcols(grObj[-idx])[,n]
    }
    df.bed<-data.frame( seqnames=seqnames(grObj[-idx]),
                        starts=as.integer(start(grObj[-idx])-1),ends=as.integer(end(grObj[-idx])-1),
                        scores=scores)
    fileConn<-file(paste0(bedPrefix,n,".bedGraph"))
    writeLines(paste0("track type=bedGraph name=",n),fileConn)
    write.table(df.bed,file=fileConn,append=T,quote=F,sep="\t",row.names=F,col.names=F)
  }
}

#########################################
### raw data tracks with no smoothing ###
#########################################

gr2bedGraph(cgMeth,sampleNames,paste0(path,"/bed/CG_"),oneMinusScore=T)
gr2bedGraph(gcMeth,sampleNames,paste0(path,"/bed/GC_"),oneMinusScore=T)
gr2bedGraph(meth_all,sampleNames,paste0(path,"/bed/allCs_"),oneMinusScore=T)

# create chrom sizes file needed by bedtobigwig command
#genome<-readDNAStringSet(genomeFile)
#df.chrom.sizes<- data.frame(chrNames=names(genome),chrSize=width(genome))
#write.table(df.chrom.sizes, paste0(path,"/bed/Ce11.chrom.sizes"),quote=F,sep="\t",row.names=F,col.names=F)

chromSizesFile<-"~/Documents/MeisterLab/GenomeVer/chromSizes/ws235.chrom.sizes"
#unix command to sort bedGraph files (required for bedtobigwig command)
sortCmd<-"LC_COLATE=C sort -k1,1 -k2,2n"

#names of input bedGraph files
inBG<-list.files(paste0(path,"/bed"),pattern="_M.bedGraph",full.names=T)

#names of output sorted bedGraph files
sortBG<-paste0(path,"/bed/", sub(".bedGraph","_s.bedGraph",basename(inBG)))

sortCmd<-paste(sortCmd,inBG, ">", sortBG)
#sorting bedGraph files
lapply(sortCmd,system)

#output bigwig filenames
outBW<-paste0(path,"/bed/", sub("bedGraph","bw",basename(inBG)))

#location of bedGraphToBigWig program
bg2bw<-"~/mySoftware/ucscUtils/bedGraphToBigWig.dms"

bwCmd<-paste(bg2bw, sortBG, chromSizesFile, outBW)

lapply(bwCmd,system)

system(paste0("rm ",path,"/bed/*.bedGraph"))

# #cutting the genome in equally sized bins for plotting
# seqLengths<-as.vector(df.chrom.sizes$chrSize)
# names(seqLengths)<-df.chrom.sizes$chrNames
# bins<-tileGenome(seqLengths,tilewidth=100,cut.last.tile.in.chrom=T)
#
# v009<-mcolAsRleList(gcMeth,"F2_DE_gwV009_M")
# avr<-binnedAverage(bins,v009,"F2_DE_gwV009_M",na.rm=T)
#
# winWidth=10
# ir=IRanges(start=seq(from=1,to=15072434-winWidth,by=winWidth),width=winWidth)
# agg<-aggregate(v009[[1]],ir,FUN=mean,na.rm=T)
#
# vi<-Views(v009[[1]],start=seq(from=1,to=1500000,by=10000),width=15098)
#
#
# vi<-Views(v009, as(bins,"RangesList"))
# mean(vi,na.rm=T)
# # few bugs in:
# # genome<-readDNAStringSet(genomeFile)
# # binViews<-Views(Celegans,bins)
# # gcCount<-dinucleotideFrequency(binViews, "GC")




##################################
### produce smoothed tracks   ###
#################################

removeAllNAs<-function(grObj,dataCols=c(1)) {
  # function to remove lines in a GRanges object where the given data columns all
  # have NA values. Input should be a GRanges object, and a vector with the names of the metadata
  # columns that you want to check for NAs (default is first mcols column).
  numCols<-length(dataCols)
  data2check<-mcols(grObj)[,dataCols]
  i<-which(rowSums(is.na(data2check))==numCols)
  return(grObj[-i])
}


length(meth_all)
#3131525 5918563
meth_all_s<-removeAllNAs(meth_all,sampleNames)
length(meth_all_s)
#3010924 5692560 5744817


smootheGRdata<-function(grObj,sampleNames,winSize=10) {
  #function to smoothe mcol data in a GRanges object with a window size of WinSize.
  #returns GRanges object with smooothed values (might be smaller as all NA only rows were removed)
  # first remove all rows that have only NAs in all the sample columns of interest
  grObj_noNA<-removeAllNAs(grObj,sampleNames)
  #cycle through the samples and smooothe the rle and convert back to GRanges
  for (s in sampleNames){
    dataCol<-as.vector(mcols(grObj_noNA)[,s])
    smData<-rollapply(dataCol,width=winSize,FUN=mean,na.rm=T)
    #pad smData ends with first and last values to get same length vector as input
    mcols(grObj_noNA)[,s]<-padEnds(smData,winSize)
  }
  return(grObj_noNA)
}


padEnds<-function(vec,winSize){
  #for roll apply you always loose winSize-1 values. Here i assume the data was centered
  #and i pad with end values
  padSize<-winSize-1
  leftPad<-padSize%/%2
  rightPad<-padSize-leftPad
  paddedVec<-c(rep(vec[1],leftPad),vec,rep(vec[length(vec)],rightPad))
  return(paddedVec)
}
smGRmethCG<-smootheGRdata(cgMeth,sampleNames,winSize=10)
smGRmethGC<-smootheGRdata(gcMeth,sampleNames,winSize=10)
smGRmethAllC<-smootheGRdata(meth_all,sampleNames,winSize=10)

#now make bigwig files from them
# create chrom sizes file needed by bedtobigwig command
#genome<-readDNAStringSet(genomeFile)
#df.chrom.sizes<- data.frame(chrNames=names(genome),chrSize=width(genome))
#write.table(df.chrom.sizes, "./bed/Ce11.chrom.sizes",quote=F,sep="\t",row.names=F,col.names=F)

chromSizesFile<-"~/Documents/MeisterLab/GenomeVer/chromSizes/ws235.chrom.sizes"
#unix command to sort bedGraph files (required for bedtobigwig command)
sortCmd<-"LC_COLATE=C sort -k1,1 -k2,2n"


#make input bedGraph files with smoothed data
gr2bedGraph(smGRmethCG,sampleNames,paste0(path,"/bed/CG_"),oneMinusScore=T)
gr2bedGraph(smGRmethGC,sampleNames,paste0(path,"/bed/GC_"),oneMinusScore=T)
gr2bedGraph(smGRmethAllC,sampleNames,paste0(path,"/bed/allCs_"),oneMinusScore=T)

inBG<-list.files(paste0(path,"/bed"),pattern="_M.bedGraph",full.names=T)

#names of output sorted bedGraph files
sortBG<-paste0(path,"/bed/", sub(".bedGraph","_s.bedGraph",basename(inBG)))

sortRun<-paste(sortCmd,inBG, ">", sortBG)
#sorting bedGraph files
lapply(sortRun,system)

#output bigwig filenames
outBW<-paste0(path,"/bed/sm_", sub("bedGraph","bw",basename(inBG)))

#location of bedGraphToBigWig program
bg2bw<-"~/mySoftware/ucscUtils/bedGraphToBigWig.dms"

bwCmd<-paste(bg2bw, sortBG, chromSizesFile, outBW)

lapply(bwCmd,system)

system(paste0("rm ",path,"/bed/*.bedGraph"))


# ##this function is not really what i want because i want a sliding window.
# binRanges<-function(grObj,winSize=10) {
#   # function to merge winSize number of ranges together and creating a new GRanges object that
#   # has the start and end of those 10 ranges. Need to deal with +-strand. Need to deal with
#   # no binning across chromosomes. Return GR object that spans winSize ranges from original grObj
#   allChr<-GRanges()
#   seqinfo(allChr)<-seqinfo(grObj)
#   for (chr in seqlevels(grObj)) {
#     oneChr<-grObj[seqnames(grObj)==chr]
#     strand(oneChr)<-"*"
#     startBins<-c(1,seq(winSize,length(oneChr),by=(winSize)))
#     endBins<-startBins[2:length(startBins)]
#     chrGR<-GRanges(seqnames=c(chr),
#                    ranges= IRanges(start=start(oneChr[startBins]),
#                                 end=c(start(oneChr[endBins])-1,end(oneChr[length(oneChr)]))),
#                    strand="*")
#     allChr<-append(allChr,chrGR)
#   }
#   return(allChr)
# }

#cg_bins<-binRanges(cgMeth_s)
#cgMeth_s$F2_DE_gwV008_M<-viewMeans(Views(rle,as(cg_bins, "RangesList"))) # this command doesn't work
