#library('BSgenome')
#library('BSgenome.Celegans.UCSC.ce11')
#Çlibrary("BSgenome.Ecoli.NCBI.20080805")
library(Rsamtools)
library(RColorBrewer)

sl=function(x){seq(length(x))}

filterCov=function(frag, cutOff){
  return(frag[elementMetadata(frag)$biscount/elementMetadata(frag)$nb.CG>cutOff] )
                     }
filternbCG=function(frag, cutOff){
  return(frag[elementMetadata(frag)$Ntot>cutOff] )
                     }

panel.jet <- function(...) {
                smoothScatter(..., nrpoints=0, add=TRUE, colramp=jet.colors) }


panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
     usr <- par("usr"); on.exit(par(usr))
     par(usr = c(0, 1, 0, 1))
     r <- (cor(x, y))
     txt <- format(c(r, 0.123456789), digits=digits)[1]
     txt <- paste(prefix, txt, sep="")
     if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

     test <- cor.test(x,y)
     # borrowed from printCoefmat
     Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))

     text(0.5, 0.5, txt, cex=0.6*cex)
    # text(.8, .8, Signif, cex=cex, col=2) 
}

jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



colgradPlot <- function(x, y, cor=FALSE, drawLegend=FALSE, colramp=colorRampPalette(c("darkblue", "lightblue","yellow"), space="Lab"), transformation = function(x) x^0.25, nbin=128, bandwidth, add=FALSE, ...)
	{
  
  colorVector <- MydensCols(cbind(x,y), colramp=colramp, transformation=transformation, 
nbin=nbin, bandwidth=bandwidth)

  if(add==FALSE){
    plot(x,y, col=colorVector, ...)
  }
  else{
    points(x,y, col=colorVector, ...)
  }

  if(cor==TRUE){
    legend(x="topleft", bty="n", legend=sprintf("R=%.3f", cor(x, y)))
    drawLegend <- match.arg(drawLegend)
  }
  
  xy.range<- par("usr")
  
  if (drawLegend) {
    color.legend(xl=xy.range[2]/2
                 ,xr=xy.range[2]/1.5
                 ,yb=xy.range[3]/2,
                 yt= xy.range[3]/1.3,
                 legend=c("low","high"),
                 rect.col=colramp(256),
                 gradient="y", #vertical, x=horizontal gradient
                 align="rb")
  }
}

MydensCols <- function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(blues9[-(1:3)]), transformation){
    xy <- xy.coords(x, y)
    select <- is.finite(xy$x) & is.finite(xy$y)
    x <- cbind(xy$x, xy$y)[select, ]
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
    xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
    dens <- transformation(map$fhat)[cbind(xbin, ybin)]
    dens[is.na(dens)] <- 0
    colpal <- cut(dens, length(dens), labels = FALSE)
    cols <- rep(NA_character_, length(select))
    cols[select] <- colramp(length(dens))[colpal]
    cols
}




 sigModelCG=function(CGdensity){
   100 / (1 + exp(-co[1] * (CGdensity - co[2])))
  }
 #co= readRDS('/work2/gschub/arnaud/DNAlib/meta_analysis/model/sigModelCoef.rds')


#useful palettes
	colSet=colorRampPalette(brewer.pal(9,"Set1"))(9)
#

#concerts all C to T in a non CG context
bisConv=function(in.seq){
  C.pos=vmatchPattern('C', in.seq)
  CG.pos=vmatchPattern('CG', in.seq)
  do.call(c,
    lapply(seq(length(in.seq)),function(i){
      Cind=start(C.pos[[i]])
      CGind=start(CG.pos[[i]])
      convind=Cind[!Cind %in% CGind]
      DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
      })
    )
  }
#concerts all Cto T in a non CG or GCcontext
  
bisConvGC=function(in.seq){
  C.pos=vmatchPattern('C', in.seq)
  CG.pos=vmatchPattern('CG', in.seq)
  GC.pos=vmatchPattern('GC', in.seq)
  
  do.call(c,
    lapply(seq(length(in.seq)),function(i){
      Cind=start(C.pos[[i]])
      CGind=start(CG.pos[[i]])
      GCind=start(GC.pos[[i]])
      convind1=!Cind %in% CGind  
      convind2=! Cind %in% (GCind+1)
      convind= Cind[convind1 & convind2]
      DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
      })
    )
  }

 #########################
grToNames =function(gr){
paste(as.character(seqnames(gr)),':',as.character(start(gr)),'-',as.character(end(gr)), sep='')
}



bamTOgr=function(IN){
  INbam=scanBam(IN, param=ScanBamParam(what=c("rname","pos", "strand", "qwidth")))
  matched=!is.na(INbam[[1]]$rname)
  fw=INbam[[1]]$strand[matched]=='+'
  rv=INbam[[1]]$strand[matched]=='-'
  grf=GRanges(seqnames=INbam[[1]]$rname[matched][fw], ranges=IRanges(INbam[[1]]$pos[matched][fw],(INbam[[1]]$pos[matched][fw]+INbam[[1]]$qwidth[matched][fw])), strand=INbam[[1]]$strand[matched][fw])
  grr=GRanges(seqnames=INbam[[1]]$rname[matched][rv], ranges=IRanges(INbam[[1]]$pos[matched][rv]-INbam[[1]]$qwidth[matched][rv],(INbam[[1]]$pos[matched][rv])), strand=INbam[[1]]$strand[matched][rv])
  INgr=c(grf, grr)
  INgr
 # names(INgr)=INbam[[1]]$qname
  }

bamTOgrChIP=function(IN,fs){
  INbam=scanBam(IN, param=ScanBamParam(what=c("rname","pos", "strand", "qwidth")))
  matched=!is.na(INbam[[1]]$rname)
  fw=INbam[[1]]$strand[matched]=='+'
  rv=INbam[[1]]$strand[matched]=='-'
  grf=GRanges(seqnames=INbam[[1]]$rname[matched][fw], ranges=IRanges(INbam[[1]]$pos[matched][fw],(INbam[[1]]$pos[matched][fw]+INbam[[1]]$qwidth[matched][fw])), strand=INbam[[1]]$strand[matched][fw])
  grr=GRanges(seqnames=INbam[[1]]$rname[matched][rv], ranges=IRanges(INbam[[1]]$pos[matched][rv]-INbam[[1]]$qwidth[matched][rv],(INbam[[1]]$pos[matched][rv])), strand=INbam[[1]]$strand[matched][rv])
  INgr=c(grf, grr)
  INgr
 # names(INgr)=INbam[[1]]$qname
  }
bamTOgrMNase=function(IN){
	INbam=scanBam(IN, param=ScanBamParam(what=c("rname","pos", "strand", "qwidth")),isPaired=TRUE)
	fwi=seq(1,length(INbam[[1]]$rname),2)
	matched=!is.na(INbam[[1]]$rname[fwi]) & !is.na(INbam[[1]]$rname[fwi+1])
  	fwo=INbam[[1]]$strand[fwi][matched]=='+'
  	rvo=INbam[[1]]$strand[fwi][matched]=='-'
# 	grf=GRanges(seqnames=INbam[[1]]$rname[fwi][matched][fwo], ranges=IRanges(INbam[[1]]$pos[fwi][matched][fwo],(INbam[[1]]$pos[fwi][matched][fwo]+INbam[[1]]$qwidth[fwi][matched][fwo])), strand=INbam[[1]]$strand[fwi][matched][fwo])
  grf=GRanges(seqnames=INbam[[1]]$rname[fwi][matched][fwo], 
  	ranges=IRanges(INbam[[1]]$pos[fwi][matched][fwo],
  	(INbam[[1]]$pos[fwi+1][matched][fwo]+INbam[[1]]$qwidth[fwi+1][matched][fwo])), 
  	strand=INbam[[1]]$strand[fwi][matched][fwo])
 
   grrv=GRanges(seqnames=INbam[[1]]$rname[fwi][matched][rvo], 
   	ranges=IRanges(INbam[[1]]$pos[fwi+1][matched][rvo],
 	(INbam[[1]]$pos[fwi][matched][rvo]+INbam[[1]]$qwidth[fwi][matched][rvo])), 
  	strand=INbam[[1]]$strand[fwi][matched][rvo])
  	gr=c( grf,grrv)
  	gr
  	}

bamTOgrDHS=function(fh){
	gr=bamTOgr(fh)
	cutgr=GRanges(seqnames(gr),IRanges(
	ifelse(strand(gr)=='+',start(gr),end(gr)),ifelse(strand(gr)=='+',start(gr),end(gr))
	))
	cutgr
	}


string.split=function(string,sep,pos){unlist(lapply(string,function(x){lapply(strsplit(x,sep),'[',pos)}))}


getAverageAmplicons=function(amplicon.list,coverage,minnbCG){
	samples=getSampleNames(amplicon.list)
	amplicons=getAmpliconNames(amplicon.list)
	metMat=matrix(NA,nrow=length(amplicons),ncol=length(samples))
	colnames(metMat)=samples
	rownames(metMat)=amplicons
	primersgr=getAmpliconRanges(a.list)
	Amet=for(sp in sl(samples)){
		meanMet=unlist(lapply(sl(primersgr),function(i){
			CGcov=getCoverage(a.list, samples[sp],amplicons[i])
			id1=CGcov >= coverage
			id2=length(CGcov[id1]) >= minnbCG
			if(id2){
			mean(getAvMeth(a.list, samples[sp],amplicons[i])[id1])
			}else{NA}
	}))
	metMat[,sp]=meanMet
	}
	metMat
}


tileCounts=function(regions,reads,wind,st){
	mclapply(sl(regions), function(i){
		bins=seq(start(regions[i]),end(regions[i]),st)
		gr=GRanges(rep(as.character(seqnames(regions[i])),length(bins)),IRanges(bins-wind,bins+wind))
		counts=countOverlaps(gr,reads)
		names(counts)=bins
		counts
	},mc.cores=10)}

 e=function(fgr,bgr,pc,Tfg,Tbg){

   nfg=((fgr+pc)/(Tfg/1e6))
   nbg=((bgr+pc)/(Tbg/1e6))
   nfg/nbg
    }
#require(AmpliconBiSeq)
require(QuasR)
require(data.table)
#' get C methylation matrix for a given amplicon
# arguments:
# proj: qProject object
# range: GRanges object with ONE!!!! range
# samp: sample.name
#
getCMethMatrix<-function(proj,range,samp){
   Cs=qMeth(proj, query=range,mode="allC",reportLevel="alignment")
   # use data.table to get a 1,0 matrix of methylation profiles
   all.cids=unique(Cs[[samp]]$Cid) # get all possible C locations
   # make the data.table object
   dt=data.table(meth=Cs[[samp]]$meth ,aid=Cs[[samp]]$aid ,cid=Cs[[samp]]$Cid)
   # this function converts cids to columns
   myfun2<-function(x,all.cids){
      vec=rep(-1,length(all.cids))
      names(vec)=as.character(all.cids)
      b=as.list((vec))
      b[ as.character(x$cid)]=as.double(x$meth)
      return(b)
   }
   dtm=dt[,myfun2(.SD,all.cids), by=aid]
   ronames=dtm$aid
   dtm[,aid:=NULL] # remove unwanted row
   CpGm=as.matrix(dtm)
   CpGm[CpGm == -1]=NA # put NAs
   rownames(CpGm)=ronames
   # remove columns if they are outside the window of interest
   CpGm=CpGm[,start(range) <= as.numeric(colnames(CpGm)) & end(range) >= as.numeric(colnames(CpGm)) ]
   return(CpGm)
}
 
#' Get GC matrix from the C matrix
#'
#' function filters the matrix based on conversion rate and gets the matrix for GC
#'
#' @param matC is the C matrix with colnames corresponding to bp positions
#' @param chr chromosome names
#' @param genome
#' @param conv.rate conversion rate, default 80, a value between 0 and 100.
#' @param destrand, if TRUE the reads are destranded so that minus strand Cs
#' are converted to plus strand in their coordinates
#'
#' @return list(matGC=matC[Convi,GCi],matCG=matC[Convi,CGi])
getGCMatrix<-function(matC=matC.list[[1]],chr="chrI",genome=Celegans,conv.rate=80
,destrand=TRUE){
   Cpos=as.numeric(colnames(matC))
   rGR=GRanges(rep(chr,length(Cpos)),IRanges(Cpos,Cpos))
   rSeq=getSeq(genome,resize(rGR,3,fix='center'))
   GCpos=vcountPattern('GC',rSeq)
   CGpos=vcountPattern('CG',rSeq)
   GCGpos=vcountPattern('GCG',rSeq)
   GCi= GCpos & !GCGpos
   CGi= CGpos & !GCGpos
   Ci= Cpos & !GCpos & !CGpos
   #CmatC=cbind(-1*(Cmatfu[,GCi]) , (Cmatfu[,CGi]))
   #pos=c(Cpos[GCi],Cpos[CGi])
   # filter based on conversion
   ConvRate=100-(rowMeans( matC[,Ci,drop=FALSE],na.rm=T)*100)
   Convi=ConvRate>conv.rate
   matGC=matC[Convi,GCi,drop=FALSE] # get matrices
   matCG=matC[Convi,CGi,drop=FALSE]
   # bases G or C
   # will help determine the strands
   bps=as.character(getSeq(genome,rGR))[GCi]
   if(destrand){
      
      # column names that are for plus and minus strand
      colMinus=(colnames(matGC))[bps=="G"]
      colPlus= (colnames(matGC))[bps=="C"]
      # get the minus strand columns
      # remove empty rows
      # order by column names
      # add +1 to column names
      matGCM=matGC[, colnames(matGC) %in% colMinus,drop=FALSE]
      matGCM=matGCM[rowSums(!is.na(matGC))>0,,drop=FALSE]
      matGCM=matGCM[,order(as.numeric(colnames(matGCM))),drop=FALSE]
      colnames(matGCM)=as.numeric(colnames(matGCM))+1
      # get the plus strand columns
      # remove empty rows
      # order by column names
      matGCP=matGC[ ,colnames(matGC) %in% colPlus,drop=FALSE]
      matGCP=matGCP[rowSums(!is.na(matGCP))>0,,drop=FALSE]
      matGCP=matGCP[,order(as.numeric(colnames(matGCP))),drop=FALSE] # reorder columns
      # temp matrices for plus and minus strand reads
      resp=matrix(NA,ncol=length(unique(c(colnames(matGCM),colnames(matGCP)))),
                  nrow=nrow(matGCP))
      resm=matrix(NA,ncol=length(unique(c(colnames(matGCM),colnames(matGCP)))),
                  nrow=nrow(matGCM))
      # get col names as base numbers
      cols=unique(c(colnames(matGCM),colnames(matGCP)))
      cols=cols[order(as.numeric(cols))]
      colnames(resp)=cols
      colnames(resm)=cols
      # populate temporary matrices
      resp[,colnames(resp) %in% colnames(matGCP)]=matGCP
      resm[,colnames(resm) %in% colnames(matGCM)]=matGCM
      matGC=rbind(resp,resm)
      matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE]
   }else{
      matGC=matGC[,order(as.numeric(colnames(matGC))),drop=FALSE] # reorder columns
      matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE]
   }
   matGC=matGC[,colSums(!is.na(matGC))>0,drop=FALSE] # remove no coverage columns
   matCG=matCG[,colSums(!is.na(matCG))>0,drop=FALSE]
   matGC=matGC[rowSums(!is.na(matGC))>0,,drop=FALSE] # remove no coverage rows
   matCG=matCG[rowSums(!is.na(matCG))>0,,drop=FALSE]
   list(matGC=matGC,matCG=matCG)
}

grToIgv<-function(gr,dataset,fileName){
	scores=100-round((elementMetadata(gr)[,dataset]*100))
	NAi=is.na(scores) #remove uncovered
	grf=gr[!NAi]
	df <- data.frame(Chromosome=seqnames(grf),
  	Start=start(grf)-1,
  	End=end(grf),
  	Feature=c(rep(dataset, length(grf))),
  	R1=scores[!NAi])
	
	write.table(df, file=paste(fileName,'.igv',sep=''), quote=F, sep="\t", row.names=F, col.names=T)

# 	system(paste('sort -k1,1 -k2,2n ',fileName,'.igv ', '> ',fileName,'S.igv',sep=''))
	}





jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
 
panel.jet <- function(...) {
  smoothScatter(..., nrpoints=0, add=TRUE, colramp=jet.colors) }


panel.hist <- function(x, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
     }

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y, use="pairwise.complete.obs"))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex.cor <- 2.5/strwidth(txt)
         text(0.5, 0.5, txt)
     }


