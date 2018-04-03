# .libPaths("/work2/gschub/altuna/RlibsDev")
# .libPaths("/work2/gschub/arnaud/Rpackages")
# .libPaths('/work2/gschub/darko/myRlibs/')
# library(AmpliconBiSeq)
require(data.table)
require(QuasR)
# arguments:
# proj: qProject object
# range: GRanges object with ONE!!!! range
# samp: sample.name
#
getCMethMatrix<-function(proj,range,samp){
  # function returns methylation matrix for an amplicon (or GR)
  # rows are individual reads, columns are methylation positions
  # input is a QuasR project, and GRange and sample names for the qProject
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
	return(CpGm)
}

#.libPaths("/work2/gschub/altuna/RlibsDev")
#.libPaths("/work2/gschub/arnaud/Rpackages")

#  require(AmpliconBiSeq)
require(QuasR)
require(data.table)
library(grid)

fuseReadMat <- function(resp,resm){ #assumes that reads are no on different strands
			uReads=unique(c(rownames(resp),rownames(resm)))
			uCs=unique(c(colnames(resp),colnames(resm)))
			matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
			colnames(matCGo)=uCs
			rownames(matCGo)=uReads
			matCGo[rownames(resp),colnames(resp)]=resp #fill up positive reads
			#add negative reads
			NAi=apply(resm,1,function(x){sum(is.na(x))==length(x)})
			matCGo[rownames(resm[!NAi,]),colnames(resm[!NAi,])]=resm[!NAi,]
			matCGo
		}
# 	GC_CG_mat=CmatS2CGs
mergeGC_CGmat=function(GC_CG_mat){
			CGmat=GC_CG_mat$matCG
			GCmat=GC_CG_mat$matGC
			uReads=unique(c(rownames(CGmat),rownames(GCmat)))
			uCs=sort(unique(c(colnames(CGmat),colnames(GCmat))))
			matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
			colnames(matCGo)=uCs
			rownames(matCGo)=uReads
			matCGo[rownames(CGmat),colnames(CGmat)]=CGmat#get CGs in
			matCGo[rownames(GCmat),colnames(GCmat)]=GCmat#get GCs in
	 		matCGo
}



getGCMatrix<-function(matC=matC.list[[1]],chr="chr2L",genome=Celegans,conv.rate=80,destrand=TRUE){
	Cpos=as.numeric(colnames(matC))
	rGR=GRanges(rep(chr,length(Cpos)),IRanges(Cpos,Cpos))
	rSeq=getSeq(genome,resize(rGR,3,fix='center'))
	GCpos=vcountPattern('GC',rSeq)==1
	CGpos=vcountPattern('CG',rSeq)==1
# 	GCGpos=vcountPattern('GCG',rSeq)
	GCi= GCpos #& !GCGpos
	CGi= CGpos #& !GCGpos
	Ci= Cpos & !GCpos & !CGpos #& !GCGpos

	########
	#Filter bases
	########

	# filter based on conversion
	ConvRate=100-(rowMeans( matC[,Ci,drop=FALSE],na.rm=T)*100)
	Convi=ConvRate>conv.rate
	matGC=matC[Convi,GCi,drop=FALSE] # get matrices
	matCG=matC[Convi,CGi,drop=FALSE]
	########
	#Destrand
	########
	# bases G or C
	# will help determine the strands
	bpsGC=as.character(getSeq(genome,rGR))[GCi]
	bpsCG=as.character(getSeq(genome,rGR))[CGi]

	if(destrand){
 		#Destrand GC
		# column names that are for plus and minus strand
		colMinus=(colnames(matGC))[bpsGC=="G"]
		colPlus= (colnames(matGC))[bpsGC=="C"]

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
		rownames(resp)=rownames(matGCP)
		resm=matrix(NA,ncol=length(unique(c(colnames(matGCM),colnames(matGCP)))),
		nrow=nrow(matGCM))
		rownames(resm)=rownames(matGCM)
		# get col names as base numbers
		cols=unique(c(colnames(matGCM),colnames(matGCP)))
		cols=cols[order(as.numeric(cols))]
		colnames(resp)=cols
		colnames(resm)=cols

		# populate temporary matrices
		resp[,colnames(resp) %in% colnames(matGCP)]=matGCP
		resm[,colnames(resm) %in% colnames(matGCM)]=matGCM
		rownames(resp)=rownames(matGCP)
		rownames(resm)=rownames(matGCM)
		#add the read names


		matGC=fuseReadMat(resp,resm)#cbind(resp,resm)

		matGC=matGC[,order(as.numeric(colnames(matGC))),drop=FALSE] # reorder columns

		#Destrand CG
		# column names that are for plus and minus strand
		colMinus=(colnames(matCG))[bpsCG=="G"]
		colPlus= (colnames(matCG))[bpsCG=="C"]

		# get the minus strand columns
		# remove empty rows
		# order by column names
		# add -1 to column names
		matCGM=matCG[, colnames(matCG) %in% colMinus,drop=FALSE]
		matCGM=matCGM[rowSums(!is.na(matCG))>0,,drop=FALSE]
		matCGM=matCGM[,order(as.numeric(colnames(matCGM))),drop=FALSE]
		colnames(matCGM)=as.numeric(colnames(matCGM))-1

		# get the plus strand columns
		# remove empty rows
		# order by column names
		matCGP=matCG[ ,colnames(matCG) %in% colPlus,drop=FALSE]
		matCGP=matCGP[rowSums(!is.na(matCGP))>0,,drop=FALSE]
		matCGP=matCGP[,order(as.numeric(colnames(matCGP))),drop=FALSE] # reorder columns

		# temp matrices for plus and minus strand reads
		resp=matrix(NA,ncol=length(unique(c(colnames(matCGM),colnames(matCGP)))),nrow=nrow(matCGP))
		resm=matrix(NA,ncol=length(unique(c(colnames(matCGM),colnames(matCGP)))),nrow=nrow(matCGM))
		rownames(resp)=rownames(matCGP)
		rownames(resm)=rownames(matCGM)
		# get col names as base numbers
		cols=unique(c(colnames(matCGM),colnames(matCGP)))
		cols=cols[order(as.numeric(cols))]
		colnames(resp)=cols
		colnames(resm)=cols

		# populate temporary matrices
		resp[,colnames(resp) %in% colnames(matCGP)]=matCGP
		resm[,colnames(resm) %in% colnames(matCGM)]=matCGM
		rownames(resp)=rownames(matCGP)
		rownames(resm)=rownames(matCGM)
		#merge the two matrices

		matCG=fuseReadMat(resp,resm)

		matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE] # reorder columns

		}else{

			matGC=matGC[,order(as.numeric(colnames(matGC))),drop=FALSE] # reorder columns
			matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE]
		}
	########
	#Remove the emty row and collumns
	########


	matGC=matGC[,colSums(!is.na(matGC))>0,drop=FALSE] # remove no coverage columns
	matCG=matCG[,colSums(!is.na(matCG))>0,drop=FALSE]
	matGC=matGC[rowSums(!is.na(matGC))>0,,drop=FALSE] # remove no coverage rows
	matCG=matCG[rowSums(!is.na(matCG))>0,,drop=FALSE]
	list(matGC=matGC,matCG=matCG)
}

miObj2gr<-function(miObj) {
#function to convert Mindex object (obtained from applying GR on DNAstringset) to genomic ranges
	allGR<-GRanges()
	seqlevels(allGR)<-names(miObj)
	for (n in names(miObj)) {
		grObj<-GRanges(seqnames=Rle(c(n),length(miObj[[n]])),
			ranges=miObj[[n]], strand="*")
		allGR<-append(allGR,grObj)
	}
	return(allGR)
}


call_context_methylation=function(meth_gr,cO,genome=Celegans){
  #call_context_methylation returns list of two GRanges objects "CG" and "GC"
  # in each of these , the V1 column contains fraction methylation
  # rather than counts and 'type' column with C context
			if (class(genome)=="BSgenome"){
				fastaFlag=FALSE
			} else if (is.character(genome)){
				genome<-readDNAStringSet(genome)
				fastaFlag=TRUE
			} else {
				print("genome must be BSgenome of path to fasta file")
			}		
			genome_CGs <- vmatchPattern("CG",genome)
			genome_GCs <- vmatchPattern("GC",genome)

			if (fastaFlag==FALSE) {
				#if the genome is a BSgenome object then you have to get rid
				#of duplicate matches on positive and negative strand
				genome_CGs <- genome_CGs[strand(genome_CGs)=="+",]
				strand(genome_CGs) <- "*"

				genome_GCs <- genome_GCs[strand(genome_GCs)=="+",]
				strand(genome_GCs) <- "*"
			} else {
				# if the genome if fasta file, we need to convert Mindex object
				#  to GRanges object
				genome_CGs<-miObj2gr(genome_CGs)
				genome_GCs<-miObj2gr(genome_GCs)
			}

			sel_CGs <- meth_gr %over% genome_CGs
			sel_GCs <- meth_gr %over% genome_GCs

			meth_CGs_gr <- meth_gr[sel_CGs]
			meth_GCs_gr <- meth_gr[sel_GCs]

		#### check arrays have even number of rows and if not,
		#### add emptry row to avoid error when adding matrices
			if(length(meth_CGs_gr)%%2==1) {
			   lastRange<-meth_CGs_gr[length(meth_CGs_gr)]
			   emptyRange<-lastRange
			   mcols(emptyRange)<-matrix(rep(0,length(mcols(emptyRange))),nrow=1)
			   colnames(mcols(emptyRange))<-colnames(mcols(lastRange))
			   meth_CGs_gr<-c(meth_CGs_gr,emptyRange)
			}
			if(length(meth_GCs_gr)%%2==1) {
			   lastRange<-meth_GCs_gr[length(meth_GCs_gr)]
			   emptyRange<-lastRange
			   mcols(emptyRange)<-matrix(rep(0,length(mcols(emptyRange))),nrow=1)
			   colnames(mcols(emptyRange))<-colnames(mcols(lastRange))
			   meth_GCs_gr<-c(meth_GCs_gr,emptyRange)
			}

		##################
		# collapse strands
		##################
			meth_CGsCol_gr <- meth_CGs_gr[seq(1,length(meth_CGs_gr),by=2)]
			end(meth_CGsCol_gr) <- end(meth_CGsCol_gr)+1
			values(meth_CGsCol_gr) <- as.matrix(values(meth_CGs_gr[seq(1,length(meth_CGs_gr),by=2)]))+as.matrix(values(meth_CGs_gr[seq(2,length(meth_CGs_gr),by=2)]))

			meth_GCsCol_gr <- meth_GCs_gr[seq(2,length(meth_GCs_gr),by=2)]
			start(meth_GCsCol_gr) <- start(meth_GCsCol_gr)-1
			values(meth_GCsCol_gr) <- as.matrix(values(meth_GCs_gr[seq(1,length(meth_GCs_gr),by=2)]))+as.matrix(values(meth_GCs_gr[seq(2,length(meth_GCs_gr),by=2)]))

		#####################
		#filter for coverage
		####################
			Tcounts=grep('_T\\>',colnames(elementMetadata(meth_CGsCol_gr)))
			Mcounts=grep('_M\\>',colnames(elementMetadata(meth_CGsCol_gr)))

			######
			#CGs
			######

			CG.met.mat=as.matrix(elementMetadata(meth_CGsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])
			#filter for coverage
			CovFilter=as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])>cO
			for (i in sl(CG.met.mat[1,])){CG.met.mat[!CovFilter[,i],i]=NA}
				#bind the GRanges with the scores
			CG.met=meth_CGsCol_gr
			elementMetadata(CG.met)=CG.met.mat


			######
			#GCs
			######

			GC.met.mat=as.matrix(elementMetadata(meth_GCsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])
			#filter for coverage
			CovFilter=as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])>cO
			for (i in sl(GC.met.mat[1,])){GC.met.mat[!CovFilter[,i],i]=NA}
			#bind the GRanges with the scores

			GC.met=meth_GCsCol_gr
			elementMetadata(GC.met)=GC.met.mat

			###########################
			#create an unified object
			###########################
			oGCG=as.matrix(findOverlaps(resize(GC.met,1,fix='end'),resize(CG.met,1,fix='start'),type='equal'))
			GC.met$type='GCH'
			CG.met$type='CGH'
			GC.met$type[oGCG[,1]]='GCG'
			CG.met$type[oGCG[,2]]='GCG'

		# 	getSeq(Celegans,resize(head(GC.met[oGCG[,1]]),fix='center'))
		# 	x=GC.met[oGCG[,1]]$OSC_SS_R2_M
		# 	y=CG.met[oGCG[,2]]$OSC_SS_R2_M

		#
		# 	GCs=GC.met[!(sl( GC.met) %in% oGCG[,1])  ]
		# 	GCs$type='GCH'
		# 	CGs=CG.met[!(sl( CG.met) %in% oGCG[,2])  ]
		# 	CGs$type='HCG'
		# 	GCGs=GC.met[( oGCG[,1])  ]
		# 	GCGs$type='GCG'
			umet=list(CG.met,GC.met)
			names(umet)=c('CG','GC')
			return(umet)
		}


# horizontalHeat function derived from AmpliconBiSeq package
# horizontal heatmaps


#' make horizontal heatmap for LD block type of analysis
#'
#' @param sim.mat n-by-n symetric matrix
#' @param legend.width in lines defaults to 1
#' @param legend logical, indicating if a color legend be drawn def(legend=TRUE)
#' @param legend.text legend text
#' @param col.select vector of colors
#' @param ... some arguments that works for plot() function from base
#'
#' @examples
#' horizontalHeat(matrix(rnorm(10000),ncol=100,nrow=100),legend.width=0.5)
#' horizontalHeat(matrix(rnorm(10000),ncol=100,nrow=100),legend=FALSE)
#'
#'
#' par(mfrow=c(2,1))
#' par(mar=c(2.1,4.1,2.1,2.1))
#' plot(1:20,rnorm(20),type="b")
# horizontalHeat(matrix(rnorm(900),ncol=30,nrow=30),legend=TRUE)
horizontalHeat<-function(sim.mat,legend.width=1, legend=TRUE,legend.text="similarity",
                        col.select=colorRampPalette(c("black", "yellow","purple"))(50),...)
   {
   require(gridBase)
   # trick to scale colors from 0->1 add a 0 to the diagonal
   #sim.mat[1,1]=0
   #sim.mat[2,2]=1
   col.mat=convertToColors(sim.mat,col.select) # convert to colors
   #get range
   my.mat=sim.mat
   diag(my.mat)=NA
   rng=range(my.mat,na.rm=TRUE)
   # figure out locations of cells of the 45 degree
   # rotated heatmap
   # we will create the heatmap with hexagonal cells
   len=ncol(sim.mat)
   run.len=len-1
   Y<-X<-mX<-mY<-c()
   for(i in 1:(len-1) ){
      X=c(X,seq(i+0.5,by=0.5,len=run.len))
      #Y=c(Y,-(1:run.len) )
      Y=c(Y,-(seq(1,by=0.5,length.out=run.len)))
      mX=c(mX,rep(i,run.len))
      mY=c(mY,seq(to=len,by=1,length.out=run.len))
      run.len=run.len-1
   }
   # PART 0: PLOT HEATMAP
   #_____________________________________________________________________________
   plot(X, Y, type="n",yaxt="n",bty='n',xaxt="n",ylab=NA)
   vps <- baseViewports()
   pushViewport(vps$inner, vps$figure, vps$plot)
   #heatvp<-viewport(x = unit((legend.width+spacing), "npc"), y = unit(0.5, "npc"),
   # width = unit(1-(legend.width+spacing+flank), "npc"),
   # height = unit(1, "npc"),
   # just="left",
   # xscale = c(0.5,ncol(sim.mat)+.5 ),
   # yscale = c(-(ncol(sim.mat)+1)*0.5,0.5) )
   #pushViewport(heatvp)
   #grid.rect()
   #grid.yaxis()
   # get colors from color matrix
   my.cols=apply(cbind(mX,mY),1, function(x,col.mat) col.mat[x[2],x[1]],col.mat=col.mat )
   #plot the cells
   # +0.3 on Y controls the relative height of the heatmap within the plot
   sqpolygon(X,Y+0.3, hexC=sqcoords(dx = 0.5,dy=0.5, sep=NULL),
             border = "white", fill=my.cols)
   popViewport(3)
   #upViewport()
   if(legend){
      if(par()$mar[2]<4.1){
         warning("left margin of the plot (set by mar in par()) should not be less then 4.1 lines")
      }
      pushViewport(vps$figure) # get the plot viewport from base graphics
      #showViewport(current.viewport());current.vpTree()
      #grid.text(c("one"),
      # x=unit(1, "native"), y=unit(2, "native"),
      # just="right", rot=60)
      # make view port for the legend
      legendVp <- viewport(width=unit(legend.width, "lines"), height=unit(0.4, "npc"),
                           x = unit(3, "lines"), y = unit(0.5, "npc"),just="left")
      pushViewport(legendVp) # push the legend VP
      grid.rect()
      grid.yaxis(gp=gpar(cex=0.7) )
      #my.cols=colorRampPalette(c("black", "yellow","red","purple"))(50)
      grid.raster(rev(convertToColors(seq(rng[1],rng[2],length.out=20),col.select)), interpolate = TRUE,
                  height=unit(1,"npc"),width=1)
      #grid.raster(rev(convertToColors(seq(0,1,by=0.1),my.cols )), interpolate = TRUE)
      upViewport()
      grid.text(legend.text,x = unit(0.01,"npc"), rot=90,just="left")
      upViewport()
   }
}





convertToColors <- function(mat,col.select=colorRampPalette(c("blue", "yellow","red"))(50) ) {
   # Produce 'normalized' version of matrix, with values ranging from 0 to 1
   if(is.matrix(mat)){
      diag(mat)=NA
   }
   rng <- range(mat, na.rm = TRUE)
   m <- (mat - rng[1])/diff(rng)
   # Convert to a matrix of sRGB color strings
   m2 <- m; class(m2) <- "character"
   m2[!is.na(m2)] <- rgb(colorRamp(col.select)(m[!is.na(m)]), max = 255)
   m2[is.na(m2)] <- "transparent"
   return(m2)
}


# calculates coordinates for squares to be used in sqpolygon
#
sqcoords<-function (dx=0.5, dy = NULL, n = 1, sep = NULL)
{
   stopifnot(length(dx) == 1)
   if (is.null(dy))
      dy <- dx
   if (is.null(sep))
      list(x = rep.int(c(dx, 0,-dx, 0), n), y = rep.int(c(0, -dy, 0, dy), n), no.sep = TRUE)
   else list(x = rep.int(c(dx, 0,-dx, 0, sep), n),
             y = rep.int(c(0 ,-dy, 0, dy, sep),
                         n), no.sep = FALSE)
}

# modified version of hexypolygon from hexbin, where it plots squares
# instead of polygons
sqpolygon<-function(x, y, hexC = sqcoords(dx, dy, n = 1), dx, dy = NULL,
fill = 1, border = 0, hUnit = "native", ...)
{
   n <- length(x)
   stopifnot(length(y) == n)
   stopifnot(is.list(hexC) && is.numeric(hexC$x) && is.numeric(hexC$y))
   if (hexC$no.sep) {
      n6 <- rep.int(4:4, n)
      if (!is.null(hUnit)) {
         grid.polygon(x = unit(rep.int(hexC$x, n) + rep.int(x, n6), hUnit),
                      y = unit(rep.int(hexC$y, n) + rep.int(y, n6), hUnit), id.lengths = n6,
                      gp = gpar(col = border, fill = fill))
      }
      else {
         grid.polygon(x = rep.int(hexC$x, n) + rep.int(x,
                                                       n6), y = rep.int(hexC$y, n) + rep.int(y, n6),
                      id.lengths = n6, gp = gpar(col = border, fill = fill))
      }
   }
   else {
      n7 <- rep.int(7:7, n)
      polygon(x = rep.int(hexC$x, n) + rep.int(x, n7), y = rep.int(hexC$y,
                                                                   n) + rep.int(y, n7), ...)
   }
}



#' make horizontal heatmap for LD block type of analysis
#'
#' @param sim.mat n-by-n symetric matrix
#' @param legend.width in lines defaults to 1
#' @param legend logical, indicating if a color legend be drawn def(legend=TRUE)
#' @param legend.text legend text
#' @param col.select vector of colors
#' @param ... some arguments that works for plot() function from base
#'
#' @examples
#' horizontalHeat(matrix(rnorm(10000),ncol=100,nrow=100),legend.width=0.5)
#' horizontalHeat(matrix(rnorm(10000),ncol=100,nrow=100),legend=FALSE)
#'
#'

#' par(mfrow=c(2,1))
#' par(mar=c(2.1,4.1,2.1,2.1))
#' plot(1:20,rnorm(20),type="b")
#' horizontalHeat(matrix(rnorm(900),ncol=30,nrow=30),legend=TRUE)
