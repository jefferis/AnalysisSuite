# NeuronFunctions5.s
# #################################
# This file contains general functions for handling Neurons:
# It has been shortened since version 4 by the removal of ParseSWCTree
# ReadSWCFile and reroot which have been moved to a new file SWCFunction.s
# Plot fns need to be updated to handle MB data
# - Have done this more or less I think need to verify esp File Plots
# #################################
# GSXEJ 020629

#RELEASE
#BEGINCOPYRIGHT
###############
# R Source Code to accompany the manuscript
#
# "Comprehensive Maps of Drosophila Higher Olfactory Centers: 
# Spatially Segregated Fruit and Pheromone Representation"
# Cell (2007), doi:10.1016/j.cell.2007.01.040
# by Gregory S.X.E. Jefferis*, Christopher J. Potter*
# Alexander M. Chan, Elizabeth C. Marin
# Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo
#
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# See flybrain.stanford.edu for further details
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################
#ENDMAINCOPYRIGHT

# SegLengths<-function(MyNeuron,SegmentArray)
# seglength<-function(ThisSeg)
# ReadMyChullOutput<-function(FileName)
# SegOrders<-function(MyNeuron){
# plotneuron2d<-function(ANeuron,WithLine=T,NodesOnly=T,WithText=F)
# plotneuron3d<-function(ANeuron,UseCurPalette=F,WithContours=T,WithScale=T,
# rgbcolour<-function(colour){
# plotendpoint3d<-function(ANeuron,UseCurPalette=F,WithContours=F,WithScale=T,
# plotneurons2d<-function(NeuronRef,Ask=T,...){
# plotneurons3d<-function(NeuronRef,Ask=T,ToFile=F,...){
# plotendpoints3d<-function(NeuronRef,Ask=T,ToFile=T,Col=1,FileAppend=T,...){
# GetNeuronNum
# FirstnNeurons<-function(n=1){
# NeuronNameFromFileName<-function(FileName){
# ImageName<-function(NeuronName){
# GetCellType<-function(NeuronRef,mask=1:length(MyNeurons)){
# GetNeuronNumsofType<-function(CellTypesToMatch,mask=1:length(MyNeurons)){
# GetNeuron<-function(NeuronR	ef,mask=1:length(MyNeurons)){    
# GetNeuronName<-function(Nnum,mask=1:length(MyNeurons)){
# GetNeuronNum<-function(Nnames,mask=1:length(MyNeurons)){

# source("Greg Data:AnalysingTraceFiles:R:NeuronFunctions5.s")
# source(file.path(CodeDir,"NeuronFunctions5.s"))

if(!require(rgl) && !require(scatterplot3d)){
	stop("Please install either rgl or scatterplot3d for 3d plotting")
}   # for plotneuron3d()

is.neuron<-function(n,Strict=FALSE) {
	# If Strict is FALSE will also return TRUE
	# if n is a list which looks like a neuron
	inherits(n,"neuron") ||
		(!Strict && is.list(n) && !is.null(n$SegList))
}

is.neuronlist<-function(nl) {
	inherits(nl,"neuronlist") ||
		(is.list(nl) && length(nl)>1 && is.neuron(nl[[1]]))
}

#' Arithmetic for neuron coordinates
#'
#' If x is one number or 4-vector, multiply xyz and diameter by that
#' If x is a 3-vector, multiply xyz only
#' TODO Figure out how to document arithemtic functions in one go
#' @param n a neuron
#' @param x (a numeric vector to multiply neuron coords in neuron)
#' @return modified neuron
#' @export
#' @examples
#' n1<-MyNeurons[[1]]*2
#' n2<-MyNeurons[[1]]*c(2,2,2,2)
#' stopifnot(all.equal(n1,n2))
#' n3<-MyNeurons[[1]]*c(2,2,4)
"*.neuron" <- function(n,x) {
	# TODO look into S3 generics for this functionality
	
	nd=n$d[,c("X","Y","Z","W")]
	stopifnot(is.numeric(x))
	lx=length(x)
	if(lx==1) nd[,-4]=nd[,-4]*x
	else if(lx%in%c(3,4)) nd[,1:lx]=t(t(nd[,1:lx])*x)
	else stop("expects a numeric vector of length 1, 3 or 4")
	n$d[,colnames(nd)]=nd
	n
}

"+.neuron" <- function(n,x) {
	if(!is.numeric(x))
		stop("expects a numeric vector")
	nd=n$d[,c("X","Y","Z","W")]
	lx=length(x)
	if(lx==1) nd[,-4]=nd[,-4]+x
	else if(lx%in%c(3,4)) nd[,1:lx]=t(t(nd[,1:lx])+x)
	else stop("expects a numeric vector of length 1, 3 or 4")
	n$d[,colnames(nd)]=nd
	n
}

"-.neuron" <- function(n,x) n+(-x)
"/.neuron" <- function(n,x) n*(1/x)

#' Divide neuron coords by a factor (and optionally center)
#'
#' Note that if scale=TRUE, the neuron will be rescaled to unit sd in each axis
#' likewise if center=TRUE, the neuron will be centred around the axis means
#' @param scale 3-vector used to divide x,y,z coords
#' @param center 3-vector to subtract from x,y,z coords
#' @return neuron with scaled coordinates
#' @export
#' @seealso \code{\link{scale.default}}
#' @examples
#' n1.scaledown=scale(MyNeurons[[1]],c(2,2,3))
#' n1.scaleup=scale(MyNeurons[[1]],1/c(2,2,3))
scale.neuron<-function(n,scale,center=F){
	d=xyzmatrix(n)
	ds=scale(d,scale=scale,center=center)
	n$d[,colnames(d)]=ds
	n
}

all.equal.neuron<-function(target,current,tolerance=1e-6,check.attributes=FALSE,
	fieldsToCheck=c("NeuronName", "NumPoints", "StartPoint", "BranchPoints",
		"EndPoints", "NumSegs", "SegList", "d"), fieldsToCheckIfPresent="nTrees",
	CheckSharedFieldsOnly=FALSE, ...){
	if(length(fieldsToCheck)==1 && is.na(fieldsToCheck))
		fieldsToCheck=names(current)
		
	if(!is.neuron(target) || !is.neuron(current))
		return ("target and current must both be neurons")
	fieldsInCommon=intersect(names(target),names(current))
	# figure out which of the optional fields to check are present
	fieldsToCheckIfPresent=intersect(fieldsInCommon,fieldsToCheckIfPresent)
	# and add those to the fields to check 
	fieldsToCheck=unique(c(fieldsToCheck,fieldsToCheckIfPresent))
	if(CheckSharedFieldsOnly){
		fieldsToCheck=intersect(fieldsInCommon,fieldsToCheck)
	} else{
		# check all core fields
		missingfields=setdiff(fieldsToCheck,names(current))
		if(length(missingfields)>0)
			return(paste("Current missing fields: ",missingfields))
		missingfields=setdiff(fieldsToCheck,names(target))
		if(length(missingfields)>0)
			return(paste("Target missing fields: ",missingfields))		
	}
	all.equal(target[fieldsToCheck],current[fieldsToCheck],
		tolerance=tolerance, check.attributes=check.attributes, ...)
}

# so that you can make an empty neuronlist
neuronlist <- function(...) as.neuronlist(list(...))

as.neuron<-function(n){
	if(is.null(n)) return (NULL)
	if(!is.neuron(n,Strict=TRUE)) class(n)=c("neuron",class(n))
	n
}

# so you can just do:
# plot(ANeuron)
plot.neuron<-function(...) plotneuron2d(...)

as.neuronlist<-function(l,df,AddClassToNeurons=TRUE){
	# makes a list of neurons that can be used for 
	# coordinated plots / analysis
	# allow one neuron to be passed
	if(is.neuron(l)) {
		n<-l
		l<-list(n)
		names(l)<-n$NeuronName
	}
	if(!missing(df)) {
		if(nrow(df)!=length(l)) 
			stop("data frame must have same number of rows as there are neurons")
		attr(l,"df")=df
		if(is.null(names(l)))
			names(l)=rownames(df)
		else if(any(names(l)!=rownames(df)))
			stop("mismatch between neuronlist names and dataframe rownames")
	}
	if(!inherits(l,"neuronlist")) class(l)<-c(class(l),"neuronlist")
	if(!AddClassToNeurons) return(l)
	for(i in seq(l)){
		if(!is.neuron(l[[i]],Strict=TRUE))
			l[[i]]=as.neuron(l[[i]])
	}
	l
}

#' Subset a neuronlist returning either a new neuronlist or the names of chosen neurons
#'
#' EITHER use its attached dataframe as the basis of 
#' a subset operation. Then use rownames of the new dataframe to select
#' neuronlist entries and return that sublist
#' OR apply a function to every item in the list 
#' that returns TRUE/FASLE to determine inclusion in output list
#'
#' When ReturnList is F just return the indices into the list
#' 
#' When INDICES are specified, then use a for loop to iterate over only those
#' members of the list. This is equivalent to myneuronlist[INDICES] but is much
#' faster
#'
#' @param nl a neuronlist
#' @param INDICES optional indices to subset neuronlist (faster for big lists)
#' @param ReturnList whether to return the selected neurons (when T) or just their names
#' @param ... either a function or column names in the attached dataframe
#' @export
#' @examples
#' #Apply a 3d search function to the first 100 neurons in the neuronlist dataset
#' subset(dps[1:100],function(x) {length(subset(x,s3d))>0},ReturnList=F)
#' #The same but using INDICES, which is up to 100x faster when neuronlist is large
#' subset(dps,function(x) {length(subset(x,s3d))>0},INDICES=names(dps)[1:100])
subset.neuronlist<-function(nl, ..., INDICES=NULL, ReturnList=is.null(INDICES)){
	arglist=try(pairlist(...),silent=TRUE)
	if(!inherits(arglist,"try-error") && is.function(arglist[[1]])){
		# we are going to apply a function to every element in neuronlist 
		# and expect a return value
		if(length(arglist)>1) stop("I don't know how to handle optional function args.",
			" Use an anonymous function instead")
		if(is.null(INDICES)){
			snl=sapply(nl,arglist[[1]])
			if(ReturnList) return(nl[snl])
			else return(names(nl)[snl])
		} else {
			if(inherits(INDICES,"character")){
				snl=logical(length(INDICES))
				names(snl)=INDICES
			} else if(inherits(INDICES,"logical")){
				snl=logical(sum(INDICES))
				names(snl)=names(nl)[INDICES]
			} else if(inherits(INDICES,"integer")){
				snl=logical(length(INDICES))
				names(snl)=names(nl)[INDICES]
			}
			if(ReturnList) {
				newlist=list()
				for (n in names(snl)){
					include=arglist[[1]](nl[[n]])
					if(include) newlist[[n]]=nl[[n]]
				}
				return(newlist)
			}
			else{
				for (n in names(snl)){
					snl[n]=arglist[[1]](nl[[n]])
				}
				return(names(which(snl)))
			} 
		}
	} else {
		df=attr(nl,'df')
		sdf=subset(df,...)
	}
	if(ReturnList) nl[rownames(sdf)]
	else return(rownames(sdf))
}

"[.neuronlist" <- function(nl,inds,...) {
	attribs=attributes(nl)
	class(nl)='list'
	nl2=nl[inds,...]
	class(nl2)=attribs$class
	df=attr(nl,'df')
	if(!is.null(df)){
		attr(nl2,'df')=df[inds,,...]
	}
	nl2
}
#' 3D plots of the elements in a neuronlist, optionally using a subset expression
#'
#' @param nl a neuron list (where omitted will use MyNeurons as default)
#' @param subset - an expression passed to subset.neuronlist
#' @param ... options passed on to plot3d (such as colours, line width etc)
#' @return value of plot3d 
#' @export
plot3d.neuronlist<-function(nl,subset,...){
	if(!is.neuronlist(nl)){
		subset=nl
		nl=MyNeurons
	}
	if(!missing(subset)) nl=subset(nl,subset)
	invisible(mapply(plot3d,nl,...))
}

#------------------------------------------------------------------------#
# Returns a neuron containing the segment orders for the ordered collection
# of segments, MyNeuron$SegList 
# DON'T KNOW WHO STILL USES THIS!

SegOrders<-function(MyNeuron){
    # This assumes that the segments are in the correct order from a tree
    # traversal starting at the root MyNeuron$StartPoint
    
    #Check if input is sensible
    # i.e. >1 segment && TreeRoot is an EndPoint
    if(  MyNeuron$NumSegs<1 ){
	SegOrderArray=0
    } else if (MyNeuron$NumSegs==1){
	SegOrderArray=1
    } else {
	# Fill the result arrays with -1 
	PointOrderArray<-rep(-1,MyNeuron$NumPoints)
	SegOrderArray<-rep(-1,MyNeuron$NumSegs)
	
	#Set the tree root to order 0
	PointOrderArray[MyNeuron$StartPoint]<-0
	#Set the first segment to order 1
	SegOrderArray[1]<-1
	#Set the other points in that segment to Order 1
	PointOrderArray[MyNeuron$SegList[[1]][-1]]<-1
	
	# Now iterate through the segments setting the segment order
	# to 1 greater than the point order of the head segment
	# which already has an entry in PointOrderArray
	# because of the order the segments are listed in!
	for(i in 2:MyNeuron$NumSegs){
	    #Find the order of the segment head
	    SegHead<-MyNeuron$SegList[[i]][1]
	    SegHeadOrder<-PointOrderArray[SegHead]
	    #Really don't need this, but just in case I've cocked up
	    if(SegHeadOrder<0){
		print(i)
		print(MyNeuron$SegList[1:i])
		stop("Error in SegOrders - SegHeadOrder unknown")
	    }
	    #Set the order of the segment to 1 greater
	    SegOrderArray[i]<-SegHeadOrder+1
	    #Now set the points in that segment as well
	    PointOrderArray[MyNeuron$SegList[[i]][-1]]<-SegHeadOrder+1
	} # end of for(i in 2:MyNeuron$NumSegs
    }
    # Add or replace MyNeuron$SegOrders
    if(is.null(MyNeuron$SegOrders)){
	return(c(MyNeuron,list(SegOrders=SegOrderArray)))
    } else {
        MyNeuron$SegOrders<-SegOrderArray
	return(MyNeuron)
    }
    
} # end of SegOrders<-function(MyNeuron)

# SegLengths calculates the lengths of each segment in a
SegLengths=function(MyNeuron){
    # convert to numeric matrix without row names
    d=matrix(unlist(MyNeuron$d[,c("X","Y","Z")]),ncol=3)
    sapply(MyNeuron$SegList,function(x) seglength(d[x,]))
}

seglength=function(ThisSeg){
    #ThisSeg is an array of x,y and z data points
    #In order to calculate the length
    #Need to find dx,dy,dz
    #Then find sqrt(dx^2+...)
    #Then sum over the path    # nb this will fail for a segment with only 1 point
	if(!inherits(ThisSeg,"matrix")) ThisSeg=data.matrix(ThisSeg)
	ds=diff(ThisSeg)
    Squared.ds<-ds*ds
    sum(sqrt(rowSums(Squared.ds)))	
}

MergeUnconnectedPathsToSingleNeuron<-function(NeuronList){
	# expects a list of neurons which will be joined together
	# these are expected to be unconnected paths
	# The first neuron in the list will be the 'master' neuron
	
	# $ NeuronName   : chr "JIA5Lskeleton2"
	# $ InputFileName: chr "/GD/projects/PN2/tracings/tofixAmira/JIA5Lskeleton2.am"
	# $ CreatedAt    : POSIXct[1:1], format: "2009-02-24 12:36:48"
	# $ NodeName     : Named chr "psnl-jefferis-2.lmb.internal"
	#  ..- attr(*, "names")= chr "nodename"
	# $ InputFileStat:'data.frame':	1 obs. of  10 variables:
	# $ InputFileMD5 : Named chr "5b0179528bbb4bdebb7729c5edc3e3b5"
	#  ..- attr(*, "names")= chr "/GD/projects/PN2/tracings/tofixAmira/JIA5Lskeleton2.am"
	# $ NumPoints    : int 4662
	# $ StartPoint   : num 1
	# $ BranchPoints : int [1:181] 84 106 234 354 378 379 473 479 481 490 ...
	# $ EndPoints    : int [1:183] 1 450 526 574 612 633 649 658 677 678 ...
	# $ NumSegs      : int 373
	# $ SegList      :List of 373
	# $ nTrees       : int 2
	# $ d            :'data.frame':	4662 obs. of  9 variables:
	# $ SubTrees     :List of 2
	# $ OrientInfo   :List of 5

 	if(!is.list(NeuronList) || length(NeuronList)<2 ) stop("Expects a list of 2 or more neurons")
	
	MasterNeuron=NeuronList[[1]]
	NeuronList=NeuronList[-1]
	MasterNeuron$nTrees=1
	MasterNeuron$SubTrees=list()
	MasterNeuron$SubTrees[[1]]=MasterNeuron$SegList
	MasterNeuron$d[,"SubTree"]=1
	
	for (n in NeuronList){
		MasterNeuron$nTrees=MasterNeuron$nTrees+1
		offset=nrow(MasterNeuron$d)
		d=n$d
		d[,"PointNo"]=d[,"PointNo"]+offset
		d[d[,"Parent"]>0,"Parent"]=d[d[,"Parent"]>0,"Parent"]+offset
		d[,"SubTree"]=MasterNeuron$nTrees
		# only keep the columns that both neurons have in common
		commonCols=intersect(colnames(MasterNeuron$d),colnames(d))
		MasterNeuron$d=rbind(MasterNeuron$d[,commonCols],d[,commonCols])
		MasterNeuron$SubTrees[[length(MasterNeuron$SubTrees)+1]]=lapply(n$SegList,'+',offset)
	}
	MasterNeuron$NumPoints=nrow(MasterNeuron$d)
	MasterNeuron
}

ReadMyChullOutput<-function(FileName){
    LinesToSkip<-2
    ColumnNames<-c("Idx","X","Y","Z")

    read.table(FileName, header = FALSE, sep = "", quote = "\"'", dec = ".",
     col.names=ColumnNames, as.is = FALSE, na.strings = "NA",
	skip = LinesToSkip, check.names = TRUE, fill = FALSE,
	strip.white = TRUE, blank.lines.skip = TRUE)
}


PointsInRegion<-function(x,ANeuron,region=c("root","MB","LH"),pointTypes=c("EP","BP")){
		CandPoints=NULL
		region=toupper(region); pointTypes=toupper(pointTypes)
		
		if("MB"%in%region) CandPoints=unique(unlist(ANeuron$SegList[unlist(ANeuron$MBSegNos)]))
		if("LH"%in%region) CandPoints=c(CandPoints,unique(unlist(ANeuron$SegList[unlist(ANeuron$LHSegNos)])))
		if("ROOT"%in%region) CandPoints=c(CandPoints,ANeuron$StartPoint)
		
		if(!"ALL"%in%pointTypes){
				cand2=NULL
				if("EP"%in%pointTypes) cand2=ANeuron$EndPoints
				if("BP"%in%pointTypes) cand2=c(cand2,ANeuron$BranchPoints)
				CandPoints=intersect(CandPoints,cand2)
		}
				
		x[x%in%CandPoints]
}

####################
#                  #
#   plotneuron2d   #
#                  #
####################

plotneuron2d<-function(ANeuron,WithLine=T,NodesOnly=T,EPsOnly=F,BPsOnly=F,WithText=F,NewWindow=F,
		ScalePlot=F,  # ScalePlot is irrelevant now that I have found par(asp=1)
		UseCurPalette=F,WithContour=T, WithContours=F, WithLHEPoint=T,PlotAxes="XY",axes=TRUE,asp=1,
		MainTitle=paste(ANeuron$NeuronName,ANeuron$CellType),RSOrders=F,xlim=NULL,ylim=NULL,AxisDirections=c(1,-1,1),
		Superimpose=F,LineColour=NULL,PointAlpha=1,tck=NA,lwd=par("lwd"),...){
		# NodesOnly=T will plot dots at endpoints and branchpoints
		# BPsOnly will override this and plot only branchpoints
		# PointAlpha specifies alpha value to use for all points
		
		
		
		#    WithLine<-NodesOnly<-WithLHEPoint<-ScalePlot<-WithContour<-T
		#    WithText<-NewWindow<-UseCurPalette<-WithContours<-F;WithLHEPoint<-T
		# PlotAxes<-c("X","Y")
		if (is.character(ANeuron)){
				ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
		}
		if (is.numeric(ANeuron)){
				ANeuron<-MyNeurons[[ANeuron]]
		}
		
		if (!is.list(ANeuron)){
				warning("Cannot understand passed neuron")
				return(F)
		}
		
		# ImageJ etc use the top left hand corner as the origin
		# whereas R uses cartesian co-ords
		# so default is to invert Y
		# Simplest just to multiply
		# data (contours and points by 1 or -1)
		# at the outset
		if(any(AxisDirections!=1)){
				ANeuron$d[,c("X","Y","Z")]=t(t(ANeuron$d[,c("X","Y","Z")])*AxisDirections)
		}
		
		if(PlotAxes=="XY") {PlotAxes<-c("X","Y");NumPlotAxes<-c(1,2)} else
		if(PlotAxes=="YZ") {PlotAxes<-c("Y","Z");NumPlotAxes<-c(2,3)} else
		if(PlotAxes=="XZ") {PlotAxes<-c("X","Z");NumPlotAxes<-c(1,3)} else 
		if(PlotAxes=="ZY") {PlotAxes<-c("Z","Y");NumPlotAxes<-c(3,2)}
		
		OldPalette<-palette()
		if(!UseCurPalette){
				palette(c("black",rainbow(6)))
		}
		
		if (NewWindow) macintosh(w=9.5,h=9.5)
		
		# This is so that the contour points don't fall off the plot
		if((WithContours || WithContour) ){  # Removed  & !is.null(ANeuron$c)
				myxlims<-range(c(ANeuron$d[PlotAxes[1]],ANeuron$c$d[PlotAxes[1]],
								ANeuron$LH$d[PlotAxes[1]],ANeuron$MB$d[PlotAxes[1]]),na.rm = TRUE)
				myylims<-range(c(ANeuron$d[PlotAxes[2]],ANeuron$c$d[PlotAxes[2]],
								ANeuron$LH$d[PlotAxes[2]],ANeuron$MB$d[PlotAxes[2]]),na.rm = TRUE)
		} else {
				myxlims<-range(ANeuron$d[PlotAxes[1]],na.rm = TRUE)
				myylims<-range(ANeuron$d[PlotAxes[2]],na.rm = TRUE)
		}
		
		# if xlim and ylim were set, just use them
		if (!is.null(xlim)){
				myxlims=xlim
		}
		if (!is.null(ylim)){
				myylims=ylim
		}    
		
		if(ScalePlot){
				deltax<-diff(myxlims)
				deltay<-diff(myylims)
				if (deltax>deltay){
						myylims[1]<-myylims[1]-(deltax-deltay)/2
						myylims[2]<-myylims[2]+(deltax-deltay)/2
				} else {
						myxlims[1]<-myxlims[1]-(deltay-deltax)/2
						myxlims[2]<-myxlims[2]+(deltay-deltax)/2
				}
		}
		if(EPsOnly){
				NodesOnly=setdiff(ANeuron$EndPoints,ANeuron$StartPoint)
				# NB Remove start Point  - not really sure in the end whether this is 
				# always treated as an end point or not
				mycols<-rep(rgb(0,1,0,PointAlpha),length(ANeuron$EndPoints))
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else if(BPsOnly){
				NodesOnly<-ANeuron$BranchPoints
				mycols<-rep("red",length(ANeuron$BranchPoints))
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else if(NodesOnly){
				NodesOnly<-c(ANeuron$BranchPoints,ANeuron$EndPoints,ANeuron$StartPoint)
				mycols<-c(rep("red",length(ANeuron$BranchPoints)),
						rep("green",length(ANeuron$EndPoints)),"purple" )
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else { # end if(NodesOnly)
				mycols<-rep("black",ANeuron$NumPoints)
				mycols[ANeuron$BranchPoints]<-"red"
				mycols[ANeuron$EndPoints]<-"green"
				mycols[ANeuron$StartPoint]<-"purple"
				PlottedPoints<-ANeuron$d[,c("PointNo",PlotAxes)]
		} # end if(BPsOnly)
		
		
		# Add the LHAnchorPoint that is the "major" branching point
		# if it exists
		if(!is.null(ANeuron$LHAnchorPoint)){
				PlottedPoints<-rbind(PlottedPoints,
						ANeuron$d[ANeuron$LHAnchorPoint,c("PointNo",PlotAxes)])
				mycols<-c(mycols,"blue")
		}
		
		# Add the AxonLHEP that is the point on the axon
		# closest to the lateral horn entry point defined
		# by Takaki
		if(!is.null(ANeuron$AxonLHEP)){
				PlottedPoints<-rbind(PlottedPoints,
						ANeuron$d[ANeuron$AxonLHEP,c("PointNo",PlotAxes)])
				mycols<-c(mycols,"purple")
		}
		
		#NOW DO THE PLOT
		#cat("myylims =",myylims,"\n")
		if(Superimpose) points(PlottedPoints[,PlotAxes],col=mycols,pch=20,asp=asp,...) 
		else plot(PlottedPoints[,PlotAxes],col=mycols,pch=20,xlim=myxlims,ylim=myylims,
	#		main=MainTitle,asp=1,axes=all(AxisDirections==1),tck=tck,...) 
			main=MainTitle,asp=asp,axes=axes && all(AxisDirections==1),tck=tck,...) 
#		if(!all(AxisDirections==1) && !Superimpose){
		if(axes && !all(AxisDirections==1) && !Superimpose){
				# need to provide special treatment for axes
				box()
				if(AxisDirections[NumPlotAxes][1]!=1){
						axis(1, at=axTicks(1),tck=tck,
								label=axTicks(1)*AxisDirections[NumPlotAxes][1])
				} else {
						axis(tck=tck,1)
				}
				if(AxisDirections[NumPlotAxes][2]!=1){
						axis(2, at=axTicks(2),tck=tck,
								label=axTicks(2)*AxisDirections[NumPlotAxes][2])
				} else {
						axis(tck=tck,2)
				}
		}
		
		if(WithText){
				text(PlottedPoints[,PlotAxes[1]],PlottedPoints[,PlotAxes[2]],
						PlottedPoints[,"PointNo"],col=mycols,pos=3)
		}
		
		if (WithLine){
				# All Black lines
				if(!is.null(LineColour) && length(LineColour==1)){
						MyCols=rep(LineColour,ANeuron$NumSegs)
				} else {
						
						MyCols<-rep(1,ANeuron$NumSegs)
						# Coloured by Orders from root
						if(!is.null(ANeuron$SegOrders)){
								MyCols<-ANeuron$SegOrders
						}
						# Coloured by Segment Types (Axon,LH,MB)
						if(!is.null(ANeuron$SegTypes)){
								MyCols<-ANeuron$SegTypes
						}
						# Coloured according to the reverse Strahler procedure
						if(RSOrders && !(is.null(ANeuron$RSOrders))){
								MyCols[ANeuron$LHSegNos]<-ANeuron$RSOrders+1
						}
						if(!is.null(LineColour)){
							MyCols=LineColour[MyCols]
						}
						
				}
				# NOW PLOT THE LINES!	    
				#for(i in 1:ANeuron$NumSegs){
				#    lines(ANeuron$d[ANeuron$SegList[[i]],PlotAxes]
				#	,col=MyCols[i])
				#}  
				############## 
				# Found a quicker way to do this using NAs
				# to interrupt plotting and start a new line
				# Unfortunately have to do this for each new line colour
				# since lines can only do one line colour at a time
				for(thisCol in unique(MyCols)){
						SegsToPlot=ANeuron$SegList[MyCols==thisCol]
						LinesToPlot=unlist(sapply(SegsToPlot,function(x) c(x,NA)))	
						lines(ANeuron$d[LinesToPlot,PlotAxes],col=thisCol,lwd=lwd)
				}
		} # end if (WithLine) 
		
		
		
		#See if any contours exist
		for(ContSet in list(ANeuron[c("c","LH","MB")])){
				if(is.null(ContSet$d)) next  # Seems to be optional, but better safe
				#Plot all contours if required
				if(WithContours){
						
						ContsToPlot=sort(unique(ContSet$d$ContourID))
						PointsToPlot=rep(0,nrow(ContSet$d)+length(ContsToPlot))
						PointsToPlot=NULL
						for(j in ContsToPlot){
								#ThisContourPoints<-which(ContSet$d$ContourID==j)
								# polygon will join the points (only works in 2D)
								#polygon(ContSet$d[ThisContourPoints,PlotAxes],
								#    type="l",lty="dotted",col="black")
								#
								PointsToPlot=c(PointsToPlot,which(ContSet$d$ContourID==j),NA)
						}
						#TRY A SPEEDIER VERSION
						polygon(ContSet$d[PointsToPlot,PlotAxes],lty="dotted",type='l')
						
				}#endif(WithContours){
				# Plot a boundary contour (2D convex hull)
				if(WithContour){
						polygon(ContSet$d[chull(ContSet$d[,PlotAxes]),PlotAxes])
				}
		}#end for
		
		if(WithLHEPoint & !is.null(ANeuron$MarkerPoints$LHEPoint)){
				points(ANeuron$MarkerPoints$LHEPoint[NumPlotAxes[1]],
						ANeuron$MarkerPoints$LHEPoint[NumPlotAxes[2]],col="purple",pch=22)
		}
		if(!is.null(ANeuron$c$GrandCent)){
				points(ANeuron$c$GrandCent[NumPlotAxes[1]],
						ANeuron$c$GrandCent[NumPlotAxes[2]],col="red",bg="red",pch=22)
		}
		
		palette(OldPalette)
		#nb makes an invisible return - doesn't print on command line
		#but can be stored in a variable.
		invisible(PlottedPoints)
}

# Util function for making a line for RGL to draw
makerglline=function(df){
    if(is.null(dim(df))){
    ldf=length(df)
	if(ldf >2)
	    rep(df,rep(2, ldf))[-c(1, ldf*2)]
	else df
    } else {
    	nrdf=nrow(df)
	if(nrdf>2)
	    df[ rep(1: nrdf,rep(2, nrdf)) [-c(1, nrdf*2)],]
	else df
    }
}


####################
#                  #
#   plotneuron3d   #
#                  #
####################
plotneuron3d<-function(ANeuron,UseCurPalette=F,WithContours=T,WithScale=T,
    # note that rgl is default display mode if rgl library is available
    JustLH=F,ToFile=F,ScaleRotater=T,Colour=NULL,UseRGL=NULL,AxisDirections=c(1,-1,-1),
    # These only work with RGL at the moment
    WithLine=T,WithNodes=T,WithAllPoints=F,WithText=F,PlotSubTrees=T,ClearRGL=T,NeuronList=MyNeurons,...){    
	if (is.null(UseRGL)) {
		UseRGL=require(rgl)
	} else if(any(UseRGL)) {
		if(!require(rgl)) cat("Unable to load RGL library, switching to scatterplot3d\n")
	} 
	
    if(UseRGL && ClearRGL) rgl.clear() # clear RGL buffer
   	OldPalette<-palette()
   	if(!is.null(Colour)){
		# tricky thing is what to do if I just want colour to come out
		# as a number 
		if(!is.numeric(Colour)) {
		    # plot everything in this colour
		    palette(rep(Colour,5)) # nb palette expects >1 colour names
		}
    } else if(!UseCurPalette){
		# RGL background is black
		if(UseRGL) palette(c("white",rainbow(6)))
		else palette(c("black",rainbow(6)))
    }

    if (is.character(ANeuron))
		ANeuron<-NeuronList[[GetNeuronNum(ANeuron)]]

    if (is.numeric(ANeuron))
		ANeuron<-NeuronList[[ANeuron]]
    
    if (!is.list(ANeuron)){
		warning("Cannot understand passed neuron")
		return(F)
    }
    
    if(any(AxisDirections!=1))
		ANeuron$d[,c("X","Y","Z")]=t(t(ANeuron$d[,c("X","Y","Z")])*AxisDirections)
    
    
    # Check to see if we want to produce a rotater file
    if(ToFile!=F){
	# could either come in as true in which case a default name
	# is given to the file or as a string in which case a file
	# of that name is created in RotDir
	if(is.character(ToFile)){
	    OutFile<-file.path(RotDir,ToFile)
	} else {
	    OutFile<-file.path(RotDir,paste(ANeuron$CellType,sep="",".",ImageName(ANeuron$NeuronName),
		    ifelse(ScaleRotater,"_scl.rot",""),ifelse(WithContours,"_wc.rot","_woc.rot")))
	}
	
	# Try creating file
	if(!file.create(OutFile)) stop(paste("Couldn\'t create file",Outfile))
	# Write out some header information
	cat("#",basename(OutFile),"\n",file=OutFile,append=T)
	cat("# created on",date(),"\n",file=OutFile,append=T)
	cat("# Neuron",ImageName(ANeuron$NeuronName),"CellType",ANeuron$CellType,"\n",file=OutFile,append=T)
	cat("# NumPoints",length(ANeuron$d$X),"\n",file=OutFile,append=T)
	cat("# NumContours",ANeuron$c$ContInfo$NumContours,"\n",file=OutFile,append=T)
	# Set this flag so that later routines know to write to file
	ToFile<-T
	
	# For rotater, it's worth setting some useful point as the
	# zero position
    
	if (ScaleRotater){
	    if(!is.null(ANeuron$Scl)){
			ZeroPos<-unlist(ANeuron$c$GrandCent)
			names(ZeroPos)<-c("X","Y","Z")
			RotScl<-ANeuron$Scl
	    }
	    else{
			cat("Can't scale rotater output since",ANeuron$Name,"has no scale information")
			stop("Try sourcing SpatialAnalysis.s to update MyNeurons")
	    }
	}
	
	# Definitions for 16 bit rotator
	# If (Rotater=="Original") {
	DrawDot<--1;DrawMove<-0;DrawLine<-1
	RotDot<-function(colour){return(c(rgbcolour(colour),DrawDot))}
	RotMove<-function(colour){return(c(rgbcolour(colour),DrawMove))}
	RotLine<-function(colour){return(c(rgbcolour(colour),DrawLine))}

	# Original Rotater commands
	DrawDot<--1;DrawMove<-0;DrawLine<-1
	simpCol<-function(colour){
	    switch(colour,'red'=1,'green'=2,'blue'=3,'yellow'=4,'purple'=5,'cyan'=6,7)
	}
	# just revert to using colour numbers
	# My improved version of Rotater can handle up to 17
	simpCol<-function(colour){
	    # nb 0 => no dot, just move
	    if(colour>=0 && colour <=19) return (colour)
	    return (7)
	}    
	
	# NB here are the colour definitions from GJRotater5
########################################################################
#                                                                      #
#           case 1: R=255; G=0; B=0; break; //red                      #
#           case 2: R=0; G=255; B=0; break; //green                    #
#           case 3: R=0; G=0; B=255; break; //blue                     #
#           case 4: R=255; G=255; B=0; break; //yellow                 #
#           case 5: R=255; G=0; B=255; break; //purple                 #
#           case 6: R=0; G=255; B=255; break;  //cyan                  #
#           case 8: R=255; G=127; B=0; break;  //GJ: new colours!      #
#           case 9: R=255; G=0; B=127; break;  //GJ: new colours!      #
#           case 10: R=255; G=127; B=127; break;  //GJ: new colours!   #
#           case 11: R=127; G=255; B=0; break;  //GJ: new colours!     #
#           case 12: R=0; G=255; B=127; break;  //GJ: new colours!     #
#           case 13: R=127; G=255; B=127; break;  //GJ: new colours!   #
#           case 14: R=127; G=0; B=255; break;  //GJ: new colours!     #
#           case 15: R=0; G=127; B=255; break;  //GJ: new colours!     #
#           case 16: R=127; G=127; B=255; break;  //GJ: new colours!   #
#           case 17: R=127; G=255; B=255; break;  //GJ: new colours!   #
#           case 18: R=255; G=127; B=255; break;  //GJ: new colours!   #
#           case 19: R=255; G=255; B=127; break;  //GJ: new colours!   #
#           case 7:                                                    #
#           default: R=255; G=255; B=255; break; //white               #
#                                                                      #
########################################################################

	RotDot<-function(colour){return(simpCol(colour)*DrawDot)}
	RotMove<-function(colour){return(simpCol(colour)*DrawMove)}
	RotLine<-function(colour){return(simpCol(colour)*DrawLine)}

	RotWrite<-function(Points,PointCols,DrawMethod,Contours=F){
	    if(is.null(PointCols)) PointCols<-rep(Colour,length(Points))
	    if(length(PointCols)==1) PointCols<-rep(PointCols,length(Points))
	    
	    if(Contours){
		# GJ: added a routine to handle new style contours
		PointsXYZ=NULL
		for(ContSet in list(ANeuron$c,ANeuron$LH,ANeuron$MB)){
		    if(is.null(ContSet)) next  # Seems to be optional, but better safe
		    PointsXYZ<-c(PointsXYZ,ContSet)
		}
	    } else {
		PointsXYZ<-ANeuron$d
	    }
	    PointCols<<-PointCols
	    # Scale points to have origin at centre of LH
	    # and have lh height width etc as 1.
	    # should perhaps run these measurements
	    # on standard brains in order to figure out
	    # what would be appropriate ratios between these
	    # XYZ axes.
	    if(ScaleRotater){
		#PointsXYZ<-(PointsXYZ-ZeroPos)/RotScl
		PointsXYZ[,c("X","Y","Z")]<-scale(PointsXYZ[,c("X","Y","Z")],center=ZeroPos,scale=RotScl)
	    }
	    
	    write.table( 
		cbind( PointsXYZ[Points,c("X","Y","Z")],(sapply(PointCols,DrawMethod)) ),
		row.names=F,col.names=F,file=OutFile,append=T )
	}
    }
        
    Scl<-c(1,1,1)
    
    # UNFINISHED _ NEED TO PUT if(JustLH) in the plot routines
    if(!UseRGL){
	if(JustLH){
	    LHSegList<-ANeuron$SegList[ANeuron$LHSegNos]
	    LHPoints<-unique(unlist(LHSegList))	    
	    myxlims<-range(Scl[1]*c(ANeuron$d$X[LHPoints],ANeuron$c$d$X))
	    myylims<-range(Scl[2]*c(ANeuron$d$Y[LHPoints],ANeuron$c$d$Y))
	    myzlims<-range(Scl[3]*c(ANeuron$d$Z[LHPoints],ANeuron$c$d$Z)) 
	    
	} else {
	    # We want the whole of the neuron
	    if(WithContours){
		# UPDATED THIS TO CHECK FOR OLD AND NEW CONTOUR INFO
		myxlims<-range(Scl[1]*c(ANeuron$d$X,ANeuron$c$d$X,ANeuron$LH$d$X,ANeuron$MB$d$X))
		myylims<-range(Scl[2]*c(ANeuron$d$Y,ANeuron$c$d$Y,ANeuron$LH$d$Y,ANeuron$MB$d$Y))
		myzlims<-range(Scl[3]*c(ANeuron$d$Z,ANeuron$c$d$Z,ANeuron$LH$d$Z,ANeuron$MB$d$Z)) 
	    } else {
		myxlims<-range(Scl[1]*c(ANeuron$d$X))
		myylims<-range(Scl[2]*c(ANeuron$d$Y))
		myzlims<-range(Scl[3]*c(ANeuron$d$Z)) 
	    }
	}
    
	# I would have liked to use asp=1, but that parameter
	# doesn't seem to work for scatterplot3d
	if(WithScale){
	    AxisRanges<-abs(diff(cbind(myxlims,myylims,myzlims)))
	    MaxRange<-max(AxisRanges)
	    Deltas<-diff(rbind(AxisRanges,rep(MaxRange,3)))
	    myxlims<-myxlims+c(-Deltas[1]/2,+Deltas[1]/2)
	    myylims<-myylims+c(-Deltas[2]/2,+Deltas[2]/2)
	    myzlims<-myzlims+c(-Deltas[3]/2,+Deltas[3]/2)
	}
    }

    if(!ToFile || !is.numeric(Colour)){
	# Get the colours right to plot the root and EndPoints
	NodesOnly<-c(ANeuron$BranchPoints,
	    ANeuron$EndPoints[-which(ANeuron$EndPoints==ANeuron$StartPoint)],
	    ANeuron$StartPoint)
	NodeCols<-c(rep("red",length(ANeuron$BranchPoints)),
	    rep("green",length(ANeuron$EndPoints)-1),"purple" )
	
	if(!is.null(ANeuron$AxonLHEP)){
	    # Check if we've already put the LHEP into the array
	    if (any(NodesOnly==ANeuron$AxonLHEP)){
		NodeCols[which(NodesOnly==ANeuron$AxonLHEP)]<-"purple"
	    } else {
		NodesOnly<-c(NodesOnly,ANeuron$AxonLHEP)
		NodeCols<-c(NodeCols,"purple")
	    }
	}# end  if(!is.null(ANeuron$AxonLHEP)){
	
	# Add the LHAnchorPoint (i.e. main branch if it exists)
	if(!is.null(ANeuron$LHAnchorPoint)){
	    NodeCols[which(NodesOnly==ANeuron$LHAnchorPoint)]<-"blue"
	}# end  if(!is.null(ANeuron$LHAnchorPoint)){
    } else {
	# this is assuming colours come in as numbers
	NodesOnly<-c(ANeuron$BranchPoints,
	    ANeuron$EndPoints[-which(ANeuron$EndPoints==ANeuron$StartPoint)],
	    ANeuron$StartPoint)
	NodeCols<-c(rep(0,length(ANeuron$BranchPoints)),
	    rep(Colour,length(ANeuron$EndPoints)-1),0 )
	
	if(!is.null(ANeuron$AxonLHEP)){
	    # Check if we've already put the LHEP into the array
	    if (any(NodesOnly==ANeuron$AxonLHEP)){
		NodeCols[which(NodesOnly==ANeuron$AxonLHEP)]<-0
	    } else {
		NodesOnly<-c(NodesOnly,ANeuron$AxonLHEP)
		NodeCols<-c(NodeCols,0)
	    }
	}# end  if(!is.null(ANeuron$AxonLHEP)){
	
	# Add the LHAnchorPoint (i.e. main branch if it exists)
	if(!is.null(ANeuron$LHAnchorPoint)){
	    NodeCols[which(NodesOnly==ANeuron$LHAnchorPoint)]<-0
	}# end  if(!is.null(ANeuron$LHAnchorPoint)){

    }

    #Just nodes of neuron
	if(UseRGL){
		if(WithNodes){
			if(!WithLine) NodeCols=rep(Colour,length(NodeCols))
			rgl.points(Scl[1]*ANeuron$d$X[NodesOnly],
				Scl[2]*ANeuron$d$Y[NodesOnly],
				Scl[3]*ANeuron$d$Z[NodesOnly],color=NodeCols,size=3)
			if(WithText) # text labels for nodes
			rgl.texts(Scl[1]*ANeuron$d$X[NodesOnly],
				Scl[2]*ANeuron$d$Y[NodesOnly],
				Scl[3]*ANeuron$d$Z[NodesOnly],NodesOnly,color=NodeCols,adj=c(0,0.5))
			
		}
		
	} else if(!ToFile){	
		My3DPlot<-scatterplot3d(Scl[1]*ANeuron$d$X[NodesOnly],
			Scl[2]*ANeuron$d$Y[NodesOnly],
			Scl[3]*ANeuron$d$Z[NodesOnly],
			xlim=myxlims,ylim=myylims,zlim=myzlims,
			color=NodeCols,pch=20,main=paste(ANeuron$NeuronName,ANeuron$CellType))
	} else {
		RotWrite(NodesOnly,NodeCols,RotDot)
		cat("# End of Nodes\n",file=OutFile,append=T)
	}
	
	
    
    # all points (in white, for RGL)
    if(UseRGL && WithAllPoints){
	rgl.points(Scl[1]*ANeuron$d$X[-NodesOnly],
	    Scl[2]*ANeuron$d$Y[-NodesOnly],
	    Scl[3]*ANeuron$d$Z[-NodesOnly],color='white',size=2)
    }

    #Just neuron lines
	if(WithLine){
		if(is.null(ANeuron$SegTypes)) ANeuron$SegTypes=rep(1,ANeuron$NumSegs)
		if(UseRGL){
			if(PlotSubTrees && !is.null(ANeuron$nTrees) && ANeuron$nTrees>1){
				# handle plotting of subtrees in different colours
				for(i in 1:ANeuron$nTrees){
					pointIndexes=unlist(sapply(ANeuron$SubTrees[[i]],makerglline))
					rgl.lines(Scl[1]*ANeuron$d$X[pointIndexes],
						Scl[2]*ANeuron$d$Y[pointIndexes],
						Scl[3]*ANeuron$d$Z[pointIndexes],col=rainbow(ANeuron$nTrees)[i],...)
				}
			}  else {
				# just 1 tree
				pointIndexes=unlist(sapply(ANeuron$SegList,makerglline))
				if(!is.numeric(Colour) && !is.null(Colour)){
					cols=Colour
				} else {
					cols= rep(ANeuron$SegTypes,sapply(ANeuron$SegList,function(x) ifelse(length(x)>2,2*length(x)-2,length(x))))
					if(any(is.na(cols))) cols='red'
					else cols=palette()[cols]
					# note that if cols is passed as a number then rgl.lines
					# flashes up a regular graphics window for some reason
				}
				rgl.lines(Scl[1]*ANeuron$d$X[pointIndexes],
					Scl[2]*ANeuron$d$Y[pointIndexes],
					Scl[3]*ANeuron$d$Z[pointIndexes],col=cols,...)
			}

		} else {
			for (j in 1:ANeuron$NumSegs){
				ThisSegPoints<-ANeuron$SegList[[j]]
				if(UseRGL){
					rgl.lines(makerglline(ANeuron$d$X[ThisSegPoints]),
						makerglline(ANeuron$d$Y[ThisSegPoints]),
						makerglline(ANeuron$d$Z[ThisSegPoints]),col=ifelse(is.numeric(Colour),ANeuron$SegTypes[j],Colour))
				}
				else if(!ToFile){
					My3DPlot$points3d(Scl[1]*ANeuron$d$X[ThisSegPoints],
						Scl[2]*ANeuron$d$Y[ThisSegPoints],
						Scl[3]*ANeuron$d$Z[ThisSegPoints],col=ANeuron$SegTypes[j],type="l"
					)
				} else {
					RotWrite(ThisSegPoints[1],ANeuron$SegTypes[j],RotMove)
					RotWrite(ThisSegPoints[-1],ANeuron$SegTypes[j],RotLine)
				}
			}
			
		}
	} # end of if(WithLine){

    #Just Contours
    if(!UseRGL && WithContours){
	# This little extra loop turns out to be a nice way to deal
	# with the uncertainty of which type of contour information
	# will be present
		for(ContSet in list(ANeuron$c,ANeuron$LH,ANeuron$MB)){
		    if(is.null(ContSet)) next  # Seems to be optional, but better safe
		    for(j in unique(ContSet$d$ContourID)){
			ThisContourPoints<-which(ContSet$d$ContourID==j)
			# Just to join the circle
			ThisContourPoints<-c(ThisContourPoints,ThisContourPoints[1])
			if(!ToFile){
			    My3DPlot$points3d(Scl[1]*ContSet$d$X[ThisContourPoints],
				Scl[2]*ContSet$d$Y[ThisContourPoints],
				Scl[3]*ContSet$d$Z[ThisContourPoints]
				,type="l",lty="dotted")
			} else {
			    RotWrite(ThisContourPoints[1],'white',RotMove,Contours=T)
			    RotWrite(ThisContourPoints[-1],'white',RotDot,Contours=T)
			}
		    }
	    
		}# end for ContSet
    }
    palette(OldPalette)
    invisible(T)
    
}

#################
#               #
#   rgbcolour   #
#               #
#################
# function to return an RCB appropriate colour 3 vector
# if given an integer or a defined colour name string
rgbcolour<-function(colour){
    if(is.character(colour)){
	return(switch(colour, red=c(31,0,0), green=c(0,31,0),
	blue=c(0,0,31),purple=c(31,0,31),cyan=c(0,31,31),yellow=c(31,31,0),
	peach=c(31,15,15),eight=c(15,31,15),nine=c(15,15,31),
	ten=c(31,15,0),eleven=c(31,0,15),twelve=c(15,0,31),
	thirteen=c(15,31,0),fourteen=c(0,15,31),white=c(31,31,31) ))
    } else {   
	if(colour<1) return(rgbcolour('white'))
	return( 
rgbcolour(switch(colour,'red','green','blue','purple','cyan','yellow',
	    'peach','eight','nine','ten','eleven',
	    'twelve','thirteen','fourteen','white')) )
    }
}


#######################
#                     #
#   plotendpoint3d   #
#                     #
#######################
# function to plot only scaled endpoints directly
# to rotater files
# So far this only works for LH data (whether it is stored in $c or $LH)
plotendpoint3d<-function(ANeuron,UseCurPalette=F,WithContours=F,WithScale=T,
    JustLH=F,ToFile=T,ScaleRotater=T,ThisCol=1,FileAppend=F){
 
    if (is.character(ANeuron)){
	ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
    }
    if (is.numeric(ANeuron)){
	ANeuron<-MyNeurons[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
	warning("Cannot understand passed neuron")
	return(F)
    }
    
    # could either come in as true in which case a default name
    # is given to the file or as a string in which case a file
    # of that name is created in RotDir
	if(is.character(ToFile)){
		# check if ToFile is a full path or just a filename
		if(dirname(ToFile)=="."){	    
			OutFile<-file.path(RotDir,ToFile)
		} else {
			#check to see if the path is correct
			if (file.exists(dirname(ToFile))){
				OutFile<-ToFile
			} else {
				stop(paste("Couldn't find the directory referenced in the supplied ToFile",ToFile))
			}
		}
	} else {
		OutFile<-file.path(RotDir,paste(ANeuron$CellType,sep="",".",ImageName(ANeuron$NeuronName),
				".end",ifelse(ScaleRotater,".scl",""),ifelse(WithContours,".wc",".woc")))
	}
	
	
    # Try creating file if there isn't already one open
    if(!FileAppend){
	if(!file.create(OutFile)) stop(paste("Couldn\'t create file",Outfile))
    }
    
    # Write out some header information
    cat("#",basename(OutFile),"\n",file=OutFile,append=T)
    cat("# created on",date(),"\n",file=OutFile,append=T)
    cat("# Neuron",ImageName(ANeuron$NeuronName),"CellType",ANeuron$CellType,"\n",file=OutFile,append=T)
    cat("# NumPoints",length(ANeuron$d$X),"\n",file=OutFile,append=T)
    cat("# NumContours",ANeuron$c$ContInfo$NumContours,"\n",file=OutFile,append=T)
    # Set this flag so that later routines know to write to file
    ToFile<-T
    
    # For rotater, it's worth setting some useful point as the
    # zero position

    if (ScaleRotater){
	if(!is.null(ANeuron$Scl)){
	    ZeroPos<-unlist(ANeuron$c$GrandCent)
	    names(ZeroPos)<-c("X","Y","Z")
	    RotScl<-ANeuron$Scl
	}
	else{
	    cat("Can't scale rotater output since",ANeuron$Name,"has no scale information")
	    stop("Try sourcing SpatialAnalysis.s to update MyNeurons")
	}
    }
    
    # Definitions for 16 bit rotator
    DrawDot<--1;DrawMove<-0;DrawLine<-1
    RotDot<-function(colour){return(c(rgbcolour(colour),DrawDot))}
    RotMove<-function(colour){return(c(rgbcolour(colour),DrawMove))}
    RotLine<-function(colour){return(c(rgbcolour(colour),DrawLine))}

    RotWrite<-function(Points,PointCols,DrawMethod,Contours=F){

	if(length(PointCols)==1) PointCols<-rep(PointCols,length(Points))
	if(Contours){
	    PointsXYZ<-ANeuron$c$d
	} else {
	    PointsXYZ<-ANeuron$d
	}
	# Scale points to have origin at centre of LH
	# and have lh height width etc as 1.
	# should perhaps run these measurements
	# on standard brains in order to figure out
	# what would be appropriate ratios between these
	# XYZ axes.
	if(ScaleRotater){
	    #PointsXYZ<-(PointsXYZ-ZeroPos)/RotScl
	    PointsXYZ[,c("X","Y","Z")]<-scale(PointsXYZ[,c("X","Y","Z")],center=ZeroPos,scale=RotScl)
	}
	
	write.table(
	    cbind(PointsXYZ[Points,c("X","Y","Z")],
		t(sapply(PointCols,DrawMethod)) ),
	    row.names=F,col.names=F,file=OutFile,append=T)
    }
	
    # Get the colours right to plot the root and EndPoints
    LHPoints<-unique(unlist(ANeuron$SegList[ANeuron$LHSegNos]))
    LHEndPoints<-ANeuron$EndPoints[sapply(ANeuron$EndPoints,function(x){any(LHPoints==x)})]

    NodesOnly<-LHEndPoints
    NodeCols<-rep(ThisCol,length(NodesOnly))
    
    RotWrite(NodesOnly,NodeCols,RotDot)
    cat("# End of Nodes\n",file=OutFile,append=T)
     
     #Just Contours
    if(WithContours){
	# This little extra loop turns out to be a nice way to deal
	# with the uncertainty of which type of contour information
	# will be present
	for(ContSet in list(ANeuron$c,ANeuron$LH)){
	    if(is.null(ContSet)) next  # Seems to be optional, but better safe
	    for(j in unique(ContSet$d$ContourID)){
		ThisContourPoints<-which(ContSet$d$ContourID==j)
		# Just to join the circle
		ThisContourPoints<-c(ThisContourPoints,ThisContourPoints[1])
		if(!ToFile){
		    My3DPlot$points3d(Scl[1]*ContSet$d$X[ThisContourPoints],
			Scl[2]*ContSet$d$Y[ThisContourPoints],
			Scl[3]*ContSet$d$Z[ThisContourPoints]
			,type="l",lty="dotted")
		} else {
		    RotWrite(ThisContourPoints[1],'white',RotMove,Contours=T)
		    RotWrite(ThisContourPoints[-1],'white',RotDot,Contours=T)
		}
	    }

	}# end for ContSet
    }
    
    return(OutFile)
    
}



# Wrapper function to call plotneuron2d
# if given a set of neurons (names or numbers)
plotneurons2d<-function(NeuronRef,MultiPlot=F,Ask=!MultiPlot,LineColours=NULL,...){
    # If there are several neurons to plot, it makes sense to pause
    # unless we have multiplot situation when we want to superimpose
    # neurons
    oldpar<-par(ask=Ask)
    # the ... should allow any additional arguments to be passed to plotneuron2d
    if(length(LineColours)!=length(NeuronRef)) LineColours=rep(LineColours,length(NeuronRef))
    if(length(NeuronRef)>1 && MultiPlot){
	plotneuron2d(NeuronRef[1],LineColour=LineColours[1],...)
 	for(i in 2:length(NeuronRef))
	    plotneuron2d(NeuronRef[i],Superimpose=T,LineColour=LineColours[i],...)
    } else for(i in 1:length(NeuronRef)){
	plotneuron2d(NeuronRef[i],LineColour=LineColours[i],...)
    }
    par(oldpar)
}

# Wrapper function to call plotneuron3d
# if given a set of neurons (names or numbers)
plotneurons3d<-function(NeuronRef,Ask=F,ToFile=F,Colours=NULL,UseRGL=TRUE,NeuronList=MyNeurons,...){
    # If there are several neurons to plot, it makes sense to pause
    if(!ToFile && !any(UseRGL)) oldpar<-par(ask=Ask)
    
    # the ... should allow any additional arguments to be passed to plotneuron2d
    if(!is.null(Colours)){
	if(length(Colours)!=length(NeuronRef)) Colours=rep(Colours[1],length(NeuronRef))
	for(i in 1:length(Colours)) plotneuron3d(NeuronRef[i],ToFile=ToFile,Colour=Colours[i],UseRGL=UseRGL,NeuronList=NeuronList,...)
    } else {
	t<-sapply(NeuronRef,plotneuron3d,ToFile=ToFile,UseRGL=UseRGL,NeuronList=NeuronList,...)
    }
    if(!ToFile && !any(UseRGL)) par(oldpar)
}

# Wrapper function to call plotendpoint3d
# if given a set of neurons (names or numbers)
plotendpoints3d<-function(NeuronRef,Ask=T,ToFile=T,Col=1,FileAppend=T,...){
    # If there are several neurons to plot, it makes sense to pause
    # the ... should allow any additional arguments to be passed to plotneuron2d
    if(length(Col)==1){
	Col<-rep(Col,length(NeuronRef))
    }
    
    if(FileAppend==T && length(NeuronRef)>1){
	FirstOutfile<-plotendpoint3d(NeuronRef[1],ToFile=ToFile,FileAppend=T,ThisCol=Col[1],...)
	if(file.exists(FirstOutfile)){
	    # all well
	    NeuronRef<-NeuronRef[-1]
	    ToFile<-FirstOutfile
	    Col<-Col[-1]
	} else {
	    stop(paste("Error writing file",FirstOutfile))
	}
	for(i in 1:length(NeuronRef)){
	    t<-
	    plotendpoint3d(NeuronRef[i],ToFile=ToFile,FileAppend=FileAppend,ThisCol=Col[i],...)
	}
	
    }
}

# Handy little function to return the MyNeurons
# subscript for a given name - see pmatch for details
# of partial matching
GetNeuronNum<-function(Nnames,mask=1:length(MyNeurons)){
    CellNames<-NULL
    for(i in mask){
	CellNames[i]<-MyNeurons[[i]]$NeuronName
    }
    Nnum<-pmatch(toupper(Nnames),toupper(CellNames),nomatch=0)
    return(Nnum)
}
GetNeuronName<-function(Nnum,mask=1:length(MyNeurons)){
    CellNames<-NULL
    for(i in mask){
	CellNames[i]<-MyNeurons[[i]]$NeuronName
    }
    Nnames<-CellNames[Nnum]
    return(Nnames)
}


# Handy little function to return the MyNeurons
# entry for a given name or number
# nb can only return one neuron
GetNeuron<-function(NeuronRef,mask=1:length(MyNeurons)){    
    if (is.character(NeuronRef)){
	Nnum<-GetNeuronNum(NeuronRef,mask)
    } else {
	# it was already a number
	Nnum<-NeuronRef
    }
    return(MyNeurons[[mask[Nnum]]])
}

#Copied directly from the match help file
intersect <- function(x, y) y[match(x, y, nomatch = 0)]
# Find all the elements of x which have a match in y
matchingyinx <- function(x, y) which(match(x,y,nomatch=0)!=0)

GetNeuronNumsofType<-function(CellTypesToMatch,mask=1:length(MyNeurons)){
    CellTypesToMatch=as.character(CellTypesToMatch)
    CellTypes=sapply(MyNeurons[mask],function(x) x$CellType)
    matchingyinx(toupper(CellTypes),toupper(CellTypesToMatch))
}

GetCellType<-function(NeuronRef,mask=1:length(MyNeurons)){
    if (is.list(NeuronRef)){
	return(NeuronRef$CellType)
    } else {
        # assume that NeuronRef was a number or char
	if (is.character(NeuronRef)){
	    Nnum<-GetNeuronNum(NeuronRef,mask)
	} else {
	    # it was already a number
	    Nnum<-NeuronRef
	}
    }
    
    CellTypes<-NULL
    for(i in mask){
	CellTypes[i]<-MyNeurons[[i]]$CellType
    }
    return(CellTypes[Nnum])
}
ImageName<-function(NeuronName){
	return(unlist(strsplit(NeuronName,"[._]"))[1])
}


FirstnNeurons<-function(n=1){
    # return an array containing the index numbers of the first n neurons
    # of each class
    CellType<-NULL
    for(i in 1:length(MyNeurons)){
	CellType[i]<-MyNeurons[[i]]$CellType
    }
    
    t<-table(CellType)
    if(n>min(t)) stop("That\'s more than the least numerous cell type")
    as.vector(sapply(sort(unique(CellType)),function(x){which(CellType==x)[1:n]}))
}

NeuronNameFromFileName<-function(FileName){
	if(length(FileName)>1) return(sapply(FileName,NeuronNameFromFileName))
    # Get the name of the neuron NB strsplit returns a list)
    MyNeuronName<-unlist(strsplit(basename(FileName),"[._]"))[1]
    # Check that a sensible name resulted
    if(length(MyNeuronName)==0) stop(paste("Invalid neuron name generated from file",FileName))
    return(MyNeuronName)
}

# Guesses the likely input path of a neuron
# based on its input file name and Cell Type and the current setting
# of TraceFileDir (from Startup.s)
InputFilePath<-function(ANeuron){
    if (is.character(ANeuron)){
	ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
    }
    if (is.numeric(ANeuron)){
	ANeuron<-MyNeurons[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
	warning("Cannot understand passed neuron")
	return(F)
    }
    ThePath<-file.path(TraceFileDir,ANeuron$CellType,paste("traced",ANeuron$CellType))
    PathandName<-file.path(ThePath,ANeuron$InputFileName)
    return(PathandName)
}

plotall3=function(ANeuron,...){
    oldpar=par('mfrow')
    par(mfrow=c(1,3))
    plotneuron2d(ANeuron,PlotAxes='XY',...)
    plotneuron2d(ANeuron,PlotAxes='XZ',MainTitle="From above; anterior up, medial left")
    plotneuron2d(ANeuron,PlotAxes='YZ',MainTitle="From the side; anterior up, ventral left")
    par(mfrow=oldpar)
}

plot15=function(recs,...){
    if(length(recs)>15) recs=recs[1:15]
    par(mfrow=c(3,5))
    plotneurons2d(recs,...)
    par(mfrow=c(1,1))
}
plot15g=function(CellType,...){
    plot15(GetNeuronNumsofType(CellType),...)
}


# directly returns filtered vector
# doesn't fall over when pattern doesn't match anything
# 2005-02-03
mygrep=function(pattern,x,keep=TRUE,...){
    found=grep(pattern,x,...)
    if(any(found)){
	if(keep){
	    x[found]
	} else {
	    x[-found]
	}
    } else return(x)
}

getID=function(fileNames){
    if(is.factor(fileNames)) fileNames=as.character(fileNames)
    # first trim either side of first & 2nd underscores
    x=gsub("^[^_]*_([^_]*)_[^_]*.?*$","\\1",fileNames)
    # May still have something left if only 1 underscore
    x=basename(x)
    x=gsub("^([^_]*)[._].*$","\\1",x)
    x=toupper(gsub("^([A-Z]{2,3}[1-9]{1}[0-9]{0,2}[RLTB]([1-4]|LH)).*$","\\1",x,ignore.case=TRUE))
    # if the filename wasn't in a recognised format, it will likely
    # still end in 01, 02 or perhaps 03. If this is the case, remove
    # that terminal digit pair
    gsub("^(.*)0[1-3]$","\\1",x)
}
getBrain=function(fileNames){
    # trim off image number
    x=getID(fileNames)
    x=gsub("^([A-Z]{2,3}[1-9]{1}[0-9]{0,2}[RLTB])([1-4]|LH){0,1}.*$","\\1",x,ignore.case=TRUE)
    # some of Chris' ones will just end in LH
    x=gsub("^(.*)LH$","\\1",x)
	# some will still not conform to this but have a hyphenated terminus
	# eg NP6099MARCMHS35-RESIZED
	gsub("^(.?*)\\-.*$","\\1",x)

}

ReadNeuronFromSWC<-function(f){
	d=ReadSWCFile(f)
	ParseSWCTree(d,f)
}

read.neuron<-function(f, ...){
	# generic function to read in neuron from any kind of file we know about
	# should return exactly one neuron on success
	if(!file.exists(f)) stop("Unable to read file: ",f)
	ext=tolower(sub(".*\\.([^.]+$)","\\1",basename(f)))
	if(ext=="asc")
		n=ReadNeuronFromAsc(f, ...)
	else if(ext=="swc")
		n=ReadNeuronFromSWC(f, ...)
	else {
		h=readLines(f,1)
		if(regexpr("amira",h,ignore.case=TRUE)>0)
			n=ReadNeuronFromAM3D(f, ...)
		else if(regexpr("xml",h,ignore.case=TRUE)>0)
			n=ReadNeuronsFromLongairTraces(f, MergePaths=TRUE, ...)
		else if(regexpr("^;",h)>0)
			n=ReadNeuronFromAsc(f, ...)
		else stop("Unable to identify file format of file: ",f)
	}
	
	if(is.neuron(n,Strict=FALSE)) as.neuron(n)
	else n
}

read.neurons<-function(paths, patt, OmitFailures=TRUE,
	neuronnames=basename, ...){
	if(!is.character(paths)) stop("Expects a character vector of filenames")
	
	if(length(paths)==1 && file.info(paths)$isdir){
		paths=if(missing(patt))
			dir(paths,full=TRUE)
		else
			dir(paths,patt=patt,full=TRUE)
	}
	
	nl=neuronlist()
	for(f in paths){
		n=try(read.neuron(f))
		if(inherits(n,'try-error')) n=NA
		nl[[length(nl)+1]]=n
	}
	if(is.function(neuronnames))
		names(nl)=neuronnames(paths)
	else
		names(nl)=neuronnames
	if(OmitFailures) nl=nl[!is.na(nl)]
	nl
}
