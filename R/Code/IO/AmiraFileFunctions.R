# AmiraFileFunctions.R

# 2005-02-03
# Functions to parse AmiraMesh 3D format - the native
# ouput of the skeletonize plugin and to read and write the density
# data in Amira file formats.
# At the moment depends on SWCFunctions.R since the amiramesh 
# data format can be easily converted to SWC.  However there
# is some redundancy in this approach

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

# source(file.path(CodeDir,"AmiraFileFunctions.R"))

require(tools) # for md5sum

ReadAmiramesh<-function(filename,DataSectionsToRead=NULL,Verbose=FALSE,AttachFullHeader=FALSE,Simplify=TRUE,endian=NULL){
  .Deprecated('read.amiramesh','nat')
  nat::read.amiramesh(file=filename,sections=DataSectionsToRead,
    header=AttachFullHeader,simplify=Simplify,Verbose=Verbose,endian=endian)
}

ReadAmiramesh.Header<-function(con,Verbose=TRUE,CloseConnection=NULL){
  .Deprecated('read.amiramesh.header','nat')
  
  if( is.logical(CloseConnection) && xor(is.character(con), CloseConnection) )
    stop("If con is a character vector, CloseConnection must now be TRUE.",
        " If not CloseConnection must now be FALSE")
  nat::read.amiramesh.header(file=con,Verbose=Verbose)
}

ReadAM3DData<-function(filename,OmitNAs=TRUE){
	# function to read in the basic data from the 
	# files produced by the Amira Skeletonize plugin
	
	ndata=ReadAmiramesh(filename)
	required_fields=c("Coordinates", "NeighbourCount", "Radii", "NeighbourList")
	missing_fields=setdiff(required_fields,names(ndata))
	if(length(missing_fields))
		stop("Neuron: ",filename," is missing fields: ",paste(missing_fields,collapse=" "))
	
	d=data.frame(ndata$Coordinates)
	colnames(d)=c("X","Y","Z")
	d$W=ndata$Radii*2
	d$NeighbourCount=ndata$NeighbourCount
	nVertices=nrow(d)
	d$PointNo=seq(nVertices)
	d[,1:4]=round(d[,1:4],digits=3)
	# Note these numbers come in zero indexed, but I will want them 1-indexed
	Neighbours=data.frame(Neighbour=ndata$NeighbourList+1,CurPoint=rep(seq(nVertices),d$NeighbourCount))
	Origin=ndata$Origins
	if(!is.null(Origin)) Origin=Origin+1
	
	NeuronData=list(PointList=d,EdgeList=Neighbours,Origin=Origin)
	
	# Handle materials
	if(length(ndata$vertexTypeList)){
		SegmentProps=data.frame(
				PointNo=rep(seq(nVertices),ndata$vertexTypeCounter),
				Id=ndata$vertexTypeList)
		if(length(attr(ndata,'Materials'))){
			# we have named materials, so let's add them as an attribute
			attr(SegmentProps,'Materials')=attr(ndata,'Materials')
		}
		# we can only add one numeric label to the SWC format version of the neuron
		FirstSegmentProps=subset(SegmentProps,!duplicated(PointNo))
		# 0 = undefined
		NeuronData$PointList$Label=0L
		NeuronData$PointList[FirstSegmentProps$PointNo,"Label"]=FirstSegmentProps$Id
		NeuronData=c(NeuronData,list(SegmentProps=SegmentProps))
	}
	
	if(OmitNAs) return(RemoveInvalidPointsFromNeighbourList(NeuronData))
	else NeuronData
}

#' Remove invalid points from data block read from Amira skeletonize neuron
#' @param x List containing PointList, EdgeList and Origin
#' @return List with same components as input
#' @author jefferis
#' @seealso \link{\code{ReadAM3DData}}
#' @examples
RemoveInvalidPointsFromNeighbourList<-function(x){
	# check if there are actually any NAs - these should only be in XYZ
	InvalidPoints=which(apply(x$PointList[,c("X","Y","Z")],1,function(p) any(is.na(p))))
	
	if(length(InvalidPoints)>0){
		# Remove all edges that reference an invalid point
		x$EdgeList=subset(x$EdgeList,
				!(CurPoint%in%InvalidPoints | Neighbour%in%InvalidPoints) )
		# find the remaining valid points NB valid points is not necessarily
		# the exact complement of invalid points - consider an end point which 
		# remains defined just distal to a point that doesn't transform.
		ValidPoints=sort(unique(x$EdgeList$CurPoint))
		
		# Restrict d to the valid points
		x$PointList=x$PointList[ValidPoints,]
		x$PointList$OldPointNo=x$PointList$PointNo # keep track of the old point numbers
		x$PointList$PointNo=seq(len=nrow(x$PointList)) # make the new ones
		
		# TODO figure out how to find the closest end point
		# to the original origin - hmm I think a better idea
		# would be to reroot the tree on the origin
		# For the moment, just assume that the origin is not deleted
		# or make it null if it is
		dO=subset(x$PointList,x$Origin==OldPointNo)
		if(nrow(dO)==1) x$Origin=dO$PointNo else x$Origin=NULL
		
		# figure out the new number of neighbours for remaining points
		x$PointList$NeighbourCount=table(x$EdgeList$CurPoint)
		
		# Renumber Neighbours data frame
		x$EdgeList$CurPoint=rep(x$PointList$PointNo,x$PointList$NeighbourCount)
		# Get the new neighbours by looking up against the old points
		# NB could have been written d$PointNo[match(...)]
		# but unnecessary since d$PointNo=seq(nrow(d))
		NewNeighbours=match(x$EdgeList$Neighbour,x$PointList$OldPointNo)
		# Verify that we have non-zero (or NA) entries
		stopifnot(all(NewNeighbours))
		x$EdgeList$Neighbour=NewNeighbours
	}
	return(x)
}

#' Convert data Point and Edge data to the core components of a neuron
#' 
#' The Point and Edge data is essentially what is saved by the Amira 
#' skeletonize format, but the same data can also easily be generated
#' for SWC data and is also pretty much identical to what is
#' required for graph.data.frame from the igraph package.
#' @param PointData Coordinates and another information relating to Points
#' @param Neighbours Reciprocal edges for all connected points (2 col matrix)
#' @param Origin Optional root for tree
#' @param ProcessAllTrees Process all trees if graph is not connected (def T)
#' @param Verbose Whether to show detailed information about progress
#' @return a neuron including both SegList and SWC format data
#' @author jefferis
#' @export
#' @seealso \code{\link{ParseAM3DToNeuron},\link{ParseEdgeListForAllSubTrees}}
CoreNeuronFromPointAndEdgeData<-function(PointData,Neighbours,Origin=NULL,ProcessAllTrees=TRUE,Verbose=FALSE){
	nVertices=nrow(PointData)
	EndPoints=subset(PointData,NeighbourCount==1)
	BranchPoints=subset(PointData,NeighbourCount>2)
	
	# a single unbranched segment will have 2 end points
	# if a branch point has 3 neighbours then it should increase
	# the number of free ends by one.
	# if 4 neighbours then free ends += 2 etc
	
	PredictedEndPoints=sum(BranchPoints$NeighbourCount-2)+2
	if(nrow(EndPoints)!=PredictedEndPoints)
		warning("Mismatch between number of end points (",nrow(EndPoints),
				") and number predicted from branch statistics (",PredictedEndPoints,")")
	
	# Start off with the start point as the first endpoint
	# we may change our mind if we try to parse multiple trees
	if(length(Origin) && Origin%in%EndPoints$PointNo){
		StartPoint=Origin 
	} else {
		StartPoint=min(EndPoints$PointNo)
		warning("No Origin specified. Using: ",StartPoint)
	}
	
	if(ProcessAllTrees){
		SubTrees=ParseEdgeListForAllSubTrees(Neighbours,Origin=Origin,Silent=!Verbose)
		nTrees=length(SubTrees)
		nPointsParsed=length(unique(unlist(SubTrees)))
		if(nPointsParsed != nrow(PointData))
			stop("subtrees do not include all points in neuron!")
		
		# NB recurs = F to prevent us just ending up with a point list
		SegList=unlist(SubTrees,recurs=FALSE)
		if (nTrees>1) {
			SegSubTrees=rep(1:nTrees,sapply(SubTrees,length))
			# I suppose at this juncture I should choose the StartPoint with
			# the largest associated subtree - I will use complexity as the
			# measure
			SubTreeLengths=sapply(SubTrees,length)
			# Pick the first point of the first segment of the largest subtree
			LargestSubTree=order(SubTreeLengths,decreasing=TRUE)[1]
			StartPoint=SubTrees[[LargestSubTree]][[1]][1]		
			SegList=SubTrees[[LargestSubTree]]
			
			PointData$SubTree=rep(-1,nrow(PointData))
			for(i in 1:length(SubTrees)){
				PointData$SubTree[PointData$PointNo%in%unique(unlist(SubTrees[[i]]))]=i
			}
			# New (and more correct method for assigning parents)
			PointData$Parent=-1
			for (tree in SubTrees){
				for(s in tree){
					PointData$Parent[s[-1]]=s[-length(s)]
				}
			}
		}
	} else {
		# Just look for one tree starting from root
		SegList=ParseEdgeList(Neighbours,RootPoint=StartPoint)
	}
	
	# New (and more correct method for assigning parents)
	# Check, since we may already have done this for multi subtree neurons
	if(is.null(PointData$Parent)){
		PointData$Parent=-1
		for(s in SegList){
			PointData$Parent[s[-1]]=s[-length(s)]
		}		
	}
	
	if(is.null(PointData$Label)){
		# set up a default label if there was not label information
		PointData$Label=2
	}
	firstFields=c("PointNo","Label","X","Y","Z","W","Parent")
	remainingFields=setdiff(names(PointData),firstFields)
	PointData=PointData[,c(firstFields,remainingFields)]
	
	# Remove any Branch or End points that didn't make it into SegList
	PointsInSubTree=unique(unlist(SegList))
	EndPoints=subset(EndPoints,PointNo%in%PointsInSubTree)
	BranchPoints=subset(BranchPoints,PointNo%in%PointsInSubTree)
	
	CoreNeuron<-list(
			NumPoints=nrow(PointData),
			StartPoint=StartPoint, # NB I am assuming that this is always 1
			BranchPoints=BranchPoints$PointNo,
			EndPoints=EndPoints$PointNo,
			NumSegs=length(SegList),
			SegList=SegList,
			nTrees=nTrees,
			d=PointData)
	# If there are multiple subtrees then make that data available as well
	if(nTrees>1)
		CoreNeuron$SubTrees=SubTrees

	as.neuron(CoreNeuron)
}

ParseAM3DToNeuron=function(datalist,filename,
  method=c("original","igraph"),Force=FALSE,
  ProcessAllTrees=TRUE,Verbose=FALSE){
	# function to parse an amira mesh 3D
	# format skeleton tree produced by the Amira Skeletonize plugin
	# 050505 - This function has been rewritten extensively over the
	# previous week - and now works! 
	# Reading in the data and parsing the resultant edge list
	# have been delegated to a two separate functions
	# One major issue that came up with reading Amira files is that a
	# lot of them have multiple sub-trees because this is an easy mistake
	# to make in the Skeletonize plugin.  I have ended up reading the
	# largest tree into SegList and setting StartPoint to tha first node of
	# that subtree as a default - if ProcessAllTrees = FALSE
	# then only the first subtree will be read. 
	# 050512 In either case have decided that BranchPoints & EndPoints
	# should correspond only to the main subtree in SegList
	# Plotting functions will have to find the additional branch points
	# if needed
	
	# Bail out if we couldn't read any data
	if(is.null(datalist)) return(NULL)
  method=match.arg(method)
	Neighbours=datalist$EdgeList
	PointData=datalist$PointList
  if(method=='original'){
    CoreNeuron=CoreNeuronFromPointAndEdgeData(
      datalist$PointList,datalist$EdgeList,datalist$Origin,
      ProcessAllTrees = ProcessAllTrees, Verbose = Verbose)
  } else {
    CoreNeuron=CoreNeuronFromAmiraSkel(datalist, Verbose=Verbose)
  }
	
	ParsedNeuron<-c(list(NeuronName=NeuronNameFromFileName(filename),
			InputFileName=filename,
			CreatedAt=Sys.time(),
			NodeName=Sys.info()["nodename"],
			InputFileStat=file.info(filename)[1,],
			InputFileMD5=md5sum(path.expand(filename))),
			CoreNeuron)
	if(!is.null(datalist$SegmentProps)){
		# If we have materials information for the segments, keep that
		ParsedNeuron$SegmentProps=datalist$SegmentProps
	}
	return(as.neuron(ParsedNeuron))
}

ParseEdgeList<-function(Nb,Silent=TRUE,Verbose=!Silent,RootPoint=1){
	# Takes an edge list from a Skeletonize file
	# ie a listing of all edges bidirectionally (ie 2 edges for every pair 
	# of linked points) and produces a SegList
	
	# Make a ragged array lookup table for neighbours
	if(is.null(colnames(Nb))) colnames(Nb)=c("CurPoint","Neighbour")
	#Nb=as.data.frame(Nb)
	lNeighboursFromPoint=split(Nb[,"Neighbour"],Nb[,"CurPoint"])
	
	SegList=list()
	BranchPoints=unique(Nb[,"CurPoint"])[table(Nb[,"CurPoint"])>2]
	EndPoints=unique(Nb[,"CurPoint"])[table(Nb[,"CurPoint"])==1]	
	
	parseSegment<-function(HeadPoint){
		# Check if the head point has any neighbours
		cHeadPoint=as.character(HeadPoint)
		if(length(lNeighboursFromPoint[[cHeadPoint]])==0) {
			if(Verbose) cat("Bailing out because no neighbours for headpoint",HeadPoint,"\n")
			return(0)
		}
		# Make a segment starting with the head point
		SegList[[length(SegList)+1]]<<-HeadPoint
		if(Verbose) cat("Seg:",length(SegList),"Added HeadPoint",HeadPoint,"\n")
		PrevPoint=HeadPoint
		while(T){
			#str(lNeighboursFromPoint); str(SegList)
			# Find the CurPoint - by picking the prev point's 1st neighbour
			cPrevPoint=as.character(PrevPoint)
			PossCurPoints=lNeighboursFromPoint[[cPrevPoint]]
			if(Verbose) cat("Seg:",length(SegList),"PrevPoint",PrevPoint,"has",length(PossCurPoints),"Neighbours:",PossCurPoints,"\n")
			# If no CurPoint return 0
			if(length(PossCurPoints)<1){
				if(Verbose) cat("Seg:",length(SegList),"Bailing out since no more points\n")
				return (0)
			} else CurPoint=PossCurPoints[1]

			# else 
			# add the point to the current segment
			SegList[[length(SegList)]]<<-c(	SegList[[length(SegList)]],CurPoint[1])
			if(Verbose) cat("Seg:",length(SegList),"Added Point",CurPoint,"\n")
			if(Verbose) print(lNeighboursFromPoint)
			# Remove the relevant PrevPoint / CurPoint edges
			lNeighboursFromPoint[[cPrevPoint]]<<-
				lNeighboursFromPoint[[cPrevPoint]][!(lNeighboursFromPoint[[cPrevPoint]]==CurPoint)]
			cCurPoint=as.character(CurPoint)
			lNeighboursFromPoint[[cCurPoint]]<<-
				lNeighboursFromPoint[[cCurPoint]][!(lNeighboursFromPoint[[cCurPoint]]==PrevPoint)]
			# What kind of point is the CurPoint? 
			# 1: an end point - terminate
			if(any(EndPoints==CurPoint)) return (0)
			Children=lNeighboursFromPoint[[cCurPoint]]
			# 2: a branch - process children and terminate
			if(any(CurPoint==BranchPoints)){
				# if branch point (ie has multiple neighbours
				# then loop over the dependent segments
				if(Verbose) cat("Seg:",length(SegList),"Looping over",length(Children),"subtrees\n")
				for ( Child in Children){					
					# parse the dependent segment
					parseSegment(CurPoint)
				}
				return(0)
			}
			# 3: else simple point, so just loop
			PrevPoint=CurPoint
		}
		return(0)
	}
	
	parseSegment(RootPoint)
	
	UnusedPoints=setdiff(unique(Nb[,"CurPoint"]),unique(unlist(SegList)))
	if(!Silent && any(UnusedPoints)){
		cat("The following points were not used in the SegList:",UnusedPoints,"\n")
		warning("")
	}
	
	return(SegList)
}	

ParseEdgeListForAllSubTrees<-function(Nb,Origin=NULL,Silent=T){
	# Repeatedly call ParseEdgeList with different roots until all points
	# accounted for
	if(is.null(colnames(Nb))) colnames(Nb)=c("Neighbour","CurPoint")
	pointsRemaining=sort(unique(Nb[,'CurPoint']))
	EndPoints=pointsRemaining[table(Nb[,'CurPoint'])==1]
	SegLists=list()
	if(!is.null(Origin) && any(Origin==EndPoints)) StartPoint=Origin else StartPoint=min(EndPoints)
	cat("Origin=",Origin,"StartPoint=",StartPoint,"\n")
	while(any(pointsRemaining)){
		cat("SegLists",length(SegLists),"StartPoint=",StartPoint)
		SegLists[[length(SegLists)+1]]=ParseEdgeList(Nb,RootPoint=StartPoint,Silent=Silent)
		if(length(SegLists[[length(SegLists)]])==0) break
		pointsRemaining=setdiff(pointsRemaining,unique(unlist(SegLists[[length(SegLists)]])))
		if(length(pointsRemaining)==0) break
		StartPoint=min(intersect(EndPoints,pointsRemaining)	)
	}
	SegLists
}

read.neuron.amiraskel<-function(AM3DFile,...){
	datalist=ReadAM3DData(AM3DFile)
	MyNeuron<-ParseAM3DToNeuron(datalist,AM3DFile,...)
}

ReadNeuronFromAM3D<-function(AM3DFile,Components="Axon",OldNeuron=NULL,ReOrient=F,WithConvexHull=T,...){
	# Function to read in a neuron from an AmiraMesh 3D ascii file
	# this is the file format of the Skeletonize plugin
	# I have adapted this from ReadNeuronFromAsc
	# However not all the options are functional now since they are not
	# really appropriate to the new format
	# In particular there is no contour data associated with this file
	# format.
	# nb if WithConvHull ==T then use convex hull to reduce number of data 
	# points in LH and MB contours 

	# .Deprecated(read.neuron.amiraskel)
	
	if(any(Components=="Axon")){	
		# WE'RE GOING TO BUILD A NEW NEURON FROM SCRATCH
		# 2. Get the axon data from the file
		
		cat("\n",AM3DFile,":")
		datalist=ReadAM3DData(AM3DFile)
		MyNeuron<-ParseAM3DToNeuron(datalist,AM3DFile,...)
		if(is.null(MyNeuron$SegList)){
			stop("Unable to Extract neuron from",AM3DFile,"\n")
		}
		
		
		# No need for orientation data, since we now insist that everything 
		# comes in in the same orientation ie left hand side of brain (fly's
		# left) viewed from the anterior with dorsal up and therefore lateral
		# to the right.  However have left in the call for compatibility
		# NB ReOrient=FALSE because I don't have any AscData to parse anyway
		# Read in Orientation data from Neuron and figure out
		# which way up and which orientation it is in
		# Also add a set of scale information
		MyNeuron<-OrientNeuron(MyNeuron,NULL,AM3DFile,ReOrient=FALSE)
		
		# Get Marker Points eg Lateral Horn Entry Point
		if(!is.null(MyNeuron$OrientInfo)){
			# No equivalent for .am files - am I going to need this info
			# or should I be trying to calculate it on the basis of
			# proximity to lateral horn / 
			# MyNeuron<-GetMarkerPoints(MyNeuron,AscData,AM3DFile)    
		}
	} else {
	# WE'RE GOING TO USE A SUPPLIED NEURON AS A BASE 
		if(is.null(OldNeuron))
			stop("Must supply either an OldNeuron or read in AxonData in ReadNeuronFromAsc")
		if(is.null(OldNeuron$OrientInfo)) stop(
			paste("Neuron",OldNeuron$NeuronName,
				"is in an old format without $OrientInfo and must be entirely reloaded from its AmiraMesh3D file"))
		MyNeuron<-OldNeuron
	}
	if(any(Components=="LH")){
		warning(paste("Don't know how to add LH data to AmiraMesh 3D file",AM3DFile))
	}
	
	if(any(Components=="MB")){
		warning(paste("Don't know how to add LH data to AmiraMesh 3D file",AM3DFile))
	}
	
	return(MyNeuron)
}

WriteNeuronToAM<-function(ANeuron,AMFile=NULL,
	suffix="am",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE,
  WriteRadius=TRUE){
	# write out a neuron in the basic AmiraMesh format which is the native format
	# of amira for linesets (as opposed to the specialised skeletonize AM3D)
	# WriteAllSubTrees will write out all the stored subtrees in a neuron 
	# which has multiple subtrees (which is often true of ill-formed 
	# skeletonize neurons)
	
	if(is.null(AMFile))
		AMFile=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
	else if(isTRUE(file.info(AMFile)$isdir)){
		# we've been given a directory
		# we want to write a file into this directory with same name as original
		AMFile=file.path(AMFile,
			paste(sub("(.*)\\.[^.]*$","\\1",basename(ANeuron$InputFileName)),sep=".",suffix))
	}
	
	if(!Force && file.exists(AMFile) ){
		warning(paste(AMFile,"already exists; use Force=T to overwrite"))
		return()
	}
	if(!file.exists(dirname(AMFile))){
		# either bail
		if(!MakeDir){
			warning(paste(dirname(AMFile),"does not exist; use MakeDir=T to create"))
			return()
		} else {
			# or try to make a directory
			if(!dir.create(dirname(AMFile),recursive=TRUE) ){
				warning(paste("Unable to create",dirname(AMFile)))
			}
		}
	}
	if(!file.create(AMFile)){
		warning(paste("Unable to write to file",AMFile))
		return()
	}

	# if asked & nTrees is present and >1
	if(WriteAllSubTrees && !is.null(ANeuron$nTrees) && ANeuron$nTrees>1){	
		WriteAllSubTrees=TRUE 
		# nb recurs =F, so list of lists -> list (rather than vector)
		SegList=unlist(ANeuron$SubTrees,recurs=F)
	} else {
		WriteAllSubTrees=FALSE
		SegList=ANeuron$SegList
	}
	chosenVertices=sort(unique(unlist(SegList)))
	nVertices=length(chosenVertices)
	# the number of points required to define the line segments
	# including the terminating -1s (1 for each segment)
	nLinePoints=length(unlist(SegList))+length(SegList) 
	
	cat("Writing to",AMFile,"\n")
	# Write the header
	cat("# AmiraMesh ASCII 1.0\n",file=AMFile)
	fc=file(AMFile,open="at") # ie append, text mode

	cat("# Created by WriteNeuronToAM - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	
	cat("define Lines",nLinePoints,"\n",file=fc)
	cat("define Vertices", nVertices,"\n\n",file=fc)
	
	cat("Parameters {\n",file=fc)
	cat("    ContentType \"HxLineSet\"\n",file=fc)
	cat("}\n\n",file=fc)
  sectionNumbers=c(Coordinates=1,LineIdx=2)
	cat("Vertices { float[3] Coordinates } = @1\n",file=fc)
  if(WriteRadius){
    cat("Vertices { float Data } = @2\n",file=fc)
    sectionNumbers=c(Coordinates=1,Data=2,LineIdx=3)
  }
	cat("Lines { int LineIdx } = @",sectionNumbers['LineIdx'],"\n",sep="",file=fc)
	if(WriteAllSubTrees) {
    sectionNumbers=c(sectionNumbers,Data2=max(sectionNumbers)+1)
    cat("Vertices { float Data2 } =@",sectionNumbers['Data2'],"\n",sep="",file=fc)
  }
	cat("\n",file=fc)
	
	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(ANeuron$d[chosenVertices,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
  
	# Write the Radii
  if(WriteRadius){
      cat("\n@",sectionNumbers['Data']," # ",nVertices," width values\n",sep="",file=fc)
      # NB Divide width by 2
      write.table(ANeuron$d$W[chosenVertices]/2,col.names=F,row.names=F,file=fc,na='NaN')
  }
  
	# Write the segment information
	cat("\n@",sectionNumbers['LineIdx']," #",nLinePoints,"line segments\n",sep="",file=fc)
	# nb have to -1 from each point because amira is 0 indexed
	# AND add -1 to each segment as a terminator
	tmp=lapply(SegList,function(x) cat(x-1,"-1 \n",file=fc) )
	if(WriteAllSubTrees) {
    cat("\n@",sectionNumbers['Data2']," # subtrees\n",sep="",file=fc)
		if(ScaleSubTreeNumsTo1) ANeuron$d$SubTree=ANeuron$d$SubTree/max(ANeuron$d$SubTree)
		write.table(ANeuron$d$SubTree,col.names=F,row.names=F,file=fc)
	}
	close(fc)
}

WritePointsToAM<-function(d,AMFile=NULL,suffix="am",Force=F,MakeDir=T){
	# write out a set of points in the  basic AmiraMesh format

	if(is.neuron(d)){
		if(is.null(AMFile)) AMFile=paste(sub("(.*)\\.[^.]*$","\\1",d$InputFileName),sep=".",suffix)
		d=d$d
	} else if (is.null(AMFile)){
		warning("No file name specified in WritePointsToAM")
		return()
	}
	
	if(!Force && file.exists(AMFile) ){
		warning(paste(AMFile,"already exists; use Force=T to overwrite"))
		return()
	}
	if(!file.exists(dirname(AMFile))){
		# either bail
		if(!MakeDir){
			warning(paste(dirname(AMFile),"does not exist; use MakeDir=T to overwrite"))
			return()
		} else {
			# or try to make a directory
			if(!dir.create(dirname(AMFile))){
				warning(paste("Unable to create",dirname(AMFile)))
			}
		}
	}
	if(!file.create(AMFile)){
		warning(paste("Unable to write to file",AMFile))
		return()
	}
  
  # assign column names 
  if(ncol(d)==3){
    colnames(d)=c("X","Y","Z")
  } else if(is.null(colnames(d))){
		stop("Unable to identify X,Y,Z coordinates")
	} else {
      d=d[,c("X","Y","Z")]
	}
	
	cat("Writing to",AMFile,"\n")
	# Write the header
	cat("# AmiraMesh ASCII 1.0\n",file=AMFile)
	fc=file(AMFile,open="at") # ie append, text mode

	cat("# Created by WritePointsToAM - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	nVertices=nrow(d)
	cat("define Markers",nVertices,"Parameters {\nContentType \"LandmarkSet\",nSets 1\n}\n",file=fc)
# 	cat("Parameters {\n",file=fc)
# 	cat("    ContentType \"HxLineSet\"\n",file=fc)
# 	cat("}\n\n",file=fc)

	cat("Markers { float[3] Coordinates } = @1\n",file=fc)
	cat("\n",file=fc)
	
	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(d[,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
	close(fc)
}

WriteNeuronToAM3D<-function(ANeuron,AMFile=NULL,
	suffix="am3",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE,sep=NULL){
	# write out a neuron in the specialised skeletonize AM3D format 
	# (as opposed to the basic AmiraMesh format which is the native format
	# of amira for linesets)
	# WriteAllSubTrees will write out all the stores subtrees in a neuron 
	# which has multiple subtrees (which is often true of ill-formed 
	# skeletonize neurons).  It will also add a data field that can be used
	# to visualised different subtrees eg by colouring
	
	if(is.null(AMFile))
		AMFile=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
	else if(isTRUE(file.info(AMFile)$isdir)){
		# we've been given a directory
		# we want to write a file into this directory with same name as original
		AMFile=file.path(AMFile,
			paste(sub("(.*)\\.[^.]*$","\\1",basename(ANeuron$InputFileName)),sep=".",suffix))
	}
	
	if(!Force && file.exists(AMFile) ){
		warning(paste(AMFile,"already exists; use Force=T to overwrite"))
		return()
	}
	if(!file.exists(dirname(AMFile))){
		# either bail
		if(!MakeDir){
			warning(paste(dirname(AMFile),"does not exist; use MakeDir=T to overwrite"))
			return()
		} else {
			# or try to make a directory
			if(!dir.create(dirname(AMFile))){
				warning(paste("Unable to create",dirname(AMFile)))
			}
		}
	}
	if(!file.create(AMFile)){
		warning(paste("Unable to write to file",AMFile))
		return()
	}

	# if asked & nTrees is >1  (NB isTRUE handles NULL case correctly)
	if(WriteAllSubTrees && isTRUE(ANeuron$nTrees>1)){	
		WriteAllSubTrees=TRUE 
		# nb recurs =F, so list of lists -> list (rather than vector)
		SegList=unlist(ANeuron$SubTrees,recurs=F)
	} else {
		WriteAllSubTrees=FALSE
		SegList=ANeuron$SegList
	}
	chosenVertices=sort(unique(unlist(SegList)))
	nVertices=length(chosenVertices)
	cat("nVertices =",nVertices,"\n")
	# I think that restristicting to chosen vertices without
	# any renumbering of points is problematic
	# FIXME - need to renumber vertices if there are NAs
#	nVertices=nrow(ANeuron$d)
	# the number of points required to define the edge list
	# if each segment contains n points, then 2(n-1) edges
	nEdgeList=sum(sapply(SegList,length)-1)*2
	# Make EdgeList
	makeEdges=function(seg){
			lSeg=length(seg)
			if(lSeg<2) return()
			elFwd=cbind(seg[-lSeg],seg[-1])
			elRev=cbind(seg[-1],seg[-lSeg])
			df=data.frame(rbind(elFwd,elRev))
			names(df)=c("CurPoint","Neighbour")
			df[order(df$CurPoint,df$Neighbour),]
	}
	EdgeList=do.call('rbind',lapply(SegList,makeEdges))
	EdgeList=EdgeList[order(EdgeList$CurPoint,EdgeList$Neighbour),]
	
	cat("Writing to",AMFile,"\n")
	# Write the header
	cat("# AmiraMesh 3D ASCII 2.0\n",file=AMFile)
	fc=file(AMFile,open="at") # ie append, text mode

	cat("# Created by WriteNeuronToAM3D -",format(Sys.time(),usetz=T),"\n\n",file=fc)	
	cat("nVertices", nVertices,"\nnEdges",nEdgeList,"\n",file=fc)

	vertexTypeList=ifelse(WriteAllSubTrees,nVertices,0) 
	cat("define Origins 1\ndefine vertexTypeList",vertexTypeList,"\n\n",file=fc)
	
	cat("Parameters {\n",file=fc)
	cat("    ContentType \"SkeletonGraph\"\n",file=fc)
	cat("}\n\n",file=fc)

	cat("Vertices { float[3] Coordinates } @1\n",file=fc)
	cat("Vertices { int NeighbourCount } @2\n",file=fc)
	cat("Vertices { float Radii } @3\n",file=fc)
	cat("EdgeData { int NeighbourList } @4\n",file=fc)
	cat("Origins { int Origins } @5\n",file=fc)
	cat("Vertices { int vertexTypeCounter } @6\n",file=fc)
	cat("vertexTypeList { int vertexTypeList } @7\n\n",file=fc)

	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	if(is.null(sep)){
		# Amira seems fussy about having nicely aligned columns
		# using format with trim = FALSE (the default actually) 
		# and after getting rid of names results in a nicely justified table
		Coords=as.matrix(ANeuron$d[,c("X","Y","Z")])
		rownames(Coords)<-colnames(Coords)<-NULL
		write.table(format(Coords,trim=FALSE,scientific=FALSE),
			quote=F,row.names=FALSE,col.names=FALSE,file=fc)
	} else {
		# sep was explicitly specified, so use that
		write.table(ANeuron$d[chosenVertices,c("X","Y","Z")],col.names=F,row.names=F,file=fc,sep=sep)
	}
	
	# Write number of neighbours
	cat("\n@2 #",nVertices,"numbers of neighbours \n",file=fc)
	numNeighbours=integer(nVertices) # filled with 0s
	numNeighbours[sort(unique(EdgeList$CurPoint))]=table(EdgeList$CurPoint)
	write.table(numNeighbours,col.names=F,row.names=F,file=fc)
	
	# Write the Radii
	cat("\n@3 #",nVertices,"radii\n",file=fc)
	# NB Divide width by 2
	write.table(ANeuron$d$W[chosenVertices]/2,col.names=F,row.names=F,file=fc)

	# Write the edgelist information
	cat("\n@4 #",nEdgeList,"bidirectional edges\n",file=fc)
	#NB -1 since Amira is 0 indexed
	write.table(EdgeList$Neighbour-1,col.names=F,row.names=F,file=fc)

	# Write the origin information NB -1 since 0 indexed
	cat("\n@5 #n 1\n",file=fc)
	cat(ANeuron$StartPoint-1,"\n",file=fc)
	
	# Write the vertexTypeCounter information
	cat("\n@6 #",nVertices,"\n",file=fc)
	cat(paste(rep(0,nVertices),"\n"),sep="",file=fc)
	
	# nb have to -1 from each point because amira is 0 indexed
	# AND add -1 to each segment as a terminator
	#tmp=do.call("paste",ANeuron$SegList)
	#writeLines(tmp,con=fc)
	if(WriteAllSubTrees) {
		cat("\n@7 # subtrees\n",file=fc)
		if(ScaleSubTreeNumsTo1) ANeuron$d$SubTree=ANeuron$d$SubTree/max(ANeuron$d$SubTree)
		write.table(ANeuron$d$SubTree,col.names=F,row.names=F,file=fc)
	}
	cat("\n",file=fc)
	close(fc)
}


ReadMBLHFromAMSurf<-function(AMSurfFile,Components=c("MB","LH"),WithConvexHull=F,...){
	# Function to read in a neuron from an AmiraMesh 3D ascii file
	# this is the file format of the Skeletonize plugin
	# I have adapted this from ReadNeuronFromAsc
	# However not all the options are functional now since they are not
	# really appropriate to the new format
	# In particular there is no contour data associated with this file
	# format.
	# nb if WithConvHull ==T then use convex hull to reduce number of data 
	# points in LH and MB contours 

	if(length(intersect(Components,c("MB","LH")))>0 ){	
		# We have some data to get
		cat("We have some data to get\n")
		Data<-ParseAMSurfToContourList(AMSurfFile,...)
		return(Data) 
	}
	return(-1)
}


ParseAMSurfToContourList<-function(filename,RegionNames="ALL",RegionChoice="Inner",Verbose=FALSE,FallbackRegionCol="grey"){
	.Deprecated('read.hxsurf','nat')
	nat::read.hxsurf(filename,RegionNames=if(RegionNames=="ALL") NULL else RegionNames,
		RegionChoice=RegionChoice,Verbose=Verbose,FallbackRegionCol=FallbackRegionCol)
}

WriteHxSurface=function(filename,Vertices,Indices=NULL,
		material=sub("^([^.]+)\\..*","\\1",basename(filename))){
		
	cat("# HyperSurface ASCII\nParameters {\n",file=filename)
	fc=file(filename,open="at") # ie append, text mode
	#cat("\tMaterials{\n\t\tInterior {\n\t\t\tid 0\n\t\t}\n",file=fc)
	cat("\t{color 0.83562 0.78 0.06,\nName \"",sep="",material,"\"}",file=fc)
	cat("\t\t",sep="",material," {\n\t\t\tid 1\n\t\t}\n",file=fc)
	cat("\t}\n",file=fc)
	cat("\tBoundaryIds {\n\t\tId0 {\n\t\t\tId 0,Info \"undefined\",Color 0.6 0.6 0.6\n\t\t}\n\t\tname \"BoundaryConditions\"\n\t}\n}\n",file=fc)


	nVertices=nrow(Vertices)
	cat("Vertices",nVertices,"\n",file=fc)
	write.table(Vertices,col.names=F,row.names=F,file=fc)
	cat("Patches 1\n{\tInnerRegion Interior \n OuterRegion", material,"\n BoundaryID 0\n BranchingPoints 0 \n\n Triangles",round(nVertices/3),"\n",file=fc)
	#cat("Markers { float[3] Coordinates } = @1\n@1\n",file=filename,append=T)
	Indices=matrix(1:nVertices,nrow=3)
	
	cat(apply(t(Indices),1,paste,collapse=" "),sep="\n",file=fc)
	cat("}\n",file=fc)
}

Write3DDensityToAmiraRectilinear<-function(filename,d){
		# Produces a Rectilinear format file -
		# that is one with an arbitrary x,y,z grid
		cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
		fc=file(filename,open="at") # ie append, text mode
		lattice=apply(d$eval.points,2,length)
		cat("define Lattice",lattice,"\n",file=fc)
		cat("define Coordinates",sum(lattice),"\n\n",file=fc)
		cat("Parameters {CoordType \"rectilinear\"}\n\n",file=fc)
		cat("Lattice { float ScalarField } = @1\n",file=fc)
		cat("Coordinates { float xyz } = @2\n\n",file=fc)
		cat("@1\n",file=fc)
		
		cat(as.vector(d$estimate),file=fc)
		cat("@2\n",file=fc)
		cat(d$eval.points[,1],"\n",file=fc)
		cat(d$eval.points[,2],"\n",file=fc)
		cat(d$eval.points[,3],"\n",file=fc)
}

Write3DDensityToAmiraLattice<-function(filename,dens,ftype=c("binary","text","hxzip"),
	dtype=c("float","byte", "short", "ushort", "int", "double"),WriteNrrdHeader=FALSE,endian=c('big','little')){
	# Produces a lattice format file -
	# that is one with a regular x,y,z grid
	# Can also write a detached Nrrd header that points to the AmiraMesh
	# data to allow it to be opened by a nrrd reader
	ftype=match.arg(ftype)
	endian=match.arg(endian)
	if(ftype=='text') cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
	else if(endian=='little') cat("# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n",file=filename)
	else cat("# AmiraMesh 3D BINARY 2.0\n\n",file=filename)
	
	fc=file(filename,open="at") # ie append, text mode
	cat("# Created by Write3DDensityToAmiraLattice - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	if(!is.list(dens)) d=dens else d=dens$estimate
	# Find data type and size for amira
	dtype=match.arg(dtype)	
	dtypesize<-c(4,1,2,2,4,8)[which(dtype==c("float","byte", "short","ushort", "int", "double"))]
	# Set the data mode which will be used in the as.vector call at the
	# moment that the binary data is written out.
	if(dtype%in%c("byte","short","ushort","int")) dmode="integer"
	if(dtype%in%c("float","double")) dmode="numeric"
	
	
	#lattice=apply(d$eval.points,2,length)
	lattice=dim(d)
	cat("define Lattice",lattice,"\n",file=fc)

	cat("Parameters { CoordType \"uniform\",\n",file=fc)
	# note Amira's definition for the bounding box:
	# the range of the voxel centres.
	# So eval.points should correspond to the CENTRE of the
	# voxels at which the density is evaluated
	cat("\t# BoundingBox is xmin xmax ymin ymax zmin zmax\n",file=fc)
	BoundingBox=NULL
	if(!is.null(attr(dens,"BoundingBox"))){
		BoundingBox=attr(dens,"BoundingBox")
	} else if(is.list(d) && !is.null(d$eval.points)){
		BoundingBox=as.vector(apply(d$eval.points,2,range))
	}
	if(!is.null(BoundingBox)) cat("\t BoundingBox",BoundingBox,"\n",file=fc)
	cat("}\n\n",file=fc)
	
	if(ftype=="hxzip"){
		raw_data=writeBin(as.vector(d,mode=dmode),raw(),size=dtypesize,endian=endian)
		zlibdata=nat:::write.zlib(raw_data)
		cat("Lattice { ",dtype," ScalarField } = @1(HxZip,",length(zlibdata),")\n\n",sep="",file=fc)
	} else cat("Lattice {",dtype,"ScalarField } = @1\n\n",file=fc)

	cat("@1\n",file=fc)
	
	#cat(str(as.vector(d)))
	close(fc)

	# Write a Nrrd header to accompany the amira file if desired
	# see http://teem.sourceforge.net/nrrd/
	if(WriteNrrdHeader) {
		if(ftype=="hxzip") stop("Nrrd cannot cope with Amira's HxZip encoding (which is subtly different from gzip)")
		nrrdFilename=paste(filename,sep=".","nhdr")
		cat("NRRD0004\n",file=nrrdFilename)
		fc=file(nrrdFilename,open="at") # ie append, text mode
		nrrdType=ifelse(dtype=="byte","uint8",dtype)
		
		cat("encoding:", ifelse(ftype=="text","text","raw"),"\n",file=fc)
		cat("type: ",nrrdType,"\n",sep="",file=fc)
		cat("endian: ",endian,"\n",sep="",file=fc)
		# Important - this sets the offset in the amiramesh file from which
		# to start reading data
		cat("byte skip:",file.info(filename)$size,"\n",file=fc)
		cat("dimension: ",length(lattice),"\n",sep="",file=fc)
		cat("sizes:",lattice,"\n",file=fc)
		voxdims=voxdim.gjdens(dens)
		if(!is.null(voxdims)) cat("spacings:",voxdims,"\n",file=fc)
		if(!is.null(BoundingBox)){
			cat("axis mins:",matrix(BoundingBox,nrow=2)[1,],"\n",file=fc)
			cat("axis maxs:",matrix(BoundingBox,nrow=2)[2,],"\n",file=fc)
		}
		cat("data file: ",basename(filename),"\n",sep="",file=fc)
		cat("\n",file=fc)
		close(fc)
	}

	if(ftype=='text'){
		write(as.vector(d,mode=dmode),ncol=1,file=filename,append=TRUE)
	} else {
		fc=file(filename,open="ab") # ie append, bin mode
		if(ftype=="hxzip")
			writeBin(zlibdata,fc,size=1,endian=endian)
		else
			writeBin(as.vector(d,mode=dmode),fc,size=dtypesize,endian=endian)
		close(fc)
	}
}

Read3DDensityFromAmiraLattice<-function(filename,Verbose=FALSE){

	fc=file(filename,'rb')
	headerLines<-readLines(fc,1)
	if(length(grep("amiramesh",headerLines,ignore.case=T))!=1)
		stop("This is not an amiramesh file")
	while ( ( nextLine<-readLines(fc,1)) !="@1") {headerLines<-c(headerLines,nextLine)}
	
	endian='big'
	binary=FALSE
	# Figure out if the file is in binary format or not
	if(length(grep("LITTLE.ENDIAN",headerLines[1],ignore.case=TRUE))>0){
#		any(grep("^\\s*#\\s+amiramesh(\\s+3d)?\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
		binary = TRUE
		endian = 'little'
	} else if(any(grep("^\\s*#\\s+amiramesh(\\s+3d)?\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
		binary = TRUE
	}

	# Find the position of the header lines defining the lattice
	# this clearly assumes that the relevant data is in position
	# @1 in the file - which may or may not be the case.
	LatticeDefLine=grep("define\\s+lattice",headerLines,ignore.case=TRUE,perl=TRUE)
	LatticeTypeDefLine=grep("^Lattice.*}\\s*[=]{0,1}\\s*@1",headerLines,ignore.case=TRUE,perl=FALSE)
	#cat("LatticeTypeDefLine = ",LatticeTypeDefLine)
	LatticeBoundsLine=grep("^[^#]*BoundingBox",headerLines,ignore.case=TRUE,perl=TRUE)
	if(!any(LatticeDefLine)){
		warning(paste("No lattice definition line in file",filename))
		close(fc); return(-1)
	}
	
	if(!any(LatticeTypeDefLine)) {
		warning(paste("No lattice type definition line in file",filename))
		close(fc); return(-1)
	}

	strtrim<-function(str){
			str<-sub('\\s+$', '', str, perl = TRUE) ## Perl-style white space
			sub('^\\s+', '', str, perl = TRUE) ## Perl-style white space
	}
	
	# fetch the lattice dimensions
	latticeDims=sub(".*lattice([0-9 \\t]+).*","\\1",headerLines[LatticeDefLine],ignore.case=T)
	latticeDims=as.numeric(unlist(strsplit(strtrim(latticeDims),"\\s+")))
	dataLength=prod(latticeDims)
	
	# and the bounds
	latticeBounds=as.numeric(unlist(strsplit(strtrim(headerLines[LatticeBoundsLine]),"(\\s|,)"))[-1])		
	
	
	# and finally the data type
	ldl=headerLines[LatticeTypeDefLine]
	# dataTypeName=sub(".*\\{\\s*(\\w+)\\s+.?*\\s*\\}.*","\\1",ldl,ignore.case=T,perl=T)
	dataTypeName=sub(".*\\{\\s*(\\w+)\\s+.*?\\s*\\}.*","\\1",ldl,ignore.case=T,perl=T)
	# check if the data is encoded in some way
	#"Lattice { byte Labels } @1(HxByteRLE,391697)"
	postAt=strsplit(ldl,"@1")[[1]][2]
	dataEncoding=""
	if(isTRUE(grep("\\(Hx.*?\\)",postAt)==1)){
		postAt=sub(".*?\\(([^)]+)\\).*?","\\1",postAt)
		dataEncoding=toupper(strsplit(postAt,",")[[1]][1])
		compressedLength=as.integer(strsplit(postAt,",")[[1]][2])
	}
	
#     Amira docs: The primitive data types must be
# 		one of byte, short, int, float, double, or complex.  Vectors of
# 		primitive data types are allowed, aggregate structs are not, however.
	
# 		primType Returns the primitive data type of the field, i.e., the way
# 		how the values are represented internally.  A number with the following
# 		meaning is returned: 0 = bytes, 1 = 16-bit signed integers, 2 = 32-bit
# 		signed integers, 3 = 32-bit floating point values, 4 = 64-bit floating
# 		point values, 7 = 16-bit unsigned integers.
	
	# now read the data
	# note that  bytes are assumed to be unsigned
	# shorts could be either but will assume signed - don't know how
	# amiramesh specifies either way
	dataTypes=data.frame(name=I(c("byte", "ushort", "short", "int", "float", "double", "complex")),
			size=c(1,2,2,4,4,8,NA),what=I(c(rep("integer",4),rep("numeric",2),NA)),
			signed=rep(c(FALSE,TRUE),c(2,5)) )
	i=which(dataTypes$name==dataTypeName)
	if(!any(i==1:7)){
		close(fc)
		stop("Unrecognised data type")
	}
	
	if(Verbose) cat("dataLength =",dataLength,"dataType =",dataTypes$what[i],"size=",dataTypes$size[i],"\n")
	if(binary){
		if(dataEncoding=="HXBYTERLE"){
			d=readBin(fc,what=raw(0),n=compressedLength,size=1)
			d=nat:::decode.rle(d,dataLength)
			d=as.integer(d)
		} else if(dataEncoding == "HXZIP"){
		  uncompressed=nat:::read.zlib(fc, compressedLength=compressedLength)
		  d=readBin(uncompressed, n=dataLength, what=dataTypes$what[i],
		            size=dataTypes$size[i], signed=dataTypes$signed[i],
                endian=endian)
		} else if(dataEncoding==""){
			d=readBin(fc,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
				signed=dataTypes$signed[i],endian=endian)
		} else {
			stop("Unimplemented data encoding",dataEncoding,"in file",filename,"\n")
		}
		close(fc)
	} else {
		# this clearly assumes that the relevant data is in position
		# @1 in the file - which may or may not be the case.
		# cat("dataTypes$what[i]=",dataTypes$what[i],"\n")
		if(dataTypes$what[i]=='integer') whatVal=integer(0) else whatVal=double(0)
		d=scan(fc,what=whatVal,nmax=dataLength)
		close(fc)
	}
	dim(d)<-latticeDims
	if(length(latticeBounds)>0){
		attr(d,"BoundingBox")<-latticeBounds
		attr(d,"x")<-seq(latticeBounds[1],latticeBounds[2],len=latticeDims[1])
		attr(d,"y")<-seq(latticeBounds[3],latticeBounds[4],len=latticeDims[2])
		attr(d,"z")<-seq(latticeBounds[5],latticeBounds[6],len=latticeDims[3])
	} else {
		# No Bounding Box available
# 			attr(d,"BoundingBox")<-NULL
	}
	return(d)
}

ReadAmiraLandmarks<-function(filename,Verbose=FALSE,CoordinatesOnly=TRUE){
  .Deprecated("read.amiralandmarks","nat")
  nat::read.amiralandmarks(filename, CoordinatesOnly=CoordinatesOnly, 
    Verbose=Verbose)
}

WriteAmiraLandmarks<-function(filename,d){
  .Deprecated("write.amiralandmarks","nat")
  nat::write.amiralandmarks(file=filename, x=d)
}

WriteAmiraColormap<-function(filename,rgba,A=1,minmax=c(0,255)){
	if(is.vector(rgba)) rgba=t(col2rgb(rgba)/255)
	
	if(ncol(rgba)==3) {
		rgba=cbind(rgba,A)
	}
	if(ncol(rgba)!=4) stop("Colormap must have 4 columns")
	if(nrow(rgba)!=256) warning("Colormap should have 256 levels to be editable with colormap editor")
	
	cmaprange=range(rgba)
	if(cmaprange[1]<0 || cmaprange[2]>1) stop ("Colormap values must be between 0 and 1")

	cat("# AmiraMesh ASCII 1.0\n\ndefine Lattice",nrow(rgba),
	"\n\nParameters {\n\tContentType \"Colormap\",\n\tMinMax",minmax,"\n}\n",file=filename)
	cat("Lattice { float[4] Data } = @1\n",file=filename,append=T)

	cat("@1\n",file=filename,append=T)
	write.table(rgba,col.names=F,row.names=F,file=filename,append=TRUE)
	cat("\n",file=filename,append=T)
}

WriteGenericAmiramesh<-function(filename,d,ContentType){
	dataDef=attr(d,"dataDef")
	if(is.null(dataDef)) stop("Cannot write without data definition")
	
	cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
	# definitions of lengths
	dimdefs=dataDef$Dims[unique(names(dataDef$Dims))]
	cat(paste("define",names(dimdefs),dimdefs,collapse="\n"),file=filename,append=TRUE)

	# Parameters - for now just handle content type
	if(missing(ContentType)) ContentType=attr(d,'Parameters')$ContentType
	if(!is.null(ContentType)){
		cat("\nParameters {\n\tContentType \"",, "\"\n}\n\n",file=filename,append=T)
	}

	with(dataDef,
		cat(paste(names(Dims)," { ",Type," ",DataName," } @",sep="",seq(nrow(dataDef)),collapse="\n"),
			file=filename,append=TRUE))
	
	cat("\n\n",file=filename,append=TRUE)
	if(!is.list(d)) d=list(d)
	for(i in seq(length(d))){
		if(length(d[[i]])==0) next 
		cat("@",i,"\n",sep="",file=filename,append=TRUE)
		write.table(d[[i]],row.names=FALSE,col.names=FALSE,sep="\t",file=filename,append=TRUE)
		cat("\n",file=filename,append=TRUE)
	}		
}


#' Read neuron in Amira's native lineset format
#' @param amfile Path to the amiramesh file
#' @param defaultDiameter If diameter information, missing use this default
#' @return A neuron object
#' @author jefferis
#' @export
#' @seealso \code{\link{read.neuron},\link{ReadNeuronFromAM3D}
ReadNeuronFromAM<-function(amfile,defaultDiameter=NA_real_){
  .Deprecated('read.neuron.hxlineset','nat')
  nat:::read.neuron.hxlineset(amfile, defaultDiameter=defaultDiameter)
}

