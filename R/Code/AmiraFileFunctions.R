# AmiraFileFunctions.R

# 2005-02-03
# Functions to parse AmiraMesh 3D format - the native
# ouput of the skeletonize plugin and to read and write the density
# data in Amira file formats.
# At the moment depends on SWCFunctions.s since the amiramesh 
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

# require(MASS)  # for write.matrix - call removed
require(tools) # for md5sum
#source(file.path(CodeDir,"SWCFunctions.s"))
source(file.path(CodeDir,"GraphTheory.R"))

ReadAmiramesh<-function(filename,DataSectionsToRead=NULL,Verbose=TRUE){
	# attempt to write a generic amiramesh reader
	firstLine=readLines(filename,n=1)
	if(!any(grep("#\\s+amiramesh",firstLine,ignore.case=T))){
		warning(paste(filename,"does not appear to be an AmiraMesh file"))
		return(NULL)
	}
	binaryfile="binary"==tolower(sub(".*(ascii|binary).*","\\1",firstLine,ignore.case=TRUE))

	con=if(binaryfile) file(filename,open='rb') else file(filename,open='rt')
	parsedHeader=ReadAmiramesh.Header(con,Verbose=Verbose)
	
	if(is.null(DataSectionsToRead)) DataSectionsToRead=parsedHeader$DataName
	else DataSectionsToRead=intersect(parsedHeader$DataName,DataSectionsToRead)	
	if(binaryfile){
		filedata=.ReadAmiramesh.BinaryData(con,parsedHeader,DataSectionsToRead,Verbose=Verbose)
		close(con)
	} else {
		close(con)
		#print(parsedHeader)
		filedata=.ReadAmiramesh.ASCIIData(filename,parsedHeader,DataSectionsToRead,Verbose=Verbose)
		#cat(length(filedata))
	}
	filedata
}

.ReadAmiramesh.BinaryData<-function(con,df,DataSectionsToRead,Verbose=TRUE){
	#print(df)
	l=list()
	# TODO
	for(i in seq(len=nrow(df))){
		if(Verbose) cat("Current offset is",seek(con),";",df$nBytes[i],"to read\n")
		
		if(all(DataSectionsToRead!=df$DataName[i])){
			# Just skip this section
			if(Verbose) cat("Skipping data section",df$DataName[i],"\n")
			seek(con,df$nBytes[i],origin="current")
		} else {
			if(Verbose) cat("Reading data section",df$DataName[i],"\n")
			if(df$RType[i]=="integer") whatval=integer(0) else whatval=numeric(0)
			x=readBin(con,df$SimpleDataLength[i],size=df$Size[i],what=whatval,signed=df$Signed[i])
			# note that first dim is moving fastest
			dims=unlist(df$Dims[i])
			# if the individual elements have subelements
			# then put those as innermost (fastest) dim
			if(df$SubLength[i]>1) dims=c(df$SubLength[i],dims)
			ndims=length(dims)
			if(ndims>1) dim(x)=dims
			if(ndims==2) x=t(x) # this feels like a hack, but ...
			l[[df$DataName[i]]]=x
		}		
		readLines(con,n=1) # Skip return at end of section
		nextSectionHeader=readLines(con,n=1)
		if(Verbose) cat("nextSectionHeader = ",nextSectionHeader,"\n")
	}
	l
}

.ReadAmiramesh.ASCIIData<-function(filename,df,DataSectionsToRead,Verbose=TRUE){
	l=list()
	df=subset(df,DataName%in%DataSectionsToRead)
	df=df[order(df$DataPos),]

	for(i in seq(len=nrow(df))){
		if(df$DataLength[i]>0){
			if(Verbose) cat("Reading ",df$DataLength[i],"lines from line",df$LineOffsets[i],"in file",filename,"\n")
			
			if(df$RType[i]=="integer") whatval=integer(0) else whatval=numeric(0) 
			l[[df$DataName[i]]]=scan(filename,what=whatval,n=df$SimpleDataLength[i],skip=df$LineOffsets[i],quiet=!Verbose)
			if(df$SubLength[i]>1){
				l[[df$DataName[i]]]=matrix(l[[df$DataName[i]]],ncol=df$SubLength[i],byrow=TRUE)
			}
		} else {
			if(Verbose) cat("Skipping data section",df$DataName[i],"\n")
		}
	}
	return(l)		
}


ReadAmiramesh.Header<-function(con,Verbose=TRUE){
	headerLines=NULL
	while( (t<-readLines(con,1))!="@1"){
		headerLines=c(headerLines,t)
	}
	nHeaderLines=length(headerLines)
	# trim comments and blanks & convert all white space to single spaces
	headerLines=trim(sub("(.*)\#.*","\\1",headerLines))
	headerLines=headerLines[headerLines!=""]
	headerLines=gsub("[[:space:]]+"," ",headerLines)
	
	#print(headerLines)
	# parse location definitions
	LocationLines=grep("^(n|define )(\\w+) ([0-9 ]+)$",headerLines,perl=TRUE)
	Locations=headerLines[LocationLines];headerLines[-LocationLines]
	LocationList=strsplit(gsub("^(n|define )(\\w+) ([0-9 ]+)$","\\2 \\3",Locations,perl=TRUE)," ") 
	LocationNames=sapply(LocationList,"[",1)
	Locations=lapply(LocationList,function(x) as.numeric(unlist(x[-1])))
	names(Locations)=LocationNames
	
	# parse parameters
	
	# parse data definitions
	DataDefLines=grep("^(\\w+).*@(\\d+)$",headerLines,perl=TRUE)
	DataDefs=headerLines[DataDefLines];headerLines[-DataDefLines]
	
	# remove all extraneous chars altogether
	DataDefs=gsub("(@|\\}|\\{|[[:space:]])+"," ",DataDefs)
	# make a df with DataDef info
	DataDefMatrix=matrix(unlist(strsplit(DataDefs," ")),ncol=4,byrow=T)
	
	DataDefDF=data.frame(DataName=I(DataDefMatrix[,3]),DataPos=as.numeric(DataDefMatrix[,4]))

	DataDefMatrix[,1]=sub("^EdgeData$","Edges",DataDefMatrix[,1])
	# Dims will store a list of dimensions that can be used later
	DataDefDF$Dims=Locations[DataDefMatrix[,1]] 
	DataDefDF$DataLength=sapply(DataDefMatrix[,1],function(x) prod(Locations[[x]])) #  notice prod in case we have multi dim
	DataDefDF$Type=I(DataDefMatrix[,2])
	DataDefDF$SimpleType=sub("(\\w+)\\s*\\[\\d+\\]","\\1",DataDefDF$Type,perl=TRUE)
	DataDefDF$SubLength=as.numeric(sub("\\w+\\s*(\\[(\\d+)\\])?","\\2",DataDefDF$Type,perl=TRUE))
	DataDefDF$SubLength[is.na(DataDefDF$SubLength)]=1

	# Find size of binary data (if required?)
	TypeInfo=data.frame(SimpleType=I(c("float","byte", "short", "int", "double", "complex")),Size=c(4,1,2,4,8,8),
		RType=I(c("numeric",rep("integer",3),rep("numeric",2))), Signed=c(TRUE,FALSE,rep(TRUE,4)) )
	DataDefDF=merge(DataDefDF,TypeInfo,all.x=T)
	# Sort (just in case)
	DataDefDF= DataDefDF[order(DataDefDF$DataPos),]
	
	DataDefDF$SimpleDataLength=DataDefDF$DataLength*DataDefDF$SubLength
	DataDefDF$nBytes=DataDefDF$SubLength*DataDefDF$Size*DataDefDF$DataLength
	DataDefDF$LineOffsets=nHeaderLines+c(0,cumsum(DataDefDF$DataLength+2)[-nrow(DataDefDF)])+1
	
	return(DataDefDF)
}

	

	
	
ReadAM3DData.Binary<-function(filename,OmitNAs=TRUE){
	# function to read in the basic data from the 
	# files produced by the Amira Skeletonize plugin
 	# Check for header confirming file type

	headerLines=NULL
	con=file(filename,open='rb')
	while( (t<-readLines(con,1))!="@1"){
		headerLines=c(headerLines,t)
	}

	# first read the Coordinates
	rval$Coordinates=matrix(readBin(con,what='numeric',n=nVertices*3,size=4,endian="big"),ncol=3,byrow=T)
	# next read the Neighbour Count
	t=readLines(con,2);cat(t[2],"\n")
	rval$NeighbourCount=readBin(con,what="integer",size=4,n=nVertices,endian="big")
	t=readLines(con,2);cat(t[2],"\n")
	rval$Radii=readBin(con,what='numeric',n=nVertices,size=4,endian="big")
	rval$NeighbourList=readBin(con,what="integer",size=4,n=nEdges,endian="big")
	rval$Origin=NULL
	close(con)
	return(rval)
}


ReadAM3DData<-function(filename,OmitNAs=TRUE){
	# function to read in the basic data from the 
	# files produced by the Amira Skeletonize plugin
	
	# 	# Check for header confirming file type
	firstLine=readLines(filename,n=1)
	if(!any(grep("#\\s+amiramesh\\s+3d",firstLine,ignore.case=T))){
		warning(paste(filename,"does not appear to be an AmiraMesh 3D file"))
		return(NULL)
	}
	filetype=ifelse(any(grep("binary",firstLine,ignore.case=T)),"binary","ascii")
		
	# Read Header
	headerLines=NULL
	con=file(filename,open='rb')
	#	while( (thisLine<-readLines(con,1))!="@1"){
	while( !isTRUE(charmatch("@1",thisLine<-readLines(con,1))==1) ){
		headerLines=c(headerLines,thisLine)
	}

	getfield=function(fName,pos=2) unlist(strsplit(headerLines[grep(fName,headerLines)],"\\s+",perl=TRUE))[pos]

	nVertices=as.numeric(getfield("nVertices",2))
	nEdges=as.numeric(getfield("nEdges",2))

	CoordinatesDefLine=grep("Coordinates",headerLines)
	NeighbourCountDefLine=grep("NeighbourCount",headerLines)
	NeighbourListDefLine=grep("NeighbourList",headerLines)
	RadiiDefLine=grep("Radii",headerLines)
	OriginDefLine=grep("int Origins",headerLines)

	coordAt=sub("^.*@([0-9]+$)","\\1",headerLines[CoordinatesDefLine])
	NCAt=sub("^.*@([0-9]+$)","\\1",headerLines[NeighbourCountDefLine])
	NLAt=sub("^.*@([0-9]+$)","\\1",headerLines[NeighbourListDefLine])
	RadAt=sub("^.*@([0-9]+$)","\\1",headerLines[RadiiDefLine])
	OriginAt=sub("^.*@([0-9]+$)","\\1",headerLines[OriginDefLine])
	
	if(	filetype=="binary"){
		d=matrix(readBin(con,what='numeric',n=nVertices*3,size=4,endian="big"),ncol=3,byrow=T)
		d=data.frame(d)
		colnames(d)=c("X","Y","Z")
		# next read the Neighbour Count
		cat(readLines(con,2)[2],"\n")
		NeighbourCount=readBin(con,what="integer",size=4,n=nVertices,endian="big")
		cat(readLines(con,2)[2],"\n")
		d$W=readBin(con,what='numeric',n=nVertices,size=4,endian="big")*2
		cat(readLines(con,2)[2],"\n")
		d$PointNo=seq(nrow(d))
		d$NeighbourCount=NeighbourCount
		# Note these numbers come in zero indexed, but I will want them 1-indexed
		# so add 1
		Neighbours=data.frame(Neighbour=readBin(con,what="integer",size=4,n=nEdges,endian="big")+1)
		cat(readLines(con,2)[2],"\n")
		
		Origin=NULL
		close(con)
		d[,1:4]=round(d[,1:4],digits=3)
		Neighbours$CurPoint=rep(seq(nVertices),d$NeighbourCount)
	} else {
		close(con)
		t=readLines(filename)
		dataStart=length(headerLines)+1
		coordStart=grep(paste("^\\s*@",sep="",coordAt,"\\s*"),t)+1
		NCStart=grep(paste("^\\s*@",sep="",NCAt,"\\s*"),t)+1
		NLStart=grep(paste("^\\s*@",sep="",NLAt,"\\s*"),t)+1
		RadiusStart=grep(paste("^\\s*@",sep="",RadAt,"\\s*"),t)+1
		OriginStart=grep(paste("^\\s*@",sep="",OriginAt,"\\s*"),t)+1

		d=read.table(filename,skip=coordStart-1,nrows=nVertices,col.names=c("X","Y","Z"),
			colClasses="numeric",na.strings=c("NA","ERR"))
		# SWC expects width rather than radius
		d$W=read.table(filename,skip=RadiusStart-1,nrows=nVertices)$V1*2
		# round to 3dp to avoid any surprises (like v small -ve numbers)
		d[,1:4]=round(d[,1:4],digits=3)
		d$PointNo=seq(nrow(d))
		d$NeighbourCount=read.table(filename,skip=NCStart-1,nrows=nVertices,colClasses="integer")$V1


		# Note these numbers come in zero indexed, but I will want them 1-indexed
		# so add 1
		Neighbours=read.table(filename,skip=NLStart-1,nrows=nEdges,col.names="Neighbour")+1
		# Figure out the current point for each of these neighbours
		Neighbours$CurPoint=rep(seq(nVertices),d$NeighbourCount)
		
		Origin=NULL
		if(any(OriginStart)) Origin=1+scan(filename,what=integer(1),skip=OriginStart-1,nlines=1,quiet=TRUE)
	}


	# check if there are actually any NAs - these should only be in XYZ
	InvalidPoints=which(apply(d[,c("X","Y","Z")],1,function(x) any(is.na(x))))
		
	if(OmitNAs && length(InvalidPoints)>0){
		# Remove all edges that reference an invalid point
		Neighbours=subset(Neighbours,
			!(CurPoint%in%InvalidPoints | Neighbour%in%InvalidPoints) )
		# find the remaining valid points NB valid points is not necessarily
		# the exact complement of invalid points - consider an end point which 
		# remains defined just distal to a point that doesn't transform.
		ValidPoints=sort(unique(Neighbours$CurPoint))

		# Restrict d to the valid points
		d=d[ValidPoints,]
		d$OldPointNo=d$PointNo # keep track of the old point numbers
		d$PointNo=seq(len=nrow(d)) # make the new ones	

		# need to figure out how to find the closest end point
		# to the original origin - hmm I think a better idea
		# would be to reroot the tree on the origin
# 		if(!is.null(Origin)){
# 			if(Origin
# 			Origin=d$OldP
		# For the moment, just assume that the origin is not deleted
		# or make it null if it is
		dO=subset(d,Origin==OldPointNo)
		if(nrow(dO)==1) Origin=dO$PointNo else Origin=NULL
		
		# figure out the new number of neighbours for remaining points
		d$NeighbourCount=table(Neighbours$CurPoint)
		
		# Renumber Neighbours data frame
		Neighbours$CurPoint=rep(d$PointNo,d$NeighbourCount)
		# Get the new neighbours by looking up against the old points
		# NB could have been written d$PointNo[match(...)]
		# but unnecessary since d$PointNo=seq(nrow(d))
		NewNeighbours=match(Neighbours$Neighbour,d$OldPointNo)
		# Verify that we have non-zero (or NA) entries
		stopifnot(all(NewNeighbours))
		Neighbours$Neighbour=NewNeighbours
	}
	
	return(list(PointList=d,EdgeList=Neighbours,Origin=Origin))
}


ParseAM3DToNeuron=function(datalist,filename,Force=FALSE,ProcessAllTrees=TRUE,Verbose=FALSE){
	# function to parse a an amira mesh 3D
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
	Neighbours=datalist$EdgeList
	SWCData=datalist$PointList
	nVertices=nrow(SWCData)
	EndPoints=subset(SWCData,NeighbourCount==1)
	BranchPoints=subset(SWCData,NeighbourCount>2)
	
	# a single unbranched segment will have 2 end points
	# if a branch point has 3 neighbours then it should increase
	# the number of free ends by one.
	# if 4 neighbours then free ends += 2 etc
	
	PredictedEndPoints=sum(BranchPoints$NeighbourCount-2)+2
	if(!Force && nrow(EndPoints)!=PredictedEndPoints){
		cat(paste("Mismatch between number of end points (",nrow(EndPoints),") and number predicted from branch statistics (",PredictedEndPoints,")\n"))
		#return (NULL)
	}

	# Start off with the start point as the first endpoint
	# we may change our mind if we try to parse multiple trees
	if(!is.null(datalist$Origin) && datalist$Origin%in%EndPoints$PointNo){
			StartPoint=datalist$Origin 
	} else {
			StartPoint=min(EndPoints$PointNo)
	}

	if(ProcessAllTrees){
		SubTrees=ParseEdgeListForAllSubTrees(Neighbours,Origin=datalist$Origin,Silent=!Verbose)
		nTrees=length(SubTrees)
		nPointsParsed=length(unique(unlist(SubTrees)))
		if(nPointsParsed != nrow(SWCData)){
			cat("unable to find valid subtrees\n")
			return(NULL)
		}
		
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
			
			SWCData$SubTree=rep(-1,nrow(SWCData))
			for(i in 1:length(SubTrees)){
				SWCData$SubTree[SWCData$PointNo%in%unique(unlist(SubTrees[[i]]))]=i
			}				
		}
# 		# Trim down to the chosen Tree
# 		# nb sort is important, because subtrees returns vertices
# 		# in the order that they are found ie not sorted
# 		SWCData=datalist$PointList[sort(SubTrees[[which.max(lengths)]]),]
# 		Neighbours=subset(datalist$EdgeList,CurPoint%in%SWCData$PointNo)
# 		nVertices=nrow(SWCData)
# 		# Replace the old point numbers with new ones which are contiguous
# 		OldNums=SWCData$PointNo
# 		SWCData$PointNo=seq(nVertices)
# 		Neighbours$CurPoint=rep(SWCData$PointNo,SWCData$NeighbourCount)
# 		Neighbours$Neighbour=SWCData$PointNo[match(Neighbours$Neighbour,OldNums)]
	} else {
		# I believe that my new ParseEdgeLists is robust to this now
		#if(StartPoint!=1) warning(paste(filename,": StartPoint is",StartPoint,"not 1 - this may break parsing routine"))
		SegList=ParseEdgeList(Neighbours,RootPoint=StartPoint)
	}
	
	SWCData$Parent=subset(Neighbours,!duplicated(CurPoint),sel=Neighbour)[,1]
	SWCData$Label=2
	firstFields=c("PointNo","Label","X","Y","Z","W","Parent")
	remainingFields=setdiff(names(SWCData),firstFields)
	SWCData=SWCData[,c(firstFields,remainingFields)]

	# Remove any Branch or End points that didn't make it into SegList
	PointsInSubTree=unique(unlist(SegList))
	EndPoints=subset(EndPoints,PointNo%in%PointsInSubTree)
	BranchPoints=subset(BranchPoints,PointNo%in%PointsInSubTree)
	
	if(length(SegList)>0){
		#OK There's at least one segment
		ParsedNeuron<-list(NeuronName=NeuronNameFromFileName(filename),
			InputFileName=filename,
			CreatedAt=Sys.time(),
			NodeName=Sys.info()["nodename"],
			InputFileStat=file.info(filename)[1,],
			InputFileMD5=md5sum(filename),
			NumPoints=nrow(SWCData),
			StartPoint=StartPoint, # NB I am assuming that this is always 1
			BranchPoints=BranchPoints$PointNo,
			EndPoints=EndPoints$PointNo,
			NumSegs=length(SegList),
			SegList=SegList,
			nTrees=nTrees,
			d=SWCData	)
		if(nTrees>1){
			# If there are multiple subtrees then make that data available 
			# as well.
			ParsedNeuron$SubTrees=SubTrees
		}
		
		return(ParsedNeuron)
	} else {
		cat("0 length Seglist for",filename,"\n")
		return(NULL)
	}
}


ParseEdgeList<-function(Nb,Silent=FALSE,Verbose=!Silent,RootPoint=1){
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
		if(length(lNeighboursFromPoint[[HeadPoint]])==0) {
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
			PossCurPoints=lNeighboursFromPoint[[PrevPoint]]			
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
			# Remove the relevant PrevPoint / CurPoint edges
			lNeighboursFromPoint[[PrevPoint]]<<-
				lNeighboursFromPoint[[PrevPoint]][!(lNeighboursFromPoint[[PrevPoint]]==CurPoint)]
			lNeighboursFromPoint[[CurPoint]]<<-
				lNeighboursFromPoint[[CurPoint]][!(lNeighboursFromPoint[[CurPoint]]==PrevPoint)]
			# What kind of point is the CurPoint? 
			# 1: an end point - terminate
			if(any(EndPoints==CurPoint)) return (0)
			Children=lNeighboursFromPoint[[CurPoint]]
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
	
	UnusedPoints=setdiff(unique(Nb$CurPoint),unique(unlist(SegList)))
	if(!Silent && any(UnusedPoints)){
		cat("The following points were not used in the SegList:",UnusedPoints,"\n")
		warning("")
	}
	
	return(SegList)
}	

ParseEdgeListForAllSubTrees<-function(Nb,Origin=NULL,Silent=T){
	# Repeatedly call ParseEdgeList with different roots until all points
	# accounted for
	if(is.null(names(Nb))) names(Nb)=c("Neighbour","CurPoint")
	pointsRemaining=sort(unique(Nb$CurPoint))
	EndPoints=pointsRemaining[table(Nb$CurPoint)==1]
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

	
	if(any(Components=="Axon")){	
		# WE'RE GOING TO BUILD A NEW NEURON FROM SCRATCH
		# 2. Get the axon data from the file
		
		cat("\n",AM3DFile,":")
		datalist=ReadAM3DData(AM3DFile)
		MyNeuron<-ParseAM3DToNeuron(datalist,AM3DFile,...)
		if(is.null(MyNeuron$SegList)){
			cat("Unable to Extract neuron from",AM3DFile,"\n")
			return(NULL)
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
	
	class(MyNeuron)=c("neuron",class(MyNeuron))
	return(MyNeuron)

}



WriteNeuronToAM<-function(ANeuron,AMFile=NULL,
	suffix="am",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE){
	# write out a neuron in the basic AmiraMesh format which is the native format
	# of amira for linesets (as opposed to the specialised skeletonize AM3D)
	# WriteAllSubTrees will write out all the stored subtrees in a neuron 
	# which has multiple subtrees (which is often true of ill-formed 
	# skeletonize neurons)
	
	if(is.null(AMFile))
		AMFile=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
	
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
	
	cat("Parameters \{\n",file=fc)
	cat("    ContentType \"HxLineSet\"\n",file=fc)
	cat("\}\n\n",file=fc)

	cat("Vertices \{ float\[3\] Coordinates \} = @1\n",file=fc)
	cat("Vertices \{ float Data \} = @2\n",file=fc)
	cat("Lines \{ int LineIdx \} = @3\n",file=fc)
	if(WriteAllSubTrees) cat("Vertices \{ float Data2 \} =@4\n",file=fc)
	cat("\n",file=fc)
	
	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(ANeuron$d[chosenVertices,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
	# Write the Radii
	cat("\n@2 #",nVertices,"width values\n",file=fc)
	# NB Divide width by 2
	#write(ANeuron$d$W/2,ncolumns=1,file=fc)
	#write.matrix(ANeuron$d$W/2,file=fc)
	write.table(ANeuron$d$W[chosenVertices]/2,col.names=F,row.names=F,file=fc)

	# Write the segment information
	cat("\n@3 #",nLinePoints,"line segements\n",file=fc)
	# nb have to -1 from each point because amira is 0 indexed
	# AND add -1 to each segment as a terminator
	tmp=lapply(SegList,function(x) cat(x-1,"-1 \n",file=fc) )
	#tmp=do.call("paste",ANeuron$SegList)
	#writeLines(tmp,con=fc)
	if(WriteAllSubTrees) {
		cat("\n@4 # subtrees\n",file=fc)
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

	if(is.null(colnames(d))){
		if(ncol(d)!=3) stop("Unable to identify X,Y,Z coordinates")
		colnames(d)=c("X","Y","Z")
	} else {
		colnames(d)=c("X","Y","Z")
	}
	
	cat("Writing to",AMFile,"\n")
	# Write the header
	cat("# AmiraMesh ASCII 1.0\n",file=AMFile)
	fc=file(AMFile,open="at") # ie append, text mode

	cat("# Created by WritePointsToAM - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	nVertices=nrow(d)
	cat("define Markers",nVertices,"Parameters {\nContentType \"LandmarkSet\",NumSets 1\n}\n",file=fc)
# 	cat("Parameters \{\n",file=fc)
# 	cat("    ContentType \"HxLineSet\"\n",file=fc)
# 	cat("\}\n\n",file=fc)

	cat("Markers \{ float\[3\] Coordinates \} = @1\n",file=fc)
	cat("\n",file=fc)
	
	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(d[,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
	close(fc)
}


WriteNeuronToAM3D<-function(ANeuron,AMFile=NULL,
	suffix="am3",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE){
	# write out a neuron in the specialised skeletonize AM3D format 
	# (as opposed to the basic AmiraMesh format which is the native format
	# of amira for linesets)
	# WriteAllSubTrees will write out all the stores subtrees in a neuron 
	# which has multiple subtrees (which is often true of ill-formed 
	# skeletonize neurons).  It will also add a data field that can be used
	# to visualised different subtrees eg by colouring
	
	if(is.null(AMFile))
		AMFile=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
	
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

	cat("\# Created by WriteNeuronToAM3D -",format(Sys.time(),usetz=T),"\n\n",file=fc)	
	cat("nVertices", nVertices,"\nnEdges",nEdgeList,"\n",file=fc)
	cat("define Origins 1\ndefine vertexTypeList 0\n\n",file=fc)
	
	cat("Parameters \{\n",file=fc)
	cat("    ContentType \"SkeletonGraph\"\n",file=fc)
	cat("\}\n\n",file=fc)

	cat("Vertices \{ float\[3\] Coordinates \} @1\n",file=fc)
	cat("Vertices \{ int NeighbourCount } @2\n",file=fc)
	cat("Vertices \{ float Radii \} = @3\n",file=fc)
	cat("EdgeData { int NeighbourList } @4\n",file=fc)
	cat("Origins { int Origins } @5\n",file=fc)
	cat("Vertices { int vertexTypeCounter } @6\n",file=fc)
	cat("vertexTypeList { int vertexTypeList } @7\n\n",file=fc)

	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(ANeuron$d[chosenVertices,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
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
	cat("\n@4 \#",nEdgeList,"bidirectional edges\n",file=fc)
	#NB -1 since Amira is 0 indexed
	write.table(EdgeList$Neighbour-1,col.names=F,row.names=F,file=fc)

	# Write the origin information NB -1 since 0 indexed
	cat("\n@5 \#n 1\n",file=fc)
	cat(ANeuron$StartPoint-1,"\n",file=fc)
	
	# Write the vertexTypeCounter information
	cat("\n@6 \#",nVertices,"\n",file=fc)
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


ParseAMSurfToContourList<-function(filename,RegionNames="ALL",RegionChoice="Inner",Verbose=FALSE){	
	# function to parse a an amira  HxSurface file
	# nb RegionChoice is a switch to allow the inneror outer region to
	# define the name of the region
	
	# Check for header confirming file type
	firstLine=readLines(filename,n=1)
	if(!any(grep("#\\s+hypersurface\\s+[0-9.]+\\s+ascii",firstLine,ignore.case=T,perl=T))){
		warning(paste(filename,"does not appear to be an Amira HyperSurface ASCII file"))
		return(-1)
	}
	
	t=readLines(filename)
	nLines=length(t)
	if(Verbose) cat(nLines,"lines of text to parse\n")
	
	# Find the start of the Vertices
	dataStart=grep("^\s*Vertices\s*",t)[1]
	if(Verbose) cat("Data start line =",dataStart,"\n")
	headerLines=t[seq(dataStart-1)]
	trim=function(x) sub('^\\s+', '', sub('\\s+$', '', x, perl = TRUE), perl = TRUE)
	getfield=function(fName,textLines=headerLines,pos=2) unlist(strsplit(trim(textLines[grep(fName,textLines)]),"\\s+",perl=TRUE))[pos]
	nVertices=as.numeric(getfield("Vertices",t[dataStart],2))
	if(Verbose) cat("nVertices =",nVertices,"\n")

	d=list()
	d$Vertices=read.table(filename,skip=dataStart,nrows=nVertices,col.names=c("X","Y","Z"),colClasses=rep("numeric",3))

	# round to 3dp to avoid any surprises (like v small -ve numbers)
	d$Vertices=round(d$Vertices,digits=3)
	d$Vertices$PointNo=seq(nrow(d$Vertices))
	if(Verbose) cat("Finished processing Vertices\n")
	
	# Now read in Triangles that define patches:
	linesSkipped=dataStart+nVertices-1
	remainingLines=t[(dataStart+nVertices):nLines]
	PatchDefLine=grep("^\s*Patches\s*",remainingLines,perl=TRUE)
	cat("PatchDefLine =",PatchDefLine,"\n")
	nPatches=as.numeric(getfield("Patches",remainingLines[PatchDefLine],2))
	cat("nPatches =",nPatches,"\n")
	PatchStarts=grep("^\\s*{",remainingLines[PatchDefLine:length(remainingLines)],perl=TRUE)+PatchDefLine-1
	if(length(PatchStarts)>nPatches) PatchStarts=PatchStarts[1:nPatches]
	PatchEnds=grep("^\\s*}",remainingLines[PatchDefLine:length(remainingLines)],perl=TRUE)+PatchDefLine-1
	if(length(PatchEnds)>nPatches) PatchEnds=PatchEnds[1:nPatches]
	#return(d)
	TriangleDeflines<<-grep("Triangles",remainingLines)
	#myreadtable<-function(...) scan(...)
	for(i in 1:nPatches){
		if(!any(TriangleDeflines[i])){
			warning(paste("Unable to find Triangle number in patch",i,"in",filename,"\n"))
			return (-1)
		}
		if(Verbose) cat("TriangleDefline =",TriangleDeflines[i],"\n")
		PatchHeader<<-remainingLines[PatchStarts[i]:TriangleDeflines[i]]
		if(Verbose) cat("PatchHeader is",length(PatchHeader),"lines long\n")
		# note use of RegionChoice to switch naming between inner and outer
		RegionName=getfield(paste(RegionChoice,"Region",sep=""),PatchHeader,2)
		#RegionName=getfield("InnerRegion",PatchHeader,2)
		nTriangles=as.numeric(getfield("Triangles",PatchHeader,2))
		if(nTriangles<0 || nTriangles>100000){return(-1)}
		if(Verbose) cat("nTriangles =",nTriangles,"for patch =",i,"\n")
		# Check if we want to load in this region
		if( ("ALL"%in%RegionNames) || (RegionName%in%RegionNames) ){
			# check if we have already loaded a patch in this name
			if(RegionName%in%names(d)){
				#return(d)
				# add to the old patch
				if(Verbose) cat("Adding to patch name",RegionName,"\n")
				d[[RegionName]]=rbind(d[[RegionName]],read.table(filename,skip=linesSkipped+TriangleDeflines[i],nrows=nTriangles,col.names=c("V1","V2","V3")))
				#d[[RegionName]]=rbind(d[[RegionName]],read.table(filename,skip=linesSkipped+TriangleDeflines[i],nrows=nTriangles,col.names=c("V1","V2","V3")))
			} else {
				# new patch
				if(Verbose) cat("Making new patch name",RegionName,"\n")
				# scan no quicker in these circs, problem is repeated file 
				# access
				#d[[RegionName]]=as.data.frame(matrix(scan(filename,skip=linesSkipped+TriangleDeflines[i],nlines=nTriangles),ncol=3,byrow=T))
				d[[RegionName]]=read.table(filename,skip=linesSkipped+TriangleDeflines[i],nrows=nTriangles,col.names=c("V1","V2","V3"))
			}
		}
	}
	d$RegionList=setdiff(names(d),"Vertices")
	return(d)
}

WriteAmiraLandmarkSet=function(filename,d){
	cat("# AmiraMesh ASCII 1.0\ndefine Markers",nrow(d),"Parameters {\nContentType \"LandmarkSet\",NumSets 1\n}\n",file=filename)
	cat("Markers { float[3] Coordinates } = @1\n@1\n",file=filename,append=T)
	write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
}

WriteHxSurface=function(filename,Vertices,Indices=NULL,
		material=sub("^([^.]+)\\..*","\\1",basename(filename))){
		
	cat("\# HyperSurface ASCII\nParameters {\n",file=filename)
	fc=file(filename,open="at") # ie append, text mode
	#cat("\tMaterials{\n\t\tInterior {\n\t\t\tid 0\n\t\t}\n",file=fc)
	cat("\t{color 0.83562 0.78 0.06,\nName \"",sep="",material,"\"}",file=fc)
	cat("\t\t",sep="",material," {\n\t\t\tid 1\n\t\t}\n",file=fc)
	cat("\t}\n",file=fc)
	cat("\tBoundaryIds {\n\t\t\Id0 {\n\t\t\tId 0,Info \"undefined\",Color 0.6 0.6 0.6\n\t\t}\n\t\tname \"BoundaryConditions\"\n\t}\n}\n",file=fc)


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
		cat("\# AmiraMesh ASCII 1.0\n\n",file=filename)
		fc=file(filename,open="at") # ie append, text mode
		lattice=apply(d$eval.points,2,length)
		cat("define Lattice",lattice,"\n",file=fc)
		cat("define Coordinates",sum(lattice),"\n\n",file=fc)
		cat("Parameters {CoordType \"rectilinear\"}\n\n",file=fc)
		cat("Lattice { float ScalarField } = @1\n",file=fc)
		cat("Coordinates { float xyz } = @2\n\n",file=fc)
		cat("\@1\n",file=fc)
		
		cat(as.vector(d$estimate),file=fc)
		cat("\@2\n",file=fc)
		cat(d$eval.points[,1],"\n",file=fc)
		cat(d$eval.points[,2],"\n",file=fc)
		cat(d$eval.points[,3],"\n",file=fc)
}

Write3DDensityToAmiraLattice<-function(filename,dens,ftype=c("binary","text"),
	dtype=c("float","byte", "short", "int", "double"),WriteNrrdHeader=FALSE){
	# Produces a lattice format file -
	# that is one with a regular x,y,z grid
	# Can also write a detached Nrrd header that points to the AmiraMesh
	# data to allow it to be opened by a nrrd reader
	ftype=match.arg(ftype)
	if(ftype=='text') cat("\# AmiraMesh ASCII 1.0\n\n",file=filename)
	else cat("\# AmiraMesh 3D BINARY 2.0\n\n",file=filename)
		
	fc=file(filename,open="at") # ie append, text mode
	cat("# Created by Write3DDensityToAmiraLattice - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	if(!is.list(dens)) d=dens else d=dens$estimate
	# Find data type and size for amira
	dtype=match.arg(dtype)	
	dtypesize<-c(4,1,2,4,8)[which(dtype==c("float","byte", "short", "int", "double"))]
	# Set the data mode which will be used in the as.vector call at the
	# moment that the binary data is written out.
	if(dtype%in%c("byte","short","int")) dmode="integer"
	if(dtype%in%c("float","double")) dmode="numeric"
	
	
	#lattice=apply(d$eval.points,2,length)
	lattice=dim(d)
	cat("define Lattice",lattice,"\n",file=fc)

	cat("Parameters { CoordType \"uniform\",\n",file=fc)
	# note Amira's definition for the bounding box:
	# the range of the voxel centres.
	# So eval.points should correspond to the CENTRE of the
	# voxels at which the density is evaluated
	cat("\t\# BoundingBox is xmin xmax ymin ymax zmin zmax\n",file=fc)
	BoundingBox=NULL
	if(!is.null(attr(dens,"BoundingBox"))){
		BoundingBox=attr(dens,"BoundingBox")
	} else if(!is.null(d$eval.points)){
		BoundingBox=as.vector(apply(d$eval.points,2,range))
	}
	if(!is.null(BoundingBox)) cat("\t BoundingBox",BoundingBox,"\n",file=fc)
	cat("}\n\n",file=fc)
	
	cat("Lattice {",dtype,"ScalarField } = @1\n\n",file=fc)

	cat("\@1\n",file=fc)
	
	#cat(str(as.vector(d)))
	close(fc)

	# Write a Nrrd header to accompany the amira file if desired
	# see http://teem.sourceforge.net/nrrd/
	if(WriteNrrdHeader) {
		nrrdFilename=paste(filename,sep=".","nhdr")
		cat("NRRD0004\n",file=nrrdFilename)
		fc=file(nrrdFilename,open="at") # ie append, text mode
		nrrdType=ifelse(dtype=="byte","uint8",dtype)
		
		cat("encoding:", ifelse(ftype=="text","text","raw"),"\n",file=fc)
		cat("type: ",nrrdType,"\n",sep="",file=fc)
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
		writeBin(as.vector(d,mode=dmode),fc,size=dtypesize,endian='big')
		close(fc)
	}
}

Read3DDensityFromAmiraLattice<-function(filename,Verbose=FALSE){

	fc=file(filename,'rb')
	headerLines=NULL		
	while ( ( nextLine<-readLines(fc,1)) !="@1") {headerLines<-c(headerLines,nextLine)}
	
	# Figure out if the file is in binary format or not
	if(any(grep("^\\s*#\\s+amiramesh\\s+3d\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
			binary = TRUE
	} else binary=FALSE

	# Find the position of the header lines defining the lattice
	# this clearly assumes that the relevant data is in position
	# @1 in the file - which may or may not be the case.
	LatticeDefLine=grep("define\\s+lattice",headerLines,ignore.case=TRUE,perl=TRUE)
	LatticeTypeDefLine=grep("^Lattice.*}\\s*[=]{0,1}\\s*@1",headerLines,ignore.case=TRUE,perl=FALSE)
	#cat("LatticeTypeDefLine = ",LatticeTypeDefLine)
	LatticeBoundsLine=grep("^[^\#]*BoundingBox",headerLines,ignore.case=TRUE,perl=TRUE)
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
	dataTypeName=sub(".*\\{\\s*(\\w+)\\s+.?*\\s*\\}.*","\\1",ldl,ignore.case=T)
	# check if the data is encoded in some way
	#"Lattice { byte Labels } @1(HxByteRLE,391697)"
	postAt=strsplit(ldl,"@1")[[1]][2]
	dataEncoding=""
	if(isTRUE(grep("\\(Hx.?*\\)",postAt)==1)){
		postAt=sub(".?*\\(([^)]+)\\).?*","\\1",postAt)
		dataEncoding=toupper(strsplit(postAt,",")[[1]][1])
		dataLength=as.integer(strsplit(postAt,",")[[1]][2])
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
	dataTypes=data.frame(name=I(c("byte", "short", "int", "float", "double", "complex")),
			size=c(1,2,4,4,8,NA),what=I(c(rep("integer",3),rep("numeric",2),NA)),
			signed=rep(c(FALSE,TRUE),c(1,5)) )
	i=which(dataTypes$name==dataTypeName)
	if(!any(i==1:6)){
		close(fc)
		stop("Unrecognised data type")
	}
	
	if(Verbose) cat("dataLength =",dataLength,"dataType =",dataTypes$what[i],"size=",dataTypes$size[i],"\n")
	if(binary){
		if(dataEncoding=="HXBYTERLE"){
#			d=readBin(fc,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],endian='big',signed=T)
			d=readBin(fc,what=integer(0),n=dataLength,size=1)
			d=DecodeRLEBytes(d)
		} else if(dataEncoding==""){
			d=readBin(fc,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
				signed=dataTypes$signed[i],endian="big")
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

DecodeRLEBytes<-function(ba){
	# Expects an integer array
	# Structure is that every odd byte is a count
	# and every even byte is the actual data
	# So 127 0 127 0 127 0 12 0 12 1 0
	# I think that it ends with a zero count
	stopifnot(ba[length(ba)]==0)
	m=ba[-length(ba)]
	dim(m)=c(2,length(m)/2)
	# for some reason a count of 1 is encoded as -127
	m[1,m[1,]==-127]=as.integer(1)
	#nBytes=sum(m[,1])
	# Check that there are no remaining negative counts
	stopifnot(sum(m[1,]<0)==0)
	rval=rep(m[2,],times=m[1,])
	return(rval)
}

ReadRLEBytes<-function(con,length,offset=0){
	# Expects a connection + length/offset
	if(!is.connection(con)) con=file(con,open='rb')
	skip(con,offset)
	ba=readBin(con,n=length,what=integer(),size=1,signed=T)
	return(DecodeRLEBytes(ba))
}

ReadAmiraLandmarks<-function(filename){

	fc=file(filename,'rb')
	headerLines<<-NULL		
	while ( ( nextLine<-readLines(fc,1)) !="@1") {headerLines<<-c(headerLines,nextLine)}
	
	# Figure out if the file is in binary format or not
	if(any(grep("^\\s*#\\s+amiramesh\\s+3d\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
			binary = TRUE
	} else binary=FALSE

	# Find the position of the header line defining the # of markers
	MarkersDefLine=grep("define\\s+markers",headerLines,ignore.case=TRUE,perl=TRUE)
	if(!any(MarkersDefLine)) return(-1)
#		if(!any(LatticeBoundsLine)) return(-1)
	
	strtrim<-function(str){
			str<-sub('\\s+$', '', str, perl = TRUE) ## Perl-style white space
			sub('^\\s+', '', str, perl = TRUE) ## Perl-style white space
	}

	# fetch the number of markers
	markerLineComponents=unlist(strsplit(headerLines[MarkersDefLine]," "))
# 	cat("markerLineComponents=",markerLineComponents)
# 	cat(markerLineComponents[length(markerLineComponents)])
	nMarkers=as.numeric(markerLineComponents[length(markerLineComponents)])
	dataLength=nMarkers*3
# 	cat("nMarkers =",nMarkers,"dataLength =",dataLength)
	
	# get the number of data sections
	nDataSections=length(grep("@[[:digit:]]+",headerLines))
# 	cat("nDataSections =",nDataSections)
	
	if(nDataSections>1){
		d=list()
		for(i in 1:nDataSections){
			d[[i]]=matrix(scan(fc,what=numeric(),nmax=dataLength,quiet=TRUE),ncol=3,byrow=TRUE)
			readLines(fc,n=2)
		}
	} else {
		d=matrix(scan(fc,what=numeric(),nmax=dataLength,quiet=TRUE),ncol=3,byrow=TRUE)
	}			
	close(fc)
	return(d)
}
