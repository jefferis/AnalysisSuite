# ReadNeuronFromAsc2.s
# Functions that allow axon data in Neurolucida .asc files to be read
# directly into an R Neuron object.  

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

#NB This is some of the oldest and ugliest code - but it does work so I
#have not gone back and re-written it.  It is not clear how accurate all
#the notes below are.
#GJ 2007-03-22

# They still need to be integrated into the
# main program (called by some script etc.) and it will now make sense to do
# all the parsing of .asc files in one go - i.e. look for contour data etc.
# In general, I would assume that the default behaviour should be update
# and replace at this stage, but can think about that.
# I think the routines to set the LH and MB segments could be more integrated
# now:
#   - the LH can be set as the point closest to the Lateral Horn Entry Point
#   - the MB branches could perhaps be set as the 1st branch point after root
#   + perhaps the last branch point can still be set manually - or use
#   MB contours to automate.
# GJ 020623
# Parents to be set if a segment consisted on only 2 points
# Have decided to generalise these routines by splitting up into the following steps
# 1. read in ascii data  WORKS
# 2. get axon data from that  WORKS
# 3. parse axon data and create a neuron  WORKS
# 4. Get Marker Points and Orientation info and orient axon data WORKS
# 5. Get LH and MB contour data  WORKS - except for User Line 1 problem
# ----------------------------------- Dealt with that by replacing in data files 
# 6. Scale LH and MB data WORKS
# 7. return finished neuron   EASY!
# 8. Also can now update only selected attributes (Axon, LH, MB)

# Modifications 020628
# Modifications 020715 - fixed a bad bug in ParseAxonData which caused incorrect

# PRESSING ISSUES : Decide on New Neuron Object format (more or less done)
# ALTER MAIN CODE TO REFLECT NEW NEURON Object format - mostly done

# source("Greg Data:AnalysingTraceFiles:R:ReadNeuronFromAsc2.s")
# source(file.path(CodeDir,"ReadNeuronFromAsc2.s"))

# Main function which returns a full neuron based on the supplied filename
# which should include full path if necessary
# optional Components indicates which data to read in
# optional OldNeuron can pass a Neuron object to which LH or MB data is added/replaced
ReadNeuronFromAsc<-function(AscFile,Components="Axon",OldNeuron=NULL,ReOrient=F,WithConvexHull=T){
    # nb if WithConvHull ==T then use convex hull to reduce number of data 
    # points in LH and MB contours 
    # 1. read Neurolucida ascii data into a text array for faster processing.
    AscData<-ReadAscData(AscFile)
    
    if(any(Components=="Axon")){	
	# WE'RE GOING TO BUILD A NEW NEURON FROM SCRATCH
	# 2. Get the axon data from that 
	# Pass the file name just to print in error messages
	AxonData<-GetAxonData(AscData,AscFile)
	MyNeuron<-ParseAxonData(AxonData,AscFile)
	# Read in Orientation data from Neuron and figure out
	# which way up and which orientation it is in
	# Also add a set of scale information
	MyNeuron<-OrientNeuron(MyNeuron,AscData,AscFile,ReOrient=ReOrient)

	# Get Marker Points eg Lateral Horn Entry Point
	#if(!is.null(MyNeuron$OrientInfo)) MyNeuron<-GetMarkerPoints(MyNeuron,AscData,AscFile)    
    } else {
	# WE'RE GOING TO USE A SUPPLIED NEURON AS A BASE 
	if(is.null(OldNeuron)) stop(
	    "Must supply either an OldNeuron or read in AxonData in ReadNeuronFromAsc")
        #if(is.null(OldNeuron$OrientInfo)) stop(
	#   paste("Neuron",OldNeuron$NeuronName,"is in an old format without $OrientInfo and must be entirely reloaded from its Neurolucida .asc file"))
	MyNeuron<-OldNeuron
    }
    if(any(Components=="LH")){
	# Get Raw LH Contour Data (returns a data frame)
	RawLHData<-GetContourData(AscData,"LH",WithConvexHull=WithConvexHull)
	#OK Lets scale/reorient these data based on the info
	#in OrientInfo which was filled in when the axon data was read in
	RawLHData[,c("X","Y","Z")]<-InvertAndReorderData(RawLHData[,c("X","Y","Z")],MyNeuron$OrientInfo)
	# And add to my neuron
	if(is.null(MyNeuron$LH)){
	    # NB the double
	    MyNeuron<-c(MyNeuron,list(LH=list(d=RawLHData)))
	} else {
	    MyNeuron$LH$d<-RawLHData
	}
    }
    
    if(any(Components=="MB")){
	RawMBData<-GetContourData(AscData,"MB",,WithConvexHull=WithConvexHull)
	RawMBData[,c("X","Y","Z")]<-InvertAndReorderData(RawMBData[,c("X","Y","Z")],MyNeuron$OrientInfo)
	if(is.null(MyNeuron$MB)){
	    MyNeuron<-c(MyNeuron,list(MB=list(d=RawMBData)))
	} else {
	    MyNeuron$MB$d<-RawMBData
	}
    }
    
    class(MyNeuron)=c("neuron",class(MyNeuron))
    return(MyNeuron)

}

# Returns the axon data from the file AscFile in a string array
# one string for each line
ReadAscData<-function(AscFile){        
    # First off there is a problem here 
    # Readlines inserts extra blank lines for files in IBM format
    # this is presumably because of CR/LF choices for the line terminator
    # Mac format files (e.g. from Alpha) are fine
    AscData<-readLines(AscFile,n=-1)
    # Check that this file actually has some data in it!
    if(length(AscData)<2) stop(paste("No data in file",AscFile))

    # Find any blank lines (deals with problem mentioned at top)
    BlankLines<-grep("^[ ]*$",AscData)
    CommentLines<-grep("^[ ]*;",AscData)
    # Return file without blank lines
    return(AscData[-c(BlankLines,CommentLines)]) 
}



# Utility Function which does most of the work reading in contours
# Returns the data block in a string array, one string for each line,
# from the passed array AscData (which was presumably read in by ReadAscData)
# of type BlockType e.g. Axon or Contour
GetContourData<-function(AscData,ContourType,minContourLength=5, WithConvexHull=T){
    # ContourType indicates whether we want LH, MB etc.
    # Any Contour with fewer points than minContourLength is probably a dud
    # Use Convex hull to reduce the number of points in any single contour
    
    # Will be a problem if the contour type is anything other than
    # "Contour Type" (obviously!) - would like to be able to give a list
    # of alternatives since sometimes I used "User Line 1" and sometimes "LH"
    # eg ig ContourType = "LH" searches for:
    # ("LH" at the start of a line
    StartLineString<-paste("^[(]\"",ContourType,"\"",sep="")
    # If the contour is "closed" then there are 4 lines till the start
    # If it is "open" then ther are 3 and there will be a problem!
    StartLines<-grep(StartLineString,AscData)+4
    if (length(StartLines)<1) warning(paste("No contours in this file",ContourFile))
    
    # Find lines ending contours
    # -1 since this line follows the last data point
    EndLines<-grep("^)",AscData)-1
    # nb this will find EndLines for all blocks ie both LH and MB contours
    # so have to pair up the correct EndLines with thte StartLines
    # Need to find the number which is just larger than each StartLine
    EndLines<-EndLines[sapply(StartLines,function(x){min(which(EndLines>x))})]
    if (length(StartLines)!=length(EndLines)){
	stop(paste("Mismatched",ContourType,"startlines and endlines"))
    }
    
    # Figure out how many lines of data for each individual contour
    LinesToRead<-EndLines-StartLines+1
    # Remove any contours with small numbers of data points
    StartLines<-StartLines[which(LinesToRead>=minContourLength)]
    EndLines<-EndLines[which(LinesToRead>=minContourLength)]
    LinesToRead<-EndLines-StartLines+1
    
    # Initialise an empty matrix to which data will be concatenated
    ContourPoints<-matrix(nrow=0,ncol=4)
    for( i in 1:length(StartLines) ){
	# Modified 031113 to use scan
	ThisBlock<-AscData[StartLines[i]:EndLines[i]]
	#ThisBlock=gsub("[(]+(.*)[)]+.*","\\1",ThisBlock)
	ThisBlock=gsub("[()]+","",ThisBlock)
	tCon=textConnection(ThisBlock)
	allpoints=scan(tCon,quiet=T,comment.char=";",na.strings=c("NA","ERR"))
	close(tCon)
	allpoints=matrix(allpoints,ncol=4,byrow=T)
	# replace the 4th column
	allpoints[,4]=rep(i,nrow(allpoints))
	#tb<-sapply(strsplit(ThisBlock,split="[ (]+"),function(x){x[2:4]})
	#allpoints<-matrix(as.numeric(tb),,3,byrow=T)
	
	#allpoints<-cbind(allpoints,rep(i,length(allpoints[,1])))

	# added by GJ 13 Aug 2004 to prevent problems
	# arising from imperfect warping which can result in NAs close to the edge.
	allpoints=na.omit(allpoints)
	if(WithConvexHull){
	    # Find the convex hull surrounding the points on this contour
	    # that should reduce the number a bit.
	    oldlen<-length(allpoints[,1])
	    #print(oldlen)
	    BoundingPoints<-chull(allpoints[,1:2])
	    # The sort is because chull can select a point which was not the
	    # original head of the contour as the new head
	    # this can be a problem because it causes extra lines
	    # to be shown in plotneuron3d which relies on the points being
	    # in precisely the order around the circle.
	    BoundingPoints<-sort(BoundingPoints)
	    reducedpoints<-allpoints[BoundingPoints,]
	    newlen<-length(reducedpoints[,1])
	    ContourPoints<-rbind(ContourPoints,reducedpoints)
	} else {
	    ContourPoints<-rbind(ContourPoints,allpoints)
	}
    } # end for 
    
    #OK turn the acquired data into a dataframe
    ContourPoints<-data.frame(X=ContourPoints[,1],Y=ContourPoints[,2],
	Z=ContourPoints[,3],ContourID=ContourPoints[,4])
    
    # Now there's a problem: occasionally contour files have contours
    # which are not numbered serially up or down the Z axis
    # i.e. when the tracer doesn't draw the contours in order.
    # To get round this, I need to sort the contours by Z value
    NewOrder<-order(ContourPoints$Z)
    ContourPoints<-ContourPoints[NewOrder,]
    
    return(ContourPoints)
}

# This little function gets any points defined in AscFile
# such as Lateral Horn Entry Point
GetMarkerPoints<-function(ANeuron,AscData,AscFile){
    LHEPLine<-na.omit(c(grep("LHEntryPoint",AscData)+1,grep("Marker 3",AscData)+1))[1]
    LHEPoint<-as.numeric(unlist(strsplit(AscData[LHEPLine],split="[ (]+"))[2:4])

    if ( any(is.na(LHEPoint)) ) warning(
	paste("Cant find Lateral Horn Entry Point Marker in",AscFile))
    # Scale/reorder point axes if necessary
    LHEPoint<-InvertAndReorderData(LHEPoint,ANeuron$OrientInfo)
    MarkerPoints<-list(LHEPoint=LHEPoint)
    if( is.null(ANeuron$MarkerPoints) ){
	### Need to create a new object for Marker Points in this Neuron
	### nb need to add as a list
	ANeuron<-c(ANeuron,list(MarkerPoints=MarkerPoints))
    } else {
        ANeuron$MarkerPoints<-MarkerPoints
    }
    return(ANeuron)
}

# This function gets the data points from the file that will allow the
# orientation of the image/neuron to be figured out
# so that the axes can be swapped/inverted in due course
# AT THE MOMENT CONTAINS RE-ORIENTING CODE FOR AXON DATA!!
OrientNeuron<-function(ANeuron,AscData,AscFile,ReOrient=T){

   
    if(ReOrient){
	#OK Now get the Orientation Markers and the LHEntryPoint
	
	AVMLine<-na.omit(c(grep("AVM",AscData)+1,grep("Marker 1",AscData)+1))[1]
	AVMPoint<-as.numeric(unlist(strsplit(AscData[AVMLine],split="[ (]+"))[2:4])
	if ( any(is.na(AVMPoint)) ){
	    warning(
	    paste("Cant find AVM (AnteroVentroMedial) Corner in",AscFile))
	    return(ANeuron)
	}
	

	AVLLine<-na.omit(c(grep("AVL",AscData)+1,grep("Marker 2",AscData)+1))[1]
	AVLPoint<-as.numeric(unlist(strsplit(AscData[AVLLine],split="[ (]+"))[2:4])
	if ( any(is.na(AVLPoint)) ){
	    warning(
	    paste("Cant find AVL (AnteroVentroLateral) Corner",AscFile))
	    return(ANeuron)
	}
	# OK Now get the new axes nb the largest and smallest values on each 
	# the axes of ANeuron$d (the array describing the points in the axon tree)
	# are also passed to GetRealAxes
	MyPointList<-c(list(AVMPoint=AVMPoint,AVLPoint=AVLPoint),
	    lapply(ANeuron$d[c("X","Y","Z")],range))
     
	NewAxisInfo<-GetRealAxes(MyPointList)
	# OK Now that we've figured out how to treat the axes, rescale/order axon data axes
	# and the AVM/L Marker Points
	ANeuron$d[,c("X","Y","Z")]<-InvertAndReorderData(ANeuron$d[,c("X","Y","Z")],NewAxisInfo)
	AVMPoint<-InvertAndReorderData(AVMPoint,NewAxisInfo)
	AVLPoint<-InvertAndReorderData(AVLPoint,NewAxisInfo)
    # NB AxonOriented is a flag to indicate that we have rotated the axon data
    # NB making AVMPoint and AVLPoint into their own list is so that these 2 arrays 
    # don't get unlisted into AVMPoint1,AVMPoint2 etc.
	OrientInfo<-c(AxonOriented=T,list(AVMPoint=AVMPoint,AVLPoint=AVLPoint),NewAxisInfo)
    } else {
	# GJ: 28 March 2004
	# I don't want to rescale the new trace files since they should
	# all be in the correct orientation
	# This is a hack to prevent re-orientation of trace files
        NewAxisInfo=list(Scl=c(1,1,1),NewAxes=c(1,2,3))
	OrientInfo<-c(AxonOriented=T,list(AVMPoint=NA,AVLPoint=NA),NewAxisInfo)
    }
    
    if( is.null(ANeuron$OrientInfo) ){
	### Need to create a new object for orientation info since this Neuron
	  # doesn't yet have any orientation info
	ANeuron<-c(ANeuron,list(OrientInfo=OrientInfo))
    } else {
        ANeuron$OrientInfo<-OrientInfo
    }

    return(ANeuron)
}

# Expects a point or a matrix of X,Y and Z values
InvertAndReorderData<-function(DataArray,OrientInfo){
    # First of all check if we have a single point or a matrix
    if(is.null(nrow(DataArray))){
	# we have a point, so turn it into an array
	DataArray<-matrix(DataArray,nrow=1,ncol=3)
    } # ELSE WE ALREADY HAVE A MATRIX
    # Rescale These Axes
    for (thisax in 1:3){
	DataArray[,thisax]<-OrientInfo$Scl[thisax]*DataArray[,thisax]
    }
    # Swap x and y axes if necessary
    if (OrientInfo$NewAxes[1]==2){
	DataArray[,c(1,2)]<-DataArray[,c(2,1)]
    }
    return(DataArray)
}

# Given appropriate input, this function will produce
# two vectors which indicate whether the current axes need to 
# be inverted (Scl)
# and then what the new axes are (and in particular whether
# x and y need to be swapped)
GetRealAxes<-function(PointList){
    # expects list containing X,Y and Z coords of Two Corner Markers
    # AVM = Anterior Ventral Medial corner
    # AVL = Anterior Ventral Lateral corner
    # and the range of a set of X Y and Z points from 
    # either the contour or the neuron in question
    
    AVM<-PointList$AVMPoint
    AVL<-PointList$AVLPoint
    RangeOldAxes<-matrix(unlist(c(PointList$X,PointList$Y,PointList$Z)),3,2,byrow=T)
    
    # OK Find the mediolateral (i.e. X) axis
    # This is going to be the one where we find the greatest difference
    # between AVM and AVL 
    DeltaAVMAVL<- diff(rbind( AVM, AVL))
    MLAxis<-which(abs(DeltaAVMAVL)==max(abs(DeltaAVMAVL)))
    # The 3rd axis is always z / PosteroAnterior, and the 
    # remaining axis which isn't ML 
    # must be ventrodorsal
    VDAxis<-ifelse(MLAxis==1,2,1)
    PAAxis<-3
    
    # OK AXES DONE!  NOW MOVING ON TO WHICH WAY THEY ARE POINTING
    
    # Now which way is the present coord system pointing
    # Medial should be < Lateral
    Scl<-c(0,0,0)
    if (AVM[MLAxis]<AVL[MLAxis]){
	Scl[MLAxis] <- +1
    } else {
	#wrong way so
      	Scl[MLAxis] <- -1
    }
    
    #OK Next, which way  is up? ( I want the posterior to be the origin)
    
    # The abs diff between the min and the anterior
    # should be greater than the abs diff between the max and the anterior 
    # nb the z axis is ALWAYS axis 3
    if( abs(RangeOldAxes[3,1]-AVM[3])>abs(RangeOldAxes[3,2]-AVM[3]) ){
	# All Well
	Scl[3]<-1
    } else {
        Scl[3]<--1
    }
    
    # OK Now for the VD axis
    # ventral should be closer to min than max
    # nb I had a nasty bug here in which I was accessing
    # AVM[2] and AVL[2] instead of AVM[VDAxis] etc.
    if( abs(RangeOldAxes[VDAxis,1]-AVM[VDAxis])<abs(RangeOldAxes[VDAxis,2]-AVM[VDAxis]) ){
	# All Well
	Scl[VDAxis]<-1
    } else {
	Scl[VDAxis]<--1
    }    
    return (list(Scl=Scl,NewAxes=c(MLAxis,VDAxis,PAAxis)))
}


# Returns the axon data in a string array, one string for each line,
# from the passed array AscData (which was presumably read in by ReadAscData)
GetAxonData<-function(AscData,AscFile){        

    # Find the start of the Axon Block
    StartLine<-grep("(Axon|Dendrite)",AscData)+1

    # Check the file seems sensible
    if (is.na(StartLine)){
	stop(paste("Doesn't seem to be any axon data in",AscFile))
    }
    if(length(StartLine)>1){
	stop(paste("I don't know how to deal with >1 blocks of axon data in",AscFile))
    }
    
    # Find the end of the Axon Block by looking for the next close bracket 
    # at the beginning of a line (Hmm is this correct, what if the neuron contains
    # dendrites, soma etc.?)
    # The line in front is chosen since there first level of brackets
    # only enclose the declaration that this is an axon, not
    # actual structural components.
    AxonBlockLength<-grep( "^)", AscData[StartLine:length(AscData)] )-1
    
    # Check that we actually found the end
    if (length(AxonBlockLength)<1 || is.na(AxonBlockLength)){
	stop(paste("Axon data block doesn't seem to terminate in",AscFile))
    }

    AxonData<-AscData[StartLine:(StartLine+AxonBlockLength)]

    # Return the block of data corresponding to the axon
    return(AxonData) 
}


# This function expects a block of axon data from ReadAxonData
# It returns a new Neuron object
ParseAxonData<-function(AxonData,AscFile){
    cat(AscFile,"")
    
    # First off need to find the "special lines"
    # These are 
    # (        # new BRANCH POINT ie start of a subtree
    # Normal   # JUST GET rid of these
    # |        # SPLIT - start of a 2nd or 3rd branch
    # )        # finish off a subtree
    
    # Find the Instruction Lines
    # This is rather peculiar but for grep "\(" doesn't work need to do "[(]"
    NewBranchLines<-grep("^[ ]*[(]+[ ]*$",AxonData)
    SplitLines<-grep("^[ ]*[|]+[ ]*$",AxonData)
    AllBranchLines<-sort(c(NewBranchLines,SplitLines))
    SubtreeEndLines<-grep("^[ ]*\\)",AxonData)
    InstructionLines<-sort(c(AllBranchLines,SubtreeEndLines))
    # Hmm actually just check that no lines are duplicated
    if( length(unique(InstructionLines)) < length(InstructionLines) ){
	stop(paste("Duplicated Instruction lines in function ParseAxonData() for file",AscFile))
    }
    # Find the lines actually containing Data
    # DataLines<-grep("[(]([0-9]|[ .-])*\\)",AxonData)
    # Modified grep 031112 to remove lines indicating spines which look like
    #  >( X Y Z W); Spine
   # DataLines<-grep("^[[:space:]]*[(]([0-9]|[ .-])*\\)",AxonData)
    DataLines<-grep("^[[:space:]]*[(][[:space:]]*([0-9.\\-]+|ERR|NA)+",AxonData)
    if (length(DataLines)<2){
	stop(paste("There must be at least two points to define an axon in",AscFile))
    }
    #cat("There are",length(DataLines),"Data lines")
    # Read them in
    tb<-sapply(strsplit(AxonData[DataLines],split="[ ()]+"),function(x){x[2:5]})
    # Convert them to numeric
    allpoints<-matrix(as.numeric(tb),,4,byrow=T)
    #cat("There are",nrow(allpoints),"points")
    
    # Store two useful lengths
    NumPoints<-nrow(allpoints)
    NumLines<-length(AxonData)
    
    # Make a dataframe to contain the raw point data
    # note that Label is set to 2 (for Axon) - the SWC standard
    # and Parent is set to 0 for the present 
    d<-data.frame(PointNo=1:NumPoints,Label=rep(2,NumPoints),X=allpoints[,1],Y=allpoints[,2],Z=allpoints[,3],W=allpoints[,4],Parent=rep(0,NumPoints))    

    # OK Big old loop to process segments
    OpenNewBranches<-0
    CurIdx<-0
    NextIdx<-1
    RootPoint<-1 # NB always 1 since in Neurolucida we started from the root
 
    while(CurIdx < length(InstructionLines)){
	if(CurIdx==0){
	    # This is the root segment i.e. CurInstruct is the first one after the root
	    d$Parent[1]<- -1
	    SegNo<-1
	    SegList<-list()
	    SegList[[1]]<-which(1:InstructionLines[NextIdx]%in%DataLines)
	    EndPoint<-max(SegList[[1]])
	    #print(1,EndPoint)
	    # Set the parents of these points
	    d$Parent[2:(EndPoint)]<-RootPoint:(EndPoint-1)
	    
	    if((EndPoint+1)<=NumPoints){
		#d$Parent[EndPoint+1]<-EndPoint
	    }
	    
	    # Add the root to the list keeping track of the endpoints
	    #EndPoints<-c(EndPoints,1)
	    # Increment the loop vars
	    CurIdx<-CurIdx+1
	    NextIdx<-1+NextIdx
	    next
	}
	# Check to see if we are at the end of a subtree
	if(InstructionLines[CurIdx]%in%SubtreeEndLines) {
	    OpenNewBranches<-OpenNewBranches[-1]

	    CurIdx<-CurIdx+1
	    NextIdx<-1+NextIdx
	    next
	    
	}
	
	# If not we are at a new branch or a split
	if(InstructionLines[CurIdx]%in%NewBranchLines){
	    # OK This is a new branch from StartPoint...EndPoint
	    SegNo<-SegNo+1
	    #StartPoint<-EndPoint+1
	    StartPoint<-EndPoint
	    RangeOfGoodDataLines<-which(InstructionLines[CurIdx]:InstructionLines[NextIdx]%in%DataLines)
	    EndPoint<-StartPoint+length(RangeOfGoodDataLines)
	    SegList[[SegNo]]<-StartPoint:EndPoint
	    # Set the parents of these points
	    # HAH - but what if the StartPoint and EndPoint are the same!
	    # This was a HORRID bug!
	    if(StartPoint!=EndPoint){
		d$Parent[(StartPoint+1):EndPoint]<-StartPoint:(EndPoint-1)
	    }	    
	    if((EndPoint+1)<=NumPoints){
		#d$Parent[EndPoint+1]<-EndPoint
	    }    
	    
	    # Push this new branch on the top of the stack
	    OpenNewBranches<-c(StartPoint,OpenNewBranches)
	    # Add the branchpoint to the list keeping track of the branchpoints
	    #BranchPoints<-c(BranchPoints,StartPoint)
	    # Increment the loop vars
	    CurIdx<-CurIdx+1
	    NextIdx<-1+NextIdx
	    next
	}
	if(InstructionLines[CurIdx]%in%SplitLines){
	    # OK This is a split, so this seg consists of BranchPoint, StartPoint...EndPoint
	    SegNo<-SegNo+1
	    StartPoint<-EndPoint+1
	    RangeOfGoodDataLines<-which(InstructionLines[CurIdx]:InstructionLines[NextIdx]%in%DataLines)
	    # NB the -1 here is because we already have a separate BranchPoint for splits
	    EndPoint<-StartPoint+length(RangeOfGoodDataLines)-1
	    BranchPoint<-OpenNewBranches[1]
	    SegList[[SegNo]]<-c(BranchPoint,StartPoint:EndPoint)
	    # Set the parents of these points
	    d$Parent[StartPoint]<-BranchPoint
	    # HAH - but what if the StartPoint and EndPoint are the same!
	    # This was a HORRID bug! Using : instead of seq() was 
	    # probably also bad
	    if(StartPoint!=EndPoint){
		d$Parent[(StartPoint+1):EndPoint]<-StartPoint:(EndPoint-1)
	    }
	    
	    if((EndPoint+1)<=NumPoints){
		#d$Parent[EndPoint+1]<-EndPoint
	    }    
	    
	    CurIdx<-CurIdx+1
	    NextIdx<-1+NextIdx
	    next
	}
	# If we've got to here either this is the last point or there's a problem
	
    }  # End of the big old while loop
    
    # OK Figure out the branch points & End Points
    # Copied from ParseSWCTree in NeuronFunctions4.s
    BranchPoints<-as.numeric(names(which(table(d$Parent)>1)))
    #Find out which PointNos occur in $Parent column
    NotEndPoints<-as.numeric(names(table(d$Parent)))
    #Get rid of any -1s
    NotEndPoints<-NotEndPoints[NotEndPoints>0]
    EndPoints<-c(RootPoint,d$PointNo[-NotEndPoints])

    
    # Check if we ended up with something sensible
    if(length(SegList)>0){
	#OK There's at least one segment
	ParsedNeuron<-list(NeuronName=NeuronNameFromFileName(AscFile),
	    InputFileName=AscFile,
	    CreatedAt=Sys.time(),
	    NodeName=Sys.info()["nodename"],
	    InputFileStat=file.info(AscFile)[1,],
	    InputFileMD5=md5sum(path.expand(AscFile)),
	    NumPoints=NumPoints,
	    StartPoint=RootPoint, # NB always 1 since coming Neurolucida
	    BranchPoints=BranchPoints,
	    EndPoints=EndPoints,
	    NumSegs=length(SegList),
	    SegList=SegList,
	    d=d	)
	return(ParsedNeuron)
    }
    else return(-1)
}

WriteAscFromNeuron<-function(ANeuron,filename=NULL,suffix="asc",Force=F,MakeDir=T,width=8,dp=3){
    # writes a neuroluida .asc representation of an asc object

    if(is.null(filename)) filename=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
    if(!Force && file.exists(filename) ){
	warning(paste(filename,"already exists; use Force=T to overwrite"))
	return()
    }
    if(!file.exists(dirname(filename))){
	# either bail
	if(!MakeDir){
	    warning(paste(dirname(filename),"does not exist; use MakeDir=T to overwrite"))
	    return()
	} else {
	    # or try to make a directory
	    if(!dir.create(dirname(filename))){
		warning(paste("Unable to create",dirname(filename)))
	    }
	}
    }
    if(!file.create(filename)){
	warning(paste("Unable to write to file",filename))
	return()
    }
    
    sprintfString=paste("(",  paste(rep(paste("%",width+1,".",dp,"f",sep=""),4)  ,collapse=""),")",sep="")
    OutFile=file(filename,"w")
    
    if(is.null(ANeuron$SegOrders)) ANeuron=SegOrders(ANeuron)
    cat("; V3 text file written for MicroBrightField products\n",file=OutFile)
    cat("; Written by WriteAscFromNeuron at",date(),"\n",file=OutFile)
    cat(";",ANeuron$NeuronName,ANeuron$CellType,
	"NumPoints =",ANeuron$NumPoints,"NumSegs =",ANeuron$NumSegs,"\n\n",file=OutFile)
    
    cat("( (Color Blue)  ; [10,21]\n  (Axon)\n",file=OutFile)
    # this will store how many branches coming off each branch we have
    # written
    NumBranchesWritten=rep(0,length(ANeuron$BranchPoints))
    # The head node of each segment
    SegHeads=sapply(ANeuron$SegList,function(x) x[1])
    SegTails=sapply(ANeuron$SegList,function(x) x[length(x)])
    # the number of dependent branches to write for each node.
    NumBranchesToWrite=table(SegHeads[SegHeads%in%ANeuron$BranchPoints])
    
    # For each seg calculate whether it is the 1st, 2nd etc sub-branch from
    # its parent segment
    SubBranchOrders=unsplit(sapply(split(SegHeads,factor(SegHeads)),seq),factor(SegHeads))
    # Initialise the SegDescription list with the root description
    SegDescriptions="R"
    # Now figure out descriptions for the remaining segments
    for(i in 2:length(ANeuron$SegList)){
	# note that because segments are nicely ordered, each segment
	# description can be built by simply appending to the description
 	# of its parent segment
 	NewSegDescription=paste(SegDescriptions[SegHeads[i]==SegTails],"-",SubBranchOrders[i],sep="")
	SegDescriptions=c(SegDescriptions,NewSegDescription)
    }
  
    for (i in 1:ANeuron$NumSegs){
	# How deep are we?
	IndentString=paste(rep("  ",ANeuron$SegOrders[i]),collapse="")
	#cat(IndentString,"seg order = ",ANeuron$SegOrders[i],sep="")
	oldLevel=ANeuron$SegOrders[SegTails==SegHeads[i]]
	headPoint=ANeuron$SegList[[i]][1]
	
	# 1) OPEN SEGMENT
	# Is this a) root, b) first branch c) another branch from an existing 
	# branch point?
	if(headPoint%in%ANeuron$BranchPoints){
	    if(NumBranchesWritten[which(ANeuron$BranchPoints==headPoint)]==0){
		# b first branch
		cat(rep("  ",oldLevel),"(\n",sep="",file=OutFile)
	    } else {
	       # c cont branch
		cat(rep("  ",oldLevel),"|\n",sep="",file=OutFile)
	    }
	}
	
	# 2) WRITE SEGMENT
	#cat(paste(IndentString,"(",do.call("paste",ANeuron$d[ANeuron$SegList[[i]],c("X","Y","Z","W")]),")",sep="",coll="\n"),
	#    sep="",file=filename)
	# if this is the root - print point 1 of this segment 
	# otherwise do not include the branch point 
	if(i==1) {
	    LinesToPrint=ANeuron$SegList[[i]]
	    Comments=c(" ; Root",paste("  ; ",SegDescriptions[i],", ",seq(length(LinesToPrint)-1),sep=""))
	} else {
	    LinesToPrint=ANeuron$SegList[[i]][-1]
	    Comments=paste("  ; ",SegDescriptions[i],", ",seq(length(LinesToPrint)),sep="")
	}

	FormattedNumbers=apply(  ANeuron$d[LinesToPrint,c("X","Y","Z","W")],1,function(x) 
	    paste(IndentString,sprintf(sprintfString,x[1],x[2],x[3],x[4]),sep="" )  )
	if(length(FormattedNumbers)!=length(Comments)) stop("Number of rows and comments do not match!")
	cat(paste(FormattedNumbers,Comments,collapse="\n",sep=""),sep="\n",file=OutFile)
	
	# figure out the index of this branch point in ANeuron$BranchPoints
	branchPointIdx=which(SegHeads[i]==ANeuron$BranchPoints)
	#cat("i =",i,"branchPointIdx =",branchPointIdx,"\n")
	if(length(branchPointIdx)>0)
	    NumBranchesWritten[branchPointIdx]=NumBranchesWritten[branchPointIdx]+1
	
	# 3) CLOSE SEGMENT
	# Are we at an end point
	if(SegTails[i]%in%ANeuron$EndPoints){
	    # note than NL3 moves Normal one space in from cur indent
	    cat(IndentString," Normal\n",sep="",file=OutFile)
	    # Now check to see if we should be closing split
	    orderOfNextBranch=ifelse(i==ANeuron$NumSegs,1,ANeuron$SegOrders[i+1])
	    if(orderOfNextBranch<ANeuron$SegOrders[i]){
		# nb +1 because the next branch will be open so we
		# don't want to close it
		# nb -1 because these things are normally printed back one
		# indent from the level of the current branch
		for(j in (ANeuron$SegOrders[i]:(orderOfNextBranch+1))-1){
		    IndentString=rep("  ",j)
		    cat(IndentString,")  ;  End of split\n",sep="",file=OutFile)
		}
		if (i == ANeuron$NumSegs) cat(")  ;  End of tree\n",sep="",file=OutFile)
	    }
	} # end of handling closing segments
    } # end of loop over segments
    close(OutFile)
}

    
    
    
