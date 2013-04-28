# More generic code acting on neurons that may be suitable for inclusion in 
# a future neuron package

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

#' Recalculate Neurons's SWCData using SegList and point information
#'
#' Uses the SegList field (indices into point array) to recalculate 
#' point numbers and parent points for SWC data field (d).
#' If Label column is missing (or ReplaceLabel=TRUE) then it is set to the
#' value of DefaultLabel (2=Axon by default).
#' Note that the order of point indices in SegList must match those in SWC.
#' @param Neuron Must contain both the SegList and d fields
#' @param RecalculateParents Whether to recalculate parent points (default T)
#' @param DefaultLabel Integer label to use for SWC data chunk
#' @param ReplaceLabel Whether to replace Label column if it already exists
#' @return return value
#' @export
#' @seealso \code{\link{ParseSWC}}
#' @examples
RecalculateSWCData<-function(Neuron,RecalculateParents=TRUE,DefaultLabel=2L,ReplaceLabel=FALSE){
  sl=Neuron$SegList
  d=Neuron$d
  if(is.null(d$PointNo))
    d$PointNo=seq(nrow(d))
  else {
    # check that points are consecutive
    if(any(d$PointNo!=seq(nrow(d))))
      stop("Points must be numbered consecutively from 1:npoints")
  }
  # NB this assumes that points are in order
  if(is.null(d$Parent) || RecalculateParents){
    d$Parent=-1L
    for(s in sl){
      d$Parent[s[-1]]=s[-length(s)]
    }
  }
  if(is.null(d$Label) || ReplaceLabel)
    d$Label=DefaultLabel
  Neuron$d=d
  Neuron
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

MergeNeuron<-function(src,dest,dryrun=FALSE,verbose=T){
	# merge additional fields of src neuron into dest
	missingFields=setdiff(names(src),names(dest))
	if(length(missingFields)>0){
		if(verbose) cat(ifelse(dryrun,"Dry Run",""),"Merging fields",missingFields,"for",dest$NeuronName,"\n")
		if(!dryrun){
			oldClass=class(dest)
			dest=c(dest,src[missingFields])
			class(dest)=oldClass
		}
	} else if(verbose) cat("No fields to merge")
	invisible(dest)
}
