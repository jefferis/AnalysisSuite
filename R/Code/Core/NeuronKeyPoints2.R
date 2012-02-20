# NeuronKeyPoints.R
# The functions in this file should allow 
# the key points of the neuron to be picked out 
# by the user.
# An additional function then segments the neuron
# into Axon (1), LH (2), MB (3), other (-1)
# These key points should then enable marking up of the neuron
# so that the axon, MB and LH regions are identified

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

# Depends on functions in NeuronFunctions3.R
# 
# Neuron Long.GetKeyPoints<-function(ANeuron2){
# Neuron GetKeyPoints<-function(ANeuron2){
# Neuron PartitionNeuron<-function(ANeuron){
# RootPoint AskForRoot<-function(ANeuron)
# MBPoints AskForMBPoints<-function(ANeuron)
# LHBranchPoint AskForLHBranchPoint<-function(ANeuron)

# source(file.path(CodeDir,"NeuronKeyPoints2.R"))


# GetKeyPointsofSomeNeurons
# Expects a vector of indices into MyNeurons
# Defaults to all neurons
# ... is to provide additional parameters to plot routines
GetKeyPointsofSomeNeurons<-function(NeuronList=1:length(MyNeurons),...){
    for(ThisNeuron in NeuronList){
	# NB don't forget <<- for global assignment
	MyNeurons[[ThisNeuron]]<<-PartitionNeuron(Short.GetKeyPoints(MyNeurons[[ThisNeuron]],...))
    }
}

# I prefer this one still
Long.GetKeyPoints<-function(ANeuron2){
    RootPoint<-AskForRoot(ANeuron2)
    tmp<-plotneuron2d(ANeuron2)
    title(main="Thank you.  Now calculating")
    ANeuron2<-reroot(ANeuron2,RootPoint)
    ANeuron2<-SegOrders(ANeuron2)

    MBPoints<-AskForMBPoints(ANeuron2)
    LHBranchPoint<-AskForLHBranchPoint(ANeuron2)
    
    if(is.null(ANeuron2$MBPoints)){
	ANeuron2<-c(ANeuron2,list(MBPoints=MBPoints))
    } else {
        ANeuron2$MBPoints<-MBPoints
    }
    if(is.null(ANeuron2$LHBranchPoint)){
	ANeuron2<-c(ANeuron2,list(LHBranchPoint=LHBranchPoint))
    } else {
        ANeuron2$LHBranchPoint<-LHBranchPoint
    }

    return(ANeuron2)
}

# Doesn't ask for Root point
# Assumes it's already been set
# Should do something about increasing the size of the plot window
# since titles are not displayed if they are too long to fit
# Decided to add a ... to pass additional parameters to the plot
# routines
Short.GetKeyPoints<-function(ANeuron2,Components=c("MB","LH"),...){

    if("MB"%in%Components){
	MBPoints<-AskForMBPoints(ANeuron2,...)
	if(is.null(ANeuron2$MBPoints)){
	    ANeuron2<-c(ANeuron2,list(MBPoints=MBPoints))
	} else {
	    ANeuron2$MBPoints<-MBPoints
	}
    }
    
    if("LH"%in%Components){
	LHBranchPoint<-AskForLHBranchPoint(ANeuron2,...)
	if(is.null(ANeuron2$LHBranchPoint)){
	    ANeuron2<-c(ANeuron2,list(LHBranchPoint=LHBranchPoint))
	    print("LHBP Didn't exist")
	} else {
	    ANeuron2$LHBranchPoint<-LHBranchPoint
	    print("LHBP Did exist")
	}
    }
   

    
    return(ANeuron2)
}


GetKeyPoints<-function(ANeuron2){
    RootPoint<-AskForRoot(ANeuron2)
    t<-plotneuron2d(ANeuron2)
    title("Please select the Root, Last MB BranchPoint & 1st LH BranchPoint")
    
    KeyPoints<-t$PointNo[identify(t$X,t$Y,labels=t$PointNo,n=3)]

    t<-plotneuron2d(ANeuron2)
    title(main="Thank you.  Now calculating")
    
    #RootPoint<-which(
    ANeuron2<-reroot(ANeuron2,RootPoint)
    ANeuron2<-SegOrders(ANeuron2)

    MBPoints<-AskForMBPoints(ANeuron2)
    LHBranchPoint<-AskForLHBranchPoint(ANeuron2)
    
    if(is.null(ANeuron2$MBPoints)){
	ANeuron2<-c(ANeuron2,list(MBPoints=MBPoints))
    } else {
        ANeuron2$MBPoints<-MBPoints
    }
    if(is.null(ANeuron2$LHBranchPoint)){
	ANeuron2<-c(ANeuron2,list(LHBranchPoint=LHBranchPoint))
	print("LHBP Didn't exist")
    } else {
        ANeuron2$LHBranchPoint<-LHBranchPoint
	print("LHBP Did exist")
    }

    return(ANeuron2)
}

# Figure out which segments are
# a) Axon (all), AxonPreMB,AxonMB,AxonPostMB
# b) MB Collaterals
# c) LH (Everything after branch point

# This works at the moment - the problem is that is not ready for neurons
# without a MB at all i.e. mACT neurons.

PartitionNeuron<-function(ANeuron,Components=c("Axon","MB","LH"))
{
    # Fill the result arrays with -1 
    PointOrder<-rep(-1,ANeuron$NumPoints)
    SegOrderArray<-rep(-1,ANeuron$NumSegs)
    #Set the tree root to order 0
    PointOrder[ANeuron$StartPoint]<-0
    #Set the first segment to order 1
    SegOrderArray[1]<-1
    #Set the other points in that segment to Order 1
    PointOrder[ANeuron$SegList[[1]][-1]]<-1
    for(i in 2:ANeuron$NumSegs){
	#Find the order of the segment head
	SegHead<-ANeuron$SegList[[i]][1]
	SegHeadOrder<-PointOrder[SegHead]
	#Set the order of the segment to 1 greater
	SegOrderArray[i]<-SegHeadOrder+1
	#Now set the points in that segment as well
	PointOrder[ANeuron$SegList[[i]][-1]]<-SegHeadOrder+1
    } # end of for(i in 2:ANeuron$NumSegs
    
    # Set up an array giving the number of the 1st seg containing any point
    FirstSegOfPoint<-NULL
    for(i in ANeuron$NumSegs:1){
	LastPointOfSeg<-ANeuron$SegList[[i]][length(ANeuron$SegList[[i]])]
	FirstSegOfPoint[ LastPointOfSeg ]<-i
    }
    
    # The above is equivalent to: Really?
    z<-rep(0,ANeuron$NumPoints)
    y<-NULL
    z[y]<-which(sapply(ANeuron$SegList,function(x){x[length(x)]})>0)
    
    #An Array containing the head (i.e. leading branch point)
    # of each segment
    HeadOfSegArray<-NULL
    HeadOfSegArray<-sapply(ANeuron$SegList,function(x){x[1]})
    
    AxonSegNos=NULL;LHSegNos=NULL;MBSegNos=NULL;MBBranchNum=0
    
    
    # OK All axon segments.
    if("Axon"%in%Components){
	# Start at the LH BranchPoint and go to root
	LHBPSegNo<-FirstSegOfPoint[ANeuron$LHBranchPoint]
	AxonSegNos<-NULL
	CurrSeg<-LHBPSegNo
	while(!is.null(CurrSeg)){
	    HeadOfSeg<-ANeuron$SegList[[CurrSeg]][1]
	    AxonSegNos<-c(CurrSeg,AxonSegNos)
	    # For debug!
	    # cat(HeadOfSeg,AxonSegNos," CurrSeg:",CurrSeg,"\n")
	    if(any(ANeuron$EndPoints==HeadOfSeg)){
		#we've got to the end
		break()
	    }
	    CurrSeg<-FirstSegOfPoint[HeadOfSeg]
	}
    }
    
    #Since MBPoints get returned just in the order of their indices
    # rather than in the order they were clicked
    # Need to sort by Order Of MBPs
    if("MB"%in%Components){
	OrderOfMBPs<-sort(PointOrder[ANeuron$MBPoints])
	MBPs<-ANeuron$MBPoints[order(PointOrder[ANeuron$MBPoints])]
    }
    
    # OK find 3 types of axon segment (PreMB, MB, PostMB)
    if(all(c("MB","Axon")%in%Components)){
	PreMBAxons<-AxonSegNos[which(ANeuron$SegOrders[AxonSegNos]<=OrderOfMBPs[1])]
	MBAxons<-AxonSegNos[which(
	    ANeuron$SegOrders[AxonSegNos]<=OrderOfMBPs[2]
	    & ANeuron$SegOrders[AxonSegNos]>OrderOfMBPs[1])]
	PostMBAxons<-AxonSegNos[-(1:length(c(MBAxons,PreMBAxons)))]
    }
    
    #LH Branches
    # Start at LH BranchPoint
    # Keep Picking segs while order>LHBPOrder
    #THIS WORKS
    if("LH"%in%Components){
	LHBPSegNo<-FirstSegOfPoint[ANeuron$LHBranchPoint]
	OrderLHBPSeg<-PointOrder[ANeuron$LHBranchPoint]
	LHSegNos<-NULL
	CurrSeg<-LHBPSegNo+1
	# This works because even if some MB segs have an order > than 
	# some LH seg there will still be an intervening seg
	# in the main SegList with lower order.
	while(ANeuron$SegOrders[CurrSeg]>OrderLHBPSeg & CurrSeg<=ANeuron$NumSegs){
	    LHSegNos<-c(LHSegNos,CurrSeg)
	    CurrSeg<-CurrSeg+1
	}
    }
    
    #OK MB Collaterals now 
    if(all(c("MB","Axon")%in%Components)){
	MBSegNos<-list()
	MBBranchNum<-0
	CurrSeg<-MBPs[1]   # MBPoint closer to root
	#The Branch points to check are the heads of the
	# axon segments determined to be in the MB region
	# + the foot of the last one.
	BPsToCheck<-c(HeadOfSegArray[MBAxons],MBPs[2])
	# Find the segs whose head is one of those branchpoints
	SegsToStartFrom<-unlist(sapply(BPsToCheck,function(x){which(HeadOfSegArray==x)}))
	# Now need to trace from these segments
	for(StartSeg in SegsToStartFrom){
	    # ignoring those segments which are already determined to be axons
	    if (all(AxonSegNos!=StartSeg)){
		# OK This start seg isn't an axon
		MBBranchNum<-MBBranchNum+1
		CurrSeg<-StartSeg
		ThisBranch<-CurrSeg
		CurrSeg<-CurrSeg+1
		# Keep going down the main segment list 
		# until we get to one of the segs to start from, an axon segment or the end
		while(!any(SegsToStartFrom==CurrSeg) & !any(AxonSegNos==CurrSeg) & CurrSeg<=ANeuron$NumSegs){
		    ThisBranch<-c(ThisBranch,CurrSeg)
		    CurrSeg<-CurrSeg+1
		}
		MBSegNos[[MBBranchNum]]<-ThisBranch
	    } # end if (!any(AxonSegNos==StartSeg))
	} # end for(StartSeg in SegsToStartFrom)
    }
    
    SegTypes<-rep(-1,ANeuron$NumSegs)
    if("Axon"%in%Components) SegTypes[AxonSegNos]<-1
    if("LH"%in%Components) SegTypes[LHSegNos]<-2
    if(all(c("MB","Axon")%in%Components)) {
	SegTypes[unlist(MBSegNos)]<-3
	
	#Decided to turn AxonSegNos into a list in the end:
	AxonSegNos<-list(PreMBAxons=PreMBAxons,MBAxons=MBAxons,
	    PostMBAxons=PostMBAxons)
	# unlist(AxonSegNos will return the vector.
    }
    
    #OK Now lets construct the return neuron if required
    if(is.null(ANeuron$AxonSegNos)){
	ANeuron<-c(ANeuron,list(
	    SegTypes=SegTypes,AxonSegNos=AxonSegNos,
	    LHSegNos=LHSegNos,MBSegNos=MBSegNos, 
	    NumMBBranches=MBBranchNum))
    } else {
        ANeuron$SegTypes<-SegTypes
	ANeuron$AxonSegNos<-AxonSegNos
	ANeuron$LHSegNos<-LHSegNos
	ANeuron$MBSegNos<-MBSegNos
	ANeuron$NumMBBranches<-MBBranchNum
    }
    #This 2 point array has the MB Points in the correct order
    #(Remember that thing with identify)
    if("MB"%in%Components)  ANeuron$MBPoints<-MBPs

    return(ANeuron)
}


    

# Function AskForRoot(ANeuron)
# asks the user to identify the root point of the Neuron
# Turns out there is a problem with identify - it doesn't return points
# in the order they were clicked
AskForRoot<-function(ANeuron,...){
    RootPoint<-c(0,0)
    while(length(RootPoint)!=1){
	t<-plotneuron2d(ANeuron,
	    MainTitle=paste(ANeuron$NeuronName,"Click on root, then press return. On error click again"),...)
	RootPoint<-t$PointNo[identify(t[,2],t[,3],labels=t$PointNo,n=2)]
    }
    return(RootPoint)
}

AskForMBPoints<-function(ANeuron,PlotAxes="XZ",...){
    MBPoints<-0
    while(length(MBPoints)!=2){
	#Plot out the neuron
	#Decided that PlotAxes
	t<-plotneuron2d(ANeuron,MainTitle=paste(ANeuron$NeuronName,"Please select the first and last MB points.  If only 1 click root"),PlotAxes=PlotAxes,...)
	MBPoints<-t$PointNo[identify(t[,2],t[,3],labels=t$PointNo,n=2)]
	if(any(MBPoints==ANeuron$StartPoint)) {
	    MBPoints<-rep(MBPoints[which(MBPoints!=ANeuron$StartPoint)],2)
	}
	
    }
    return(MBPoints)
}

AskForLHBranchPoint<-function(ANeuron,PlotAxes="XZ",...){
    LHBranchPoint<-c(0,0)
    while(length(LHBranchPoint)!=1){
	t<-plotneuron2d(ANeuron,
	    MainTitle=paste(ANeuron$NeuronName,"Click on 1st branch of LH"),PlotAxes=PlotAxes,...)
	LHBranchPoint<-t$PointNo[identify(t[,2],t[,3],labels=t$PointNo,n=1)]
    }
    return(LHBranchPoint)
}


AskForLHAnchorPoint<-function(ANeuron,PlotAxes="XZ",...){
	LHAnchorPoint<-c(0,0)
	while(length(LHAnchorPoint)!=1){
	t<-plotneuron2d(ANeuron,
		MainTitle=paste(ANeuron$NeuronName,"Click on LH Anchor Point"),PlotAxes=PlotAxes,...)
	LHAnchorPoint<-t$PointNo[identify(t[,2],t[,3],labels=t$PointNo,n=1)]
	}
	return(LHAnchorPoint)
}

