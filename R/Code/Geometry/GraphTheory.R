# GraphTheory.R
# Some functions to test out ideas of graph theory
# on axon trees - used by some Neuron reading routines

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

# Started this after listening to Hermann Cuntz' talk
# at the Hausser lab meeting

# Looks like there is some interesting stuff in the ggm package
# according to their definitions, the adjacency matrix has
# (i,j)=1 if there is a connection from i to j
# while the Edge Matrix has
# (j,i)=1 if there is a connection from i to j
# and also has 1s along the diagonal

# the igraph package looks like it may be a more general purpose graph
# package

require(igraph)

# Adjacency matrix
AdjacencyMatrixFromSegList<-function(SegList,Undirected=FALSE){
	#ps=sort(unique(unlist(SegList)))
	#nps=length(ps)
	nps=max(unlist(SegList))
	if(is.na(nps) || nps<1) stop("No valid points in passed SegList")
	A=matrix(0,ncol=nps,nrow=nps)
	for (i in seq(len=length(SegList))){
		s=SegList[[i]]
		for(j in seq(len=length(s)-1)){
			# Have now adopted the defintion of ggm 
			# ie for adjacency matrix
			# A(i,j)=1 if i->j
			A[s[j],s[j+1]]=1
			if(Undirected) A[s[j+1],s[j]]=1
			#A[s[j],s[j]]=1
		}			
		#c=cbind(s[-length(s)],s[-1])
		#apply(c,1,function(x) A[x[2],x[1]]=1)
	}
	A
}

ReducedAdjacencyMatrixFromSegList<-function(SegList,Undirected=FALSE){
	# omit everything but the head and tail of each segment
	nps=max(unlist(SegList))
	A=matrix(0,ncol=nps,nrow=nps)
	nSegs=length(SegList)
	for (i in seq(len=nSegs)){
		s=SegList[[i]]
		i=s[1];j=s[length(s)]
		# Have now adopted the defintion of ggm 
		# ie for adjacency matrix
		# A(i,j)=1 if i->j
		A[i,j]=1
		if(Undirected) A[j,i]=1
		#A[s[j],s[j]]=1
	}	
	rownames(A)=1:nrow(A)
	colnames(A)=1:nrow(A)
	toKeep=which(rowSums(A)!=0 | colSums(A)!=0)
	A[toKeep,toKeep]
}

RerootNeuron<-function(ANeuron,origin=1,...){
  gam=Neuron2Graph(ANeuron)
  # will use a depth first search to reorder tree starting from origin
  dfs=graph.dfs(gam,origin,father=TRUE)
  coreneuron=SegListFromFullGraph(gam,dfs=dfs)
  ANeuron[names(coreneuron)]<-coreneuron
  # Now fix the SWC data chunk as well
  d=ANeuron$d
  d$Parent=dfs$father
  # SWC says that the root will have parent -1
  d$Parent[d$Parent==0]=-1L
  
  ANeuron$d=d
  ANeuron
}

#' Construct SegList (+ other core fields) from graph of all nodes and origin
#'
#' Uses a depth first search on the tree to reorder using the given origin.
#' @details Either one of origin and dfs must be specified.
#' @param g An igraph
#' @param origin the 1-indexed root vertex id
#' @param dfs result of graph.dfs
#' @return A list with elements:
#'  (NumPoints,StartPoint,BranchPoints,EndPoints,NumSegs,SegList)
#' @export
#' @seealso \code{\link{graph.dfs},\link{RerootNeuron}}
SegListFromFullGraph<-function(g,origin=NULL,dfs=NULL){
  if(is.null(dfs)){
    if(is.null(origin)) stop("Must specify at least one of dfs or origin")
    if(origin<1 || origin>vcount(g)) stop("invalid origin:",origin)
    dfs=graph.dfs(g,root=origin,father=TRUE)
  } else {
    # dfs seems to have origin 0-indexed (but everything else 1-indexed)
    if(is.null(origin)) origin=dfs$root+1
    else {
      # check that the origin we were given matches dfs
      stop("origin specified in dfs: ",dfs$root+1,
        ' does not match that on command line: ',origin)
    }
  }
  ncount=neighborhood.size(g,order=1)-1
  # put the first vertex into the first segment
  curseg=dfs$order[1]
  if(length(ncount)==1) sl=list(curseg)
  else {
    sl=list()
    # we have more than 1 point in graph and some work to do!
    for(i in seq.int(from=2,to=length(dfs$order))){
      curpoint=dfs$order[i]
      if(length(curseg)==0){
        # segment start, so initialise with parent
        curseg=dfs$father[curpoint]
      }
      # always add current point
      curseg=c(curseg,curpoint)
      # now check if we need to close the segment
      if(ncount[curpoint]!=2){
        # branch or end point
        sl[[length(sl)+1]]=curseg
        curseg=integer(0)
      }
    }
  }
  list(NumPoints=length(ncount),
  StartPoint=origin,
  BranchPoints=seq.int(length.out=length(ncount))[ncount>2],
  EndPoints=seq.int(length.out=length(ncount))[ncount==1],
  NumSegs=length(sl),
  SegList=sl)
}

AdjacencyMatrixFromEdgeList<-function(EdgeList,Undirected=FALSE){
	# this makes an adjacency matrix from an edge list
	# ie a Nx2 matrix where each row defines an edge
	# if it is a double edge list it will be undirected
	
	# the cols of input should be in the order from -> to
	nps=max(EdgeList)
	A=matrix(0,ncol=nps,nrow=nps)
	apply(EdgeList,1,function(x) A[x[1],x[2]]<<-1)
	A
}

ReducedAdjacencyMatrixFromEdgeList<-function(EdgeList){
	A=AdjacencyMatrixFromEdgeList(EdgeList)
	rownames(A)=1:nrow(A)
	colnames(A)=1:nrow(A)
	toKeep=which(rowSums(A)!=2 | colSums(A)!=2)
	A[toKeep,toKeep]
}

		
AdjacencyMatrix<-function(SWCData){
	# For a set of data in SWC format
	# calculate the adjacency matrix
	# This looks like
	#       Source
	#    ---------------------
	#    | 0 0 0 0
	# T  | 1
	# a  | 0
	# r  | 0
	# g  | 0
	# e  | 0
	# t  | 0
	#    | 0
	#    | 0

	
	# First question - should this be directed or not?
	# Hermann presented it as a directed routine
	# OK for SWC data
	# just iterate over points considering target only
	# for each target look up parent and put a 1 in that column

	n=nrow(SWCData)
	A=matrix(0,ncol=n,nrow=n)
	Parents=subset(SWCData,Parent>0,sel=c(PointNo,Parent))
	for(i in 1:nrow(Parents)){
		A[Parents$PointNo[i],Parents$Parent[i]]=1
	}
	A
}

# For a simple Y
#        3
# 1 -> 2 
#        4
#
colSwap=function(A,a,b){
	if(a>ncol(A) || b>ncol(A)) stop("a or b exceeds matrix dimensions")
	tcol=A[,a]
	A[,a]=A[,b]
	A[,b]=tcol
	A
}

NumTreesInEdgeList<-function(Nb){
	# returns the number of trees in a set of neighbour pairs
	nPoints=length(unique(Nb[,2]))
	nEdges=nrow(Nb)
	nPoints-nEdges/2 #nTrees
}

FindSubTreesFromEdgeList<-function(Nb){
	# if the edge list contains more than one tree
	# return a list containing each subtree
	nTrees=NumTreesInEdgeList(Nb)
	# This would look different since 
	# when nTrees>1 the result would be 
	# just the unique point numbers
	#if(nTrees==1) return(Nb)
	
	NbTrees=list()
	j=1
#	for (j in seq(nTrees)){
	while(any(Nb$CurPoint>-1)){
		cat("j =",j,"\n")
		PossPoints=unique(Nb$CurPoint[Nb$CurPoint>-1])
		NbTrees[[j]]=integer(length(PossPoints))
		NbTrees[[j]][1]=min(PossPoints)
		Nb[Nb$CurPoint==NbTrees[[j]][1],]=-1
		nThisTree=1
		while(T){
			#Find neighbours of current points
			#NewNeighbours=setdiff(Nb$CurPoint[ Nb$Neighbour%in%NbTrees[[j]] ],NbTrees[[j]])
			# Nb unique is required because occasionally we will pick
			# up a new neighbour common to two points in the same round 
			NewNeighbours=unique(Nb$CurPoint[ Nb$Neighbour%in%NbTrees[[j]][1:nThisTree] ])
			#cat("NN=",NewNeighbours)
			nNewNeighbours=length(NewNeighbours)
			if(nNewNeighbours==0) break
			# set the next few points of the result vector
			# ... to the NewNeighbours
			NbTrees[[j]][nThisTree+(1:nNewNeighbours)]=NewNeighbours
			nThisTree=nThisTree+nNewNeighbours
			Nb[Nb$CurPoint%in%NewNeighbours,]=-1
		}
		NbTrees[[j]]=NbTrees[[j]][NbTrees[[j]]>0]
		j=j+1
	}
	return(NbTrees)
}

# FindRootedSubTree<-function(Nb,root){
# 	Points=unique(Nb$CurPoint)
# 	BranchPoints=as.numeric(names(which(table(Nb$CurPoint)>2)))
# 	VisitedBranchPoints=NULL
# 	EndPoints=as.numeric(names(which(table(Nb$CurPoint)==1)))
# 	
# 	SegList=list()
# 	parseSeg=function(headPoint){
# 		curParent=headPoint
# 		splitPoint=0
# 		curPoint=curParent
# 		while(nrow(Nb)>0){
# 			possPoints=Nb$Neighbour[Nb$CurPoint==curPoint]
# 			# Pick the smallest point that is >0 (will use -1)
# 			nextPoint=min(PossPoints[PossPoints>0])
# 			Nb[Nb$CurPoint==nextPoint,]=-1
# 			if(any(BranchPoints==thisPoint)){
# 				#is this a new branch point?
# 				if(all(VisitedBranchPoints!=nextPoint)){
# 					VisitedBranchPoints=c(VisitedBranchPoints,nextPoint)
# 					splitPoint=curParent
# 					SegList=list(SegList,curParent)
# 					curPoint=nextPoint
# 					# parse the first seg
# 					rval=parseSeg(curParent)
# 					while(rval!=0){
# 						# parse subsequent segs
# 						rval=parseSeg(splitPoint)
# 					}
# 				} else {
# 					# a split
# 					return(curParent)
# 				}
# 			}  # not a fork
# 			else if(any(EndPoints==nextPoint)){
# 				return(0)
# 		}
# 		
# 			
	
	

SWCFromAdjacencyAndVertices<-function(A,d){
	# A has cols and rows corresponding to points in d
	colnames(A)=seq(ncol(A))
	rownames(A)=seq(nrow(A))
	# leftmost col
	# Pick first element
	# cross off that element
	# and its reciprocal
	# Branch
}

# A=rbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),c(0,1,0,0))


AmiraDataFromSWC<-function(d){
	el=subset(d,Parent!=-1,sel=c(Parent,PointNo))
	el=DoubleFromSingleEdgeList(el)
	Origin=subset(d,Parent==-1)$CurPoint
	list(PointList=d,EdgeList=el,Origin=Origin)
}

EdgeListFromNeuron<-function(n){
	EdgeListFromSWC(n$d)
}

EdgeListFromSWC<-function(d){
	el=subset(d,Parent!=-1,sel=c(Parent,PointNo))
	DoubleFromSingleEdgeList(el)
}

CoreNeuronFromAmiraData<-function(l){
	
}

DoubleFromSingleEdgeList<-function(el){
	el=matrix(unlist(el),ncol=2)
#	as.matrix(el)
	reversed=t(apply(el,1,rev))
	#interleave
	el=rbind(el,reversed)
	names(el)=c("CurPoint","Neighbour")
	el[order(el[,1],el[,2]),]
}

SingleFromDoubleEdgeList<-function(el){
	#Hmm is this general
	apply(unique(t(apply(el,1,function(x) paste(sort(x))))),2,as.numeric)
}

Neuron2Graph<-function(x){
	.Deprecated("as.igraph.neuron",
	  'Deprecated, due to use of dense adjacency matrix')
	# returns an igraph graph object
	require(igraph)
	am=AdjacencyMatrixFromSegList(x$SegList)
	graph.adjacency(am,mode='undirected')
}

GetShortestPath.Neuron<-function(x,from,to){
	# nb igraph is 0 indexed
	g=Neuron2Graph(x)
	p=get.shortest.paths(g,from=from-1,to=to-1)
	if(length(p)!=1) stop("Unable to find unique shortest path between those points")
	p[[1]]+1
}

GetShortestPath.Neuron.ggm<-function(x,from,to){
	# just for checking
	require(ggm)
	p=findPath(am,from,to)
}

GetSubNeuron<-function(x,seglist=NULL,from,to){
	# Very simple function to extract a sub neuron from an original neuron
	# basically leaves everything intact but changes the segment list
	if(missing(to) && missing(from)){
		if(is.null(seglist)) stop("Must supply a seglist")
		if(!is.list(seglist)) 
			seglist<-list(seglist)		
	} else{
		if(missing(to)){
			to=from
			from=seglist
			seglist<-list(GetShortestPath.Neuron(x,from,to))
		}
	}
	if(!all(unlist(seglist)%in%x$d$PointNo)) stop("Seglist addresses points outside of neuron")
	xx=x
	xx$SegList=seglist
	# reset the start point
	if(!(xx$StartPoint%in%seglist)) xx$StartPoint=seglist[[1]][1]
	xx
}
