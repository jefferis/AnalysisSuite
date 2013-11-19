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

require(igraph)

# Adjacency matrix
AdjacencyMatrixFromSegList<-function(SegList,Undirected=FALSE){
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
		}
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

#' Reroot a neuron on the given origin using depth first search reordering
#'
#' Uses a depth first search on the tree to reorder using the given origin.
#' NB this does _not_ reorder the vertex ids - it changes the origin and then
#' changes the parents in SWC data block and SegList accordingly.
#' @details Either one of origin and dfs must be specified.
#' @param ANeuron Neuron to be rerooted
#' @param origin the 1-indexed root vertex id
#' @param graph.method See method argument of \code{\link{as.igraph.neuron}}
#' @return A list with elements:
#'  (NumPoints,StartPoint,BranchPoints,EndPoints,NumSegs,SegList)
#' @export
#' @seealso \code{\link{graph.dfs},\link{as.igraph.neuron}}
RerootNeuron<-function(ANeuron,origin=1,graph.method=c("swc",'seglist')){
  gam=as.igraph(ANeuron,method=graph.method)
  # will use a depth first search to reorder tree starting from origin
  dfs=graph.dfs(gam,root=origin,father=TRUE,neimode='all')
  coreneuron=CoreNeuronFromGraph(gam,dfs=dfs)
  ANeuron[names(coreneuron)]<-coreneuron
  # Now fix the SWC data chunk as well
  ANeuron$d$Parent=dfs$father
  # SWC says that the root will have parent -1
  ANeuron$d$Parent[ANeuron$d$Parent==0]=-1L
  ANeuron
}

#' Construct EdgeList matrix with start and end points from SegList
#'
#' @param SegList from a \code{neuron}
#' @return A 2 column matrix, \code{cbind(starts,ends)}
#' @export
EdgeListFromSegList<-function(SegList){
  lsl=sapply(SegList,length)
  sl=SegList[lsl>1]
  lsl=lsl[lsl>1]
  ends=unlist(lapply(sl,function(x) x[-1]))
  starts=unlist(lapply(sl,function(x) x[-length(x)]))
  cbind(starts,ends)
}

#' Construct an igraph object from a neuron
#'
#' @details note that the 'swc' and 'seglist' methods may generate different
#'  graphs, but these should always be isomorphic and this can be checked by
#'  using canonical.permutation()
#' @param x The neuron to be converted
#' @param directed Whether to produce a directed graph (default TRUE)
#' @param method Whether to use the seglist or swc data blocks to generate
#'  the graph.
#' @param prune Whether to prune the graph of any vertices not contained in
#'  the input graph (default TRUE)
#' @param keep.ids Whether to keep the original numeric ids for each vertex
#' @param ... Additional arguments, currently ignored
#' @return A matrix, \code{cbind(starts,ends)}
#' @method as.igraph neuron
#' @rdname as.igraph
#' @export
as.igraph.neuron<-function(x,directed=TRUE,method=c("swc",'seglist'),
  prune=TRUE, keep.ids=prune, ...){
  method=match.arg(method,several.ok=TRUE)
  if('swc'%in%method && !is.null(x$d$Parent) && !is.null(x$d$PointNo)){
    as.igraph.swc(x$d, directed=directed, prune=prune, keep.ids=keep.ids)
  } else {
    as.igraph.seglist(x$SegList, directed=directed, prune=prune, keep.ids=keep.ids)
  }
}

#' @rdname as.igraph
#' @export
as.igraph.seglist<-function(x, directed=TRUE, prune=TRUE, keep.ids=prune, ...){
  el=EdgeListFromSegList(x)
  as.igraph.edgelist(el, directed=directed, prune=prune, keep.ids=keep.ids, ...)
}

#' @rdname as.igraph
#' @export
as.igraph.edgelist<-function(x, directed=TRUE, prune=TRUE, keep.ids=prune, ...){
  g=graph.edgelist(x,directed=directed)
  if(keep.ids) g=set.vertex.attribute(g, 'label', value=igraph::V(g))
  if(prune) g=delete.vertices(g, setdiff(igraph::V(g), unique(as.vector(unlist(x)))))
  g
}

#' @rdname as.igraph
#' @export
as.igraph.swc<-function(x, directed=TRUE, prune=TRUE, keep.ids=prune, ...){
  el=data.matrix(EdgeListFromSWC(x))
  as.igraph.edgelist(el, directed=directed, prune=prune, keep.ids=keep.ids, ...)
}

#' Make core neuron elements from a block of SWC data
#'
#' @param swc Matrix or dataframe of swc format data
#' @rdname CoreNeuron
#' @seealso \code{\link{CoreNeuronFromGraph}}
CoreNeuronFromSWC<-function(swc,origin=NULL){
  g=as.igraph.swc(swc,directed=TRUE)
  CoreNeuronFromGraph(g)
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
#' @rdname CoreNeuron
#' @seealso \code{\link{graph.dfs},\link{RerootNeuron}}
CoreNeuronFromGraph<-function(g,origin=NULL,dfs=NULL){
  if(is.null(dfs)){
    if(is.null(origin)) {
      # if we haven't specified an origin, then we can infer this for a
      # directed graph.
      if(is.directed(g)){
        origin=rootpoints(g)
        if(length(origin)>1){
          warning("Graph has multiple origins. Using first")
          origin=origin[1]
        }
      }
      else stop("Must specify at least one of dfs or origin")
    }
    if(origin<1 || origin>vcount(g)) stop("invalid origin:",origin)
    dfs=graph.dfs(g,root=origin,father=TRUE,neimode='all')
  } else {
    # dfs seems to have origin 0-indexed (but everything else 1-indexed)
    if(is.null(origin)) origin=dfs$root+1
    else {
      # check that the origin we were given matches dfs
      stop("origin specified in dfs: ",dfs$root+1,
        ' does not match that on command line: ',origin)
    }
  }
  ncount=igraph::degree(g)
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

#' Return the root or branch points of a neuron or graph
#'
#' A neuron may have multiple subtrees and therefore multiple roots
#' @rdname rootpoints
rootpoints<-function (x, ...)
UseMethod("rootpoints")

#' @rdname rootpoints
rootpoints.neuron<-function(x, ...){
  if(x$nTrees>1) sapply(x$SubTrees, function(y) y[[1]][1])
  else x$StartPoint
}

#' @rdname rootpoints
rootpoints.igraph<-function(x, original.ids=TRUE, ...){
  if(!is.directed(x))
    stop("Cannot establish root points for undirected graph")

  # root points are those without incoming edges
  vertex_ids=igraph::V(x)[igraph::degree(x,mode='in')==0]
  if(original.ids)
    vertex_names=get.vertex.attribute(x,'label',index=vertex_ids)
  if(original.ids || is.null(vertex_names))
    as.integer(vertex_ids)
  else
    vertex_names
}

#' Return the branchpoints of a neuron or graph
#' @param x neuron or graph
#' @param ... Additional parameters
#' @export
#' @rdname rootpoints
#' @alias branchpoints
branchpoints<-function (x, ...)
UseMethod("branchpoints")

#' @rdname rootpoints
#' @detail returns a list if more than one subtree is specified
branchpoints.neuron<-function(x, subtrees=1, ...){
  if(isTRUE(subtrees==1)) x$BranchPoints
  else if(any(subtrees>x$nTrees)) stop("neuron only has ",x$nTrees," subtrees")
  else lapply(x$SubTrees[subtrees],
    function(x) branchpoints(as.igraph.seglist(x)))
}

#' @rdname rootpoints
#' @param original.ids Whether to return original point ids when available
branchpoints.igraph<-function(x, original.ids=TRUE, ...){
  vertex_ids=igraph::V(x)[igraph::degree(x)>2]
  if(original.ids)
    vertex_names=get.vertex.attribute(x,'label',index=vertex_ids)
  if(original.ids || is.null(vertex_names))
    as.integer(vertex_ids)
  else vertex_names
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

AmiraDataFromSWC<-function(d){
	el=DoubleEdgeListFromSWC(d)
	Origin=subset(d,Parent==-1)$CurPoint
	list(PointList=d,EdgeList=el,Origin=Origin)
}

EdgeListFromNeuron<-function(n){
	EdgeListFromSWC(n$d)
}

DoubleEdgeListFromNeuron<-function(n){
	DoubleFromSingleEdgeList(EdgeListFromNeuron(n))
}

EdgeListFromSWC<-function(d){
	subset(d,Parent!=-1,sel=c(Parent,PointNo))
}

DoubleEdgeListFromSWC<-function(d){
  el=EdgeListFromSWC(d)
  DoubleFromSingleEdgeList(el)
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
	# nb igraph is 1-indexed
	g=as.igraph(x)
	p=get.shortest.paths(g,from=from,to=to)
	if(length(p)!=1) stop("Unable to find unique shortest path between those points")
	p[[1]]
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
