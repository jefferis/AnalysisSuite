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

#' Return a simplified segment graph for a neuron
#' 
#' @details The resultant graph will contain all branch and endpoints of the
#'   original neuron. This will be constructed from the SegList field, or where
#'   present, the SubTrees field (containing multiple SegLists for each isolated
#'   graph in the neuron). Each edge in the output graph will match one segment
#'   in the original SegList.
#' @param x neuron
#' @param exclude.isolated Whether to eliminated isolated nodes
#' @param include.xyz Whether to include 3d location as vertex attribute
#' @return \code{igraph} object containing only nodes of neuron keeping original
#'   labels (\code{x$d$PointNo} => \code{V(g)$label}) and vertex indices 
#'   (\code{1:nrow(x$d)} => \code{V(g)$vid)}.
segmentgraph<-function(x, exclude.isolated=FALSE, include.xyz=FALSE){
  g=graph.empty()
  pointnos=x$d$PointNo
  sts<-if(is.null(x$SubTrees)) x$SegList else unlist(x$SubTrees,recursive=FALSE)
  topntail<-function(x) if(length(x)==1) x else x[c(1,length(x))]
  # just get head and tail of each segment
  simple_sts=lapply(sts,topntail)
  all_nodes=sort(unique(unlist(simple_sts)))
  # make empty graph with approriate nodes
  g=graph.empty(n=length(all_nodes))
  # store external pointnos
  igraph::V(g)$label=pointnos[all_nodes]
  # store original vertex ids
  igraph::V(g)$vid=all_nodes
  
  # handle the edges - first make list full edgelist
  el=EdgeListFromSegList(simple_sts)
  # convert from original vertex ids to vids of reduced graph
  elred=match(t(el),all_nodes)
  g=add.edges(g,elred)
  
  if(include.xyz){
    igraph::V(g)$x=x$d$X[all_nodes]
    igraph::V(g)$y=x$d$Y[all_nodes]
    igraph::V(g)$z=x$d$Z[all_nodes]
  }
  if(exclude.isolated){
    # remove any points with no neighbours
    isolated_vertices=igraph::V(g)[igraph::degree(g)==0]
    g=igraph::delete.vertices(graph=g,isolated_vertices)
  }
  g
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
  coreneuron=CoreNeuronFromGraph(gam, origin=origin)
  ANeuron[names(coreneuron)]<-coreneuron
  
  # also use depth first search to reorder SWC data starting from origin
  # note that this order will be well-defined when the tree is fully connected
  # but when there are subgraphs that are not connected to the origin the
  # numbering will be valid but have no meaning with respect to the new origin.
  dfs=graph.dfs(gam,root=origin,father=TRUE,neimode='all')
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
#' @param prune Whether to prune the graph of any vertices not connected by
#'  edges in the input graph (default TRUE)
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
    neurongraph.swc(x$d, directed=directed)
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
  .Deprecated('as.neuron.data.frame','nat')
  as.neuron(swc)
}

#' Contruct a graph to encode a neuron's connectivity
#' 
#' @details We make the following assumptions about neurons coming in
#' \itemize{ 
#'   \item They have an integer vertex label that need not start from 1 and that
#'   may have gaps 
#'   \item The edge list which defines connectivity specifies edges using pairs
#'   of vertex labels, _not_ raw vertex ids.
#' }
#' @details We make no attempt to determine the root points at this stage.
#' @details The raw vertex ids will be in the order of vertexlabels and can
#' therefore be used to index a block of vertex coordinates.
#' @param el A two columm matrix (start, end) defining edges
#' @param vertexlabels Integer labels for graph - the edge list is specified
#'   using these labels.
#' @return an \code{igraph} object with a vertex for each entry in vertexlabels,
#'   each vertex having a \code{label} attribute. All vertices are included
#'   whether connected or not.
neurongraph<-function(el, vertexlabels, directed=TRUE){
  if(any(duplicated(vertexlabels))) stop("Vertex labels must be unique!")
  # now translate edges into raw vertex_ids
  rawel=match(t(el), vertexlabels)
  g=graph(rawel, n=length(vertexlabels), directed=directed)
  V(g)$label=vertexlabels
  g
}

#' Construct neurongraph from a block of SWC format data
#' 
#' @param x Dataframe containingb block of SWC data
#' @param directed Logical indicating whether to make a directed graph
neurongraph.swc<-function(x, directed=TRUE){
  el=data.matrix(EdgeListFromSWC(x))
  neurongraph(el, x$PointNo, directed=directed)
}

#' Make SegList (and other core fields) from full graph of all nodes and origin
#' 
#' @details Uses a depth first search on the tree to reorder using the given 
#'   origin.
#' @details When the graph contains multiple subgraphs, only one will be chosen 
#'   as the master tree and used to construct the SegList of the resultant 
#'   neuron. However all subgraphs will be listed in the SubTrees element of the
#'   neuron and nTrees will be set appropriately.
#' @details When the graph vertices have a label attribute derived from PointNo,
#'   the origin is assumed to be specified with respect to the vertex labels
#'   rather than the raw vertex ids.
#' @param g An igraph
#' @param origin Root vertex, matched against labels (aka PointNo) when 
#'   available (see details)
#' @param dfs result of graph.dfs
#' @return A list with elements: 
#'   (NumPoints,StartPoint,BranchPoints,EndPoints,nTrees,NumSegs,SegList,
#'   [SubTrees])
#'   NB SubTrees will only be present when nTrees>1.
#' @export
#' @rdname CoreNeuron
#' @seealso \code{\link{graph.dfs},\link{RerootNeuron},\link{}}
CoreNeuronFromGraph<-function(g, origin=NULL, Verbose=TRUE){
  .Deprecated('as.neuron.ngraph','nat')
  if(!inherits(g,'ngraph')) class(g)=c("ngraph",class(g))
  as.neuron(g,origin=origin,Verbose=Verbose)
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

#' Construct neuron (inc SWC block) from Amira skeletonize format data
#' 
#' @param x list containing amira data with elements PointList,EdgeList,Origin
#' @return A neuron object containing both SegList and associated fields and the
#'   SWC format data block with coordinate data.
CoreNeuronFromAmiraSkel<-function(x, Verbose=FALSE){
  el=data.matrix(x$EdgeList)
  # make a doubly linked graph from this double edge list
  doubleg=neurongraph(el, x$PointList$PointNo, directed=TRUE)
  # make it undirected
  # TODO see if we can make appropriate directed graph rather than converting
  # to undirected.
  ug=as.undirected(doubleg,mode='collapse')
  # make seglist associated fields from that
  cn=CoreNeuronFromGraph(ug, origin=x$Origin, Verbose=Verbose)
  
  # now construct matching swc data block
  cn$d=cbind(PointNo=x$PointList$PointNo,Label=2,x$PointList[,c("X",'Y','Z','W')])
  # And recalculate parents
  RecalculateSWCData(as.neuron(cn))
}
