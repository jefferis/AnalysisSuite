# runitGraphTheory.R
# functions to test GraphTheory functions

# runTestFile(file.path(TestDir,"Geometry","runitGraphTheory.R"))
require(RUnit)

# make a very simple neuron

# 1 2 3b 4 5
#     6
testd=data.frame(PointNo=1:6,Label=2,
  X=c(1:5,3),Y=c(rep(1,5),2),Z=0,W=NA,
  Parent=c(-1,1:4,3))
testn=SWC2Neuron(testd,'test')

# isomorphic
# 1 2 3b 5 6
#     4
testn2=SWC2Neuron(data.frame(
  PointNo=1:6,Label=2,
  X=c(1:3,3,4,5),Y=c(rep(1,3),2,1,1),Z=0,W=NA,
  Parent=c(-1,1:3,3,5))
  ,'test')

# different by vertex position (isomorphic to 1 by connectivity)
# 1 2 3b 4
#     5
#     6
testn3=SWC2Neuron(data.frame(
  PointNo=1:6,Label=2,
  X=c(1:4,3,3),Y=c(rep(1,4),2,3),Z=0,W=NA,
  Parent=c(-1,1:3,3,5))
  ,'test')

test.as.igraph<-function(){
  g1=as.igraph(testn)
  g2=as.igraph(testn2)
  g3=as.igraph(testn3)
  checkEquals(g2,g3)
  
  g1s=as.igraph(testn,method='seglist')
  g2s=as.igraph(testn2,method='seglist')
  g3s=as.igraph(testn3,method='seglist')
  checkEquals(g2s,g3s)
  
  checkTrue(graph.isomorphic(g1,g2))
  checkTrue(graph.isomorphic(g1s,g1))
  checkTrue(graph.isomorphic(g2s,g2))
  
  # check equivalence of 2 methods for a more complicated graph
  am3da=read.neuron(file.path(TestDir,"Data","neurons","testneuron_am3d_ascii.am"))
  g2<-as.igraph(am3da,meth='swc')
  g2s<-as.igraph(am3da,meth='seglist')
  checkTrue(graph.isomorphic(g2s,g2))
}

testSWC2Neuron<-function(){
  checkEquals(testn$NumSegs,3)
  checkEquals(testn$NumPoints,6)
  checkEquals(testn$EndPoints,c(1,5,6))
  checkEquals(testn$BranchPoints,3)
  checkEquals(testn$StartPoint,1)
  checkTrue(is.neuron(testn))
}

testParseSWCTreeEquivalent<-function(){
	# ParseSWCTree has been replaced by SWC2Neuron
	# Check they produce the same result
	checkTrue(all.equal(ParseSWCTree(testd,'test'),testn))
}

testRerootNeuron<-function(){
  rn=RerootNeuron(testn)
  identicalFields=c('NumSegs','NumPoints','EndPoints','BranchPoints')
  checkTrue(is.neuron(rn))
  checkEquals(rn[identicalFields],testn[identicalFields])
  checkEquals(rn$StartPoint,1)
  checkEquals(rn$d$Parent,c(-1, 1, 2, 3, 4, 3))

  rn=RerootNeuron(testn,origin=2)
  checkEquals(rn$StartPoint,2)
  checkEquals(rn$d$Parent,c(2, -1, 2, 3, 4, 3))
  
  checkEquals(RerootNeuron(testn,origin=3)$d$Parent,c(2, 3, -1, 3, 4, 3))
  checkEquals(RerootNeuron(testn,origin=4)$d$Parent,c(2, 3, 4, -1, 4, 3))
  checkEquals(RerootNeuron(testn,origin=5)$d$Parent,c(2, 3, 4, 5, -1, 3))
  checkEquals(RerootNeuron(testn,origin=6)$d$Parent,c(2, 3, 6, 3, 4, -1))
}

test.CoreNeuronFromGraph<-function(){
  g<-as.igraph(testn)
  cn=CoreNeuronFromGraph(g,1)
  checkEquals(testn$SegList,cn$SegList)
  
  # vertex labels with gaps
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(2:4,6,7,9))
  sl=list(c(1,3,2,4),c(4,5),c(4,6))
  checkEquals(graph2seglist(g,origin=1),sl)
  checkEquals(CoreNeuronFromGraph(g,origin=2)$SegList,sl)
  # same but no origin specified (should give same result)
  checkEquals(CoreNeuronFromGraph(g)$SegList,sl)
  
  # same but different origin
  sl2=list(c(3,1),c(3,2,4),c(4,5),c(4,6))
  checkEquals(CoreNeuronFromGraph(g,origin=4)$SegList,sl2)
  
  # same connectivity but one extra (floating) point at end
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(2:4,6,7,9,10))
  n=CoreNeuronFromGraph(g,origin=4)
  checkEquals(n$SegList,sl2)
  checkEquals(n$nTrees,2)
  checkEquals(n$SubTrees,list(sl2,list(7)))
  
  # same connectivity but with extra (floating) points at start and end
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(1:4,6,7,9,10))
  # this will shift all vertex ids by 1
  sl3=lapply(sl2,'+',1)
  n=CoreNeuronFromGraph(g,origin=4)
  checkEquals(n$SegList,sl3)
  checkEquals(n$nTrees,3)
  checkEquals(n$SubTrees,list(sl3,list(1),list(8)))
  
  # 3 separate subgraphs of length 3,4,5
  g=neurongraph(c(0,1,1,2, 3,4,4,5,5,6, 7,8,8,9,9,10,10,11),vertexlabels=0:11)
  n=CoreNeuronFromGraph(g,origin=0)
  checkEquals(n$SegList,list(c(1,2,3)))
  n2=CoreNeuronFromGraph(g,origin=3)
  checkEquals(n2$SegList,list(c(4,5,6,7)))
  n3=CoreNeuronFromGraph(g,origin=7)
  checkEquals(n3$SegList,list(c(8,9,10,11,12)))
  # check that it picks largest subgraph when no origin specified
  n4=CoreNeuronFromGraph(g)
  checkEquals(n4,n3)
}

test.nodes<-function(){
  checkEquals(rootpoints(testn),1)
  checkEquals(endpoints(testn),c(1,5,6))
  checkEquals(branchpoints(testn),3)
  checkException(rootpoints(testn,subtrees=2),silent=T)
  # now check that we can cope if nTrees is not set
  # (as is sometimes true for old neurons)
  testn$nTrees=NULL
  checkEquals(rootpoints(testn),1)
  checkException(rootpoints(testn,subtrees=2),silent=T)
  
  # now more complicated neurons - isloated point 
  testd=data.frame(PointNo=1:6,Label=2,
                   X=c(1:5,3),Y=c(rep(1,5),2),Z=0,W=NA,
                   Parent=c(-1,1:4,-1))
  testn.floating=SWC2Neuron(testd,'test')
  checkEquals(rootpoints(testn.floating),1)
  checkEquals(rootpoints(testn.floating, subtrees=1:2, exclude.isolated=FALSE),c(1,6))
  
  testd2=rbind(testd,c(7,2,7,7,0,NA,6))
  testn.2trees=SWC2Neuron(testd2,'test')
  # check that we get two roots when there are indeed 2 roots
  rps=rootpoints(as.igraph(testn.2trees))
  checkEquals(rps,c(1,6))
}

test.graph2seglist<-function(){
  # simple linear graph
  g=graph(c(1, 2, 2, 3))
  sl=list(c(1, 2, 3))
  checkEquals(graph2seglist(g), sl)
  
  # simple linear graph with different vids
  g=graph(c(1, 2, 2, 3))
  igraph::V(g)$vid=3:5
  sl=list(3:5)
  checkEquals(graph2seglist(g), sl)
  
  # simple linear graph with different vids and different origin
  g=graph(c(1, 2, 2, 3))
  igraph::V(g)$vid=3:5
  sl=list(5:3)
  checkEquals(graph2seglist(g, origin=5), sl)
  
  # simple linear graph with different vids and origin at centre, resulting
  # in a branched seglist
  g=graph(c(1, 2, 2, 3))
  igraph::V(g)$vid=3:5
  sl=list(4:3,4:5)
  checkEquals(graph2seglist(g, origin=4), sl)
  
  # multiple subtrees -> exception since seglist only defined for 1 subtree
  g=graph(c(1,2,2,3,3,4,5,6))
  checkException(graph2seglist(g),silent=TRUE)
  
  # cyclic graph -> exception since seglist is undefined
  g=graph(c(1, 2, 2, 3, 3, 1))
  checkException(graph2seglist(g),silent=TRUE)
  
  # single floating point
  g=graph(NULL,n=1)
  checkEquals(graph2seglist(g),list(1))
  
  # single floating point with different vid
  igraph::V(g)$vid=4
  checkEquals(graph2seglist(g),list(4))
  
  # trifurcation
  g=graph(c(1,2, 2,3, 2,4, 2,5, 5,6, 6,7))
  sl=list(c(1,2),c(2,3),c(2,4),c(2,5,6,7))
  checkEquals(graph2seglist(g),sl)
  # undirected equivalent - nb origin must be specified
  checkEquals(graph2seglist(as.undirected(g),origin=1),sl)
  
  # rapid branching
  g=graph(c(1,2, 2,3, 2,4, 4,5, 4,6))
  sl=list( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6) )
  checkEquals(graph2seglist(g),sl)
  
  # different root
  g=graph(c(1,2, 2,3, 2,4, 4,5, 4,6))
  sl=list( c(6,4),c(4,2),c(2,1),c(2,3),c(4,5) )
  checkEquals(graph2seglist(g,origin=6),sl)
  
  # non-sequential numbering
  g=graph(c(1,2, 2,6, 2,4, 4,5, 4,3))
  sl<-list( c(1,2),c(2,4),c(4,3),c(4,5),c(2,6) )
  checkEquals(graph2seglist(g,origin=1),sl)
  
  # non-sequential numbering with vertex labels
  # in this case we imagine that there are a set of vertices with PointNo
  # 2,3,4,6,7,9
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(2:4,6,7,9))
  sl=list(c(1,3,2,4),c(4,5),c(4,6))
  checkEquals(graph2seglist(g,origin=1),sl)
  # same but with a different origin
  # NB origin is defined in terms of sequential raw vertex id so 
  # origin=3 means that the origin is the vertex with label=4
  sl2=list(c(3,1),c(3,2,4),c(4,5),c(4,6))
  checkEquals(graph2seglist(g,origin=3),sl2)
}
