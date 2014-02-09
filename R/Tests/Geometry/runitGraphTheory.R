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
  sl=seglist(c(1,3,2,4),c(4,5),c(4,6))
  checkEquals(as.seglist(g,origin=1),sl)
  checkEquals(CoreNeuronFromGraph(g,origin=2)$SegList,sl)
  # same but no origin specified (should give same result)
  checkEquals(CoreNeuronFromGraph(g)$SegList,sl)
  
  # same but different origin
  sl2=seglist(c(3,1),c(3,2,4),c(4,5),c(4,6))
  checkEquals(CoreNeuronFromGraph(g,origin=4)$SegList,sl2)
  
  # same connectivity but one extra (floating) point at end
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(2:4,6,7,9,10))
  n=CoreNeuronFromGraph(g,origin=4)
  checkEquals(n$SegList,sl2)
  checkEquals(n$nTrees,2)
  checkEquals(n$SubTrees,list(sl2,seglist(7)))
  
  # same connectivity but with extra (floating) points at start and end
  g=neurongraph(c(2,4,4,3,3,6,6,9,6,7),vertexlabels=c(1:4,6,7,9,10))
  # this will shift all vertex ids by 1
  sl3=as.seglist(lapply(sl2,'+',1))
  n=CoreNeuronFromGraph(g,origin=4)
  checkEquals(n$SegList,sl3)
  checkEquals(n$nTrees,3)
  checkEquals(n$SubTrees,list(sl3,seglist(1),seglist(8)))
  
  # 3 separate subgraphs of length 3,4,5
  g=neurongraph(c(0,1,1,2, 3,4,4,5,5,6, 7,8,8,9,9,10,10,11),vertexlabels=0:11)
  n=CoreNeuronFromGraph(g,origin=0)
  checkEquals(n$SegList,seglist(c(1,2,3)))
  n2=CoreNeuronFromGraph(g,origin=3)
  checkEquals(n2$SegList,seglist(c(4,5,6,7)))
  n3=CoreNeuronFromGraph(g,origin=7)
  checkEquals(n3$SegList,seglist(c(8,9,10,11,12)))
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
