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
}

test.rootpoints<-function(){
  checkEquals(rootpoints(testn),1)
  testd=data.frame(PointNo=1:6,Label=2,
                   X=c(1:5,3),Y=c(rep(1,5),2),Z=0,W=NA,
                   Parent=c(-1,1:4,-1))
  testn.floating=SWC2Neuron(testd,'test')
  checkEquals(rootpoints(testn.floating),1)
  testd2=rbind(testd,c(7,2,7,7,0,NA,6))
  testn.2trees=SWC2Neuron(testd2,'test')
  # check that we get two roots when there are indeed 2 roots
  rps=rootpoints(as.igraph(testn.2trees))
  checkEquals(rps,c(1,6))
}
