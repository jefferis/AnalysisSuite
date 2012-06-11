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
  checkEquals(testn$StartPoint,1)
  checkEquals(testn$d$Parent,c(-1, 1, 2, 3, 4, 3))
  rn=RerootNeuron(testn,2)
  checkEquals(rn$StartPoint,2)
}