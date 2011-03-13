# runitNeuronFunctions.R

# Some test neurons to verify that all is well for the 
# AmiraFileFunctions rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(TestDir,"IO","runitNeuronFunctions.R"))
require(RUnit)

test.read.neuron<-function(){
	fieldsToCheck=c("NeuronName", "NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
	"NumSegs", "SegList", "nTrees", "d", "OrientInfo")

	result.new=read.neuron(file.path(TestDir,"IO","Neurites.am"))
	result=ReadNeuronFromAM3D(file.path(TestDir,"IO","Neurites.am"))

	checkEquals(result[fieldsToCheck],result.new[fieldsToCheck],tol=1e-6)
	
	result.new=read.neuron(file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"))
	result=ReadNeuronsFromLongairTraces(file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"))
	
	checkException(read.neuron(file.path(TestDir,"IO","runitNeuronFunctions.R")))
}

test.all.equal.neuron<-function(){
	a=ReadNeuronFromAM3D(file.path(TestDir,"IO","Neurites.am"))
	b=ReadNeuronsFromLongairTraces(file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"))
	checkTrue(all.equal(a, a))
	checkTrue(!isTRUE(all.equal(a, b)))
	checkTrue(!isTRUE(all.equal(a, NULL)))
	checkTrue(!isTRUE(all.equal(a, 1)))
	checkTrue(!isTRUE(all.equal(a, NA)))
}
