# runitNeuronFunctions.R

# Some test neurons to verify that all is well for the 
# AmiraFileFunctions rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(TestDir,"IO","runitNeuronFunctions.R"))
require(RUnit)
test.read.neuron<-function(){
	result.new=read.neuron(file.path(TestDir,"IO","Neurites.am"))
	result=ReadNeuronFromAM3D(file.path(TestDir,"IO","Neurites.am"))
	checkTrue(is.neuron(result.new))
	checkEquals(result,result.new)
	
	result.new=read.neuron(file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"))
	result=ReadNeuronsFromLongairTraces(
		file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"),Merge=TRUE)
	checkTrue(is.neuron(result.new))
	checkEquals(result.new,result)
	
	result.new=read.neuron(file.path(TestDir,"IO","LF28R.tasc"))
	result=ReadNeuronFromAsc(file.path(TestDir,"IO","LF28R.tasc"))
	checkTrue(is.neuron(result.new))
	checkEquals(result,result.new)
	
	checkException(read.neuron(file.path(TestDir,"IO","runitNeuronFunctions.R")),
		silent=TRUE)
}
# debug(test.read.neuron)

test.all.equal.neuron<-function(){
	a=ReadNeuronFromAM3D(file.path(TestDir,"IO","Neurites.am"))
	b=ReadNeuronsFromLongairTraces(
		file.path(TestDir,"IO","SequentiallyBranchingTrace.traces"))
	checkTrue(all.equal(a, a))
	checkTrue(!isTRUE(all.equal(a, b)))
	checkTrue(!isTRUE(all.equal(a, NULL)))
	checkTrue(!isTRUE(all.equal(a, 1)))
	checkTrue(!isTRUE(all.equal(a, NA)))
}
