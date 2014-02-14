# runitNeuronFunctions.R

# Some test neurons to verify that all is well for the 
# AmiraFileFunctions rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(TestDir,"Data","neurons","runitNeuronFunctions.R"))
require(RUnit)
test.read.neuron<-function(){
	result.new=read.neuron(file.path(TestDir,"Data","neurons","Neurites.am"))
	result=read.neuron(file.path(TestDir,"Data","neurons","Neurites.am"))
	checkTrue(is.neuron(result.new))
	checkEquals(result,result.new,tol=1e-4,fieldsToExclude='d')
	# We don't both returning NeighbourCount as the 8th column
	checkEquals(result$d[,1:7],result.new$d,tol=1e-4)
	
	result.new=read.neuron(file.path(TestDir,"Data","neurons","SequentiallyBranchingTrace.traces"))
	result=ReadNeuronsFromLongairTraces(
		file.path(TestDir,"Data","neurons","SequentiallyBranchingTrace.traces"),Merge=TRUE)
	checkTrue(is.neuron(result.new))
	checkEquals(result.new,result)
	
	result.new=read.neuron(file.path(TestDir,"Data","neurons","LF28R.tasc"))
	result=ReadNeuronFromAsc(file.path(TestDir,"Data","neurons","LF28R.tasc"))
	checkTrue(is.neuron(result.new))
	checkEquals(result,result.new)
	
	checkException(read.neuron(file.path(TestDir,"Data","neurons","runitNeuronFunctions.R")),
		silent=TRUE)
}
# debug(test.read.neuron)

test.all.equal.neuron<-function(){
	a=read.neuron(file.path(TestDir,"Data","neurons","Neurites.am"))
	b=ReadNeuronsFromLongairTraces(
		file.path(TestDir,"Data","neurons","SequentiallyBranchingTrace.traces"))
	checkTrue(all.equal(a, a))
	checkTrue(!isTRUE(all.equal(a, b)))
	checkTrue(!isTRUE(all.equal(a, NULL)))
	checkTrue(!isTRUE(all.equal(a, 1)))
	checkTrue(!isTRUE(all.equal(a, NA)))
}
