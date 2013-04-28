# runitInterpolateAlongNeuron.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitInterpolateAlongNeuron.R"))
require(RUnit)


test.Interpolate<-function(){
	n=read.neuron(file.path(TestDir,"Data","neurons","Neurites.am"))
	n1=InterpolateAlongNeuron(n,stepSize=1)

	u=read.neuron(file.path(TestDir,"Data","neurons","UnbranchedNeurite.am"))
	u1=InterpolateAlongNeuron(u,stepSize=1)
	plot3d(u1)
}
