# runitDensityFunctions.R
# functions to test affine functions

# runTestFile(file.path(CodeDir,"runitDensityFunctions.R"))
require(RUnit)

test.expand.grid.gjdens<-function(){
	# make a dummy gjdens object and check
	# that my expand.grid and the default one give the
	# same answer
	dims=c(10,5,2)
	d=array(dim=dims)
	b=c(0,168.7,0,168,0,88)
	attr(d,"BoundingBox")=b
	attr(d,"x")=seq(b[1],b[2],len=dims[1])
	attr(d,"y")=seq(b[3],b[4],len=dims[2])
	attr(d,"z")=seq(b[5],b[6],len=dims[3])
	checkEqualsNumeric(expand.grid.gjdens(d),
		data.matrix(expand.grid( attr(d,"x"),attr(d,"y"),attr(d,"z") ))
	)
}

