# runitImagesToDots.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitImagesToDots.R"))
require(RUnit)
require(R.matlab)

test.DotProperties<-function(){
	dots4=t(readMat(file.path(TestDir,"Geometry","SAKW13-1_dots4.mat"))[[1]])
	props4=readMat(file.path(TestDir,"Geometry","SAKW13-1_dots4_properties.mat"))
	# reorder matrices into R forms
	props4$alpha=as.numeric(props4$alpha)
	props4$vect=t(props4$vect)
	
	props4.new<-DotProperties(dots4)
	
	checkEqualsNumeric(props4.new$alpha,props4$alpha,tol=1e-6)
	# should probably just check that dot product is +/- 1
	checkEqualsNumeric(abs(props4.new$vect),abs(props4$vect),tol=1e-6)
	checkEqualsNumeric(props4.new$vect,props4$vect,tol=1e-6)
}