# runitDensityFunctions.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Core","runitDensityFunctions.R"))
require(RUnit)

makeDummyDensityData<-function()
{
	dims=c(10,5,2)
	d=array(dim=dims)
	b=c(0,168.7,0,168,0,88)
	attr(d,"BoundingBox")=b
	attr(d,"x")=seq(b[1],b[2],len=dims[1])
	attr(d,"y")=seq(b[3],b[4],len=dims[2])
	attr(d,"z")=seq(b[5],b[6],len=dims[3])
	d
}

test.expand.grid.gjdens<-function(){
	# make a dummy gjdens object and check
	# that my expand.grid and the default one give the
	# same answer
	d=makeDummyDensityData()
	checkEqualsNumeric(expand.grid.gjdens(d),
		data.matrix(expand.grid( attr(d,"x"),attr(d,"y"),attr(d,"z") ))
	)
}

test.ijkpos.gjdens<-function(x)
{
	d=makeDummyDensityData()
	pixels=ijkpos.gjdens(d,rbind(c(95.7,60.7,0.7),c(95.7,60.7,0.7)),round=T)
	rightans=structure(c(6, 6, 2, 2, 1, 1), .Dim = 2:3)
	checkEquals(pixels,rightans)
}

test.xyzpos.gjdens<-function(x)
{
	d=makeDummyDensityData()
	corners=matrix(getBoundingBox(d),ncol=3)
	pixels=rbind(c(1,1,1),dim(d))
	checkEquals(xyzpos.gjdens(d,pixels),corners)
}
