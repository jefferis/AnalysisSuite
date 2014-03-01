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

test.getBoundingBox<-function(){
	lhfile=file.path(TestDir,'Data','LHMask.nrrd')
	lh=Read3DDensityFromNrrd(lhfile)
	checkEqualsNumeric(c(0, 68.6, 0, 68.6, 0, 68.6),getBoundingBox(lh),tol=1e-6)
	checkEquals(getBoundingBox(lh),getBoundingBox(lhfile))
	checkEqualsNumeric(structure(c(-0.7, 69.3, -0.7, 69.3, -0.7, 69.3), .Dim = 2:3),
		getBounds(lh),tol=1e-6)
	checkEquals(getBounds(lh),getBounds(lhfile))
}

test.makeScaleBar<-function(){
	lhfile=file.path(TestDir,'Data','LHMask.nrrd')
	lh=Read3DDensityFromNrrd(lhfile,ReadByteAsRaw='none')
	p=projection(lh)
	rval=image.gjdens(p)
	rval_baseline=structure(list(zlim = c(0, 54.6), nlevels.actual = 29L, nlevels.orig = 20, 
	    levels = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 
	    26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 
	    56), colors = c("#000080", "#001C8E", "#00389C", "#0055AA", 
	    "#0071B8", "#008DC6", "#00AAD4", "#00C6E2", "#00E2F0", "#00FFFF", 
	    "#1CFFE2", "#38FFC6", "#55FFAA", "#71FF8D", "#8DFF71", "#AAFF54", 
	    "#C6FF38", "#E2FF1C", "#FFFF00", "#FFE200", "#FFC600", "#FFAA00", 
	    "#FF8D00", "#FF7100", "#FF5500", "#FF3800", "#FF1C00", "#FF0000"
	    )), .Names = c("zlim", "nlevels.actual", "nlevels.orig", 
	"levels", "colors"))
	checkEquals(rval,rval_baseline)
	makeScaleBar(rval)
}
