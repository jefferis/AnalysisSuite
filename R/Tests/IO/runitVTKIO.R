# Tests for reading & writing VTK data

# runTestFile(file.path(TestDir,"IO","runitVTKIO.R"))

require(RUnit)

test.ReadWriteVTKLandmarks<-function(){
	testData=matrix(rnorm(15),ncol=3)
	tmpfile=tempfile()

	WriteVTKLandmarks(tmpfile,testData)
	testData.new=ReadVTKLandmarks(tmpfile)
	unlink(tmpfile)
	checkEqualsNumeric(testData,testData.new,tol=1e-6)	
}
