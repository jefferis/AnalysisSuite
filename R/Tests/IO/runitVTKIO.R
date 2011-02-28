# Tests for reading & writing VTK data

# runTestFile(file.path(TestDir,"IO","runitVTKIO.R"))

require(RUnit)

test.ReadWriteVTKLandmarks<-function(){
	testData=matrix(rnorm(15),ncol=3)

	tmpfile=tempfile()
	on.exit(unlink(tmpfile))

	WriteVTKLandmarks(tmpfile,testData)
	testData.new=ReadVTKLandmarks(tmpfile)

	checkEqualsNumeric(testData,testData.new,tol=1e-6)	
}
