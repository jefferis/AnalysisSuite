# runitReadIGSFiles.R

# Some test neurons to verify that all is well for the
# IGS IO routines

require(RUnit)
test.ReadWriteIGSLandmarksSingle<-function(){
	testData=matrix(rnorm(15),ncol=3)
	rownames(testData)=letters[1:5]
	tmpfile=tempfile()
	WriteIGSLandmarks(testData,tmpfile)
	testData.new=ReadIGSLandmarks(tmpfile)
	unlink(tmpfile)
	checkEquals(testData,testData.new,tol=1e-6)	
}