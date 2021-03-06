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

test.AffineToCMTKRegistration<-function(){
  affmat=matrix(c(1.09698704245255, -0.0574906547972094, -0.0575695518672382, 
                  0, 0.0794174055868413, 0.895837910136773, 0.0454677297831557, 
                  0, 0.151929624282028, -0.0103700734989691, 0.994640563641534, 
                  0, 87.462778082031, 56.3015147160819, 8.57028084743792, 1),
                ncol=4)
  reg=AffineToCMTKRegistration(affmat)
  checkEquals(attr(reg,"version"),numeric_version("2.4"))
}

test.CMTKRegFromAmira<-function(){
  # test by 1. writing an affine to an amira registration and then converting to
  #   CMTK registration
  # 2. writing the same affine direct to a CMTK registration
  # 3. comparing the two
  affmat=matrix(c(1.09698704245255, -0.0574906547972094, -0.0575695518672382, 
                  0, 0.0794174055868413, 0.895837910136773, 0.0454677297831557, 
                  0, 0.151929624282028, -0.0103700734989691, 0.994640563641534, 
                  0, 87.462778082031, 56.3015147160819, 8.57028084743792, 1),
                ncol=4)
  cmtkdirect=tempfile(fileext='_direct.list')
  on.exit(unlink(cmtkdirect,recursive=TRUE))
  cmtk.mat2dof(affmat,f=cmtkdirect)
  
  amirareg=tempfile(fileext='.hxtransform')
  on.exit(unlink(amirareg),add=TRUE)
  # nb this will write in the same form as our Amira scripts expect,
  # i.e. column-wise
  cat(affmat,file=amirareg)
  
  cmtkindirect=tempfile(fileext='_indirect.list')
  on.exit(unlink(cmtkindirect,recursive=TRUE),add=TRUE)
  reg=CMTKRegFromAmira(amirareg=amirareg,cmtkregfolder=cmtkindirect)
  
  cmtkdirect=readLines(cmtkdirect)
  cmtkindirect=readLines(cmtkindirect)
  checkEquals(cmtkdirect,cmtkindirect)
}
