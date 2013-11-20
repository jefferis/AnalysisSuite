# Tests of reading/writing of nrrd files and some other basic operations
# More sophisticated operations would belong in ImageProcessing section

test.AddOrReplaceNrrdHeaderField<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))
	lhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	
	# check we error out when providing a bad field
	checkException(
		AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
			c(fakefield="(1,2,3)")),silent=TRUE)
	
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c(`space origin`="(2,2,2)"))
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	newSpaceOrigin=c(2,2,2)
	checkEqualsNumeric(h$`space origin`,newSpaceOrigin,msg="Mismatch with replaced image origin (physical coords)",tol=1e-6)
	
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c(`sample units`="ppm"),Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	sampleunits="ppm"
	checkEquals(h$`sampleunits`,sampleunits,msg="Mismatch with added field (sample units)")

	# change both at same time
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c(`space origin`="(2,2,2)",`sample units`="ppm"),Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	checkEquals(h$`sampleunits`,sampleunits,msg="Mismatch with added field (sample units)")
	checkEqualsNumeric(h$`space origin`,newSpaceOrigin,msg="Mismatch with replaced image origin (physical coords)",tol=1e-6)
	
	# add a comment
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c("# My interesting comment"),Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	checkTrue(any(grepl("My interesting comment",attr(h,'headertext'))))

	# add two comments
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c("# My interesting comment","#Another interesting comment"),Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	checkTrue(sum(grepl("interesting comment",attr(h,'headertext'))) == 2 )
	
	# ... and a field at the same time
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		c("# My interesting comment",`space origin`="(2,2,2)"),Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	checkTrue(any(grepl("My interesting comment",attr(h,'headertext'))))
	checkEqualsNumeric(h$`space origin`,newSpaceOrigin,msg="Mismatch with replaced image origin (physical coords)",tol=1e-6)
}

test.AddOrReplaceNrrdHeaderFieldInPlace<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))
	
	origlhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	lhmaskfile=file.path(tmpdir,basename(origlhmaskfile))
	file.copy(origlhmaskfile,lhmaskfile)
	
	# Check that we error out when trying to overwrite without Force = TRUE
	checkException(
		AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=lhmaskfile,c("space origin"="(0,0,0)")),
		silent=TRUE)
	# Check replacement 
	checkException(
		      AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
		          c("fakefield"="(1,2,3)")), silent=TRUE)
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=lhmaskfile,c("space origin"="(2,2,2)"),Force=TRUE)
	h=ReadNrrdHeader(lhmaskfile)
	horig=ReadNrrdHeader(origlhmaskfile)
	newSpaceOrigin=c(2,2,2)
	checkEqualsNumeric(h$`space origin`,newSpaceOrigin,msg="Mismatch with expected image origin (physical coords)",tol=1e-6)
}

test.AddOrReplaceNrrdHeaderFieldDetached<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))
	
	origlhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	lhmaskfile=file.path(tmpdir,basename(origlhmaskfile))
	file.copy(origlhmaskfile,lhmaskfile)
	
	# Write a detached header file
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=paste(lhmaskfile,'-det.nhdr',sep=""),
      c(content=NrrdCrc(lhmaskfile)),Detached = TRUE)
  
  # minmax provides a simple way to check that file can still be read
  checkEqualsNumeric(NrrdMinMax(paste(lhmaskfile,'-det.nhdr',sep="")),c(0,1),
      "Check that nhdr can be read and produces correct minmax")
}

test.NrrdMakeDetachedHeaderForNrrd<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))
	
	origlhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	lhmaskfile=file.path(tmpdir,basename(origlhmaskfile))
	file.copy(origlhmaskfile,lhmaskfile)
	
	nhdrvec=NrrdMakeDetachedHeaderForNrrd(lhmaskfile,NA)
  checkTrue(!file.exists(file.path(tmpdir,"LHMask.nrrd.nhdr")),
      "Don't write header if we nhdr=NA")
  checkTrue(is.character(nhdrvec),"Check that we get a character vector back")
	NrrdMakeDetachedHeaderForNrrd(lhmaskfile)
  NrrdMakeDetachedHeaderForNrrd(origlhmaskfile,file.path(tmpdir,'LHMask-orig.nhdr'))
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd.nhdr"))
  h2=ReadNrrdHeader(file.path(tmpdir,"LHMask-orig.nhdr"))
  # minmax provides a simple way to check that file can still be read
  checkEqualsNumeric(NrrdMinMax(file.path(tmpdir,"LHMask.nrrd.nhdr")),c(0,1),
      "Check that nhdr can be read and produces correct minmax")
  checkEquals(h$datafile,basename(lhmaskfile),
      "Check relative path for detached header right next to data file")
  checkEquals(h2$datafile,origlhmaskfile,
      "Check absolute path for detached header not next to data file")
  
}

test.is.nrrd<-function(){
  tmpdir=tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir,recursive=TRUE))
  
  origlhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
  lhmaskfile=file.path(tmpdir,sub("\\.nrrd$",'.something',basename(origlhmaskfile)))
  file.copy(origlhmaskfile,lhmaskfile)

  checkTrue(is.nrrd(origlhmaskfile))
  checkTrue(is.nrrd(origlhmaskfile,TrustSuffix=TRUE))
  checkException(is.nrrd(origlhmaskfile,ReturnVersion=TRUE,TrustSuffix=TRUE),
                msg="Check error when asking for nrrd version using suffix",
                silent=TRUE)
  checkEqualsNumeric(is.nrrd(origlhmaskfile,ReturnVersion=TRUE),4)
  
  checkTrue(is.nrrd(lhmaskfile),msg='Can identify nrrd with funny suffix')
  checkEquals(is.nrrd(lhmaskfile,TrustSuffix=TRUE),FALSE)
}
