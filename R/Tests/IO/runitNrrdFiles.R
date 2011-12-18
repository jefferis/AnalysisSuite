# Tests of reading/writing of nrrd files and some other basic operations
# More sophisticated operations would belong in ImageProcessing section

test.AddOrReplaceNrrdHeaderField<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))
	lhmaskfile=file.path(ObjDir,"LHMask.nrrd")
	
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
	
	origlhmaskfile=file.path(ObjDir,"LHMask.nrrd")
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
