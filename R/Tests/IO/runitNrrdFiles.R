# Tests of reading/writing of nrrd files and some other basic operations
# More sophisticated operations would belong in ImageProcessing section

test.AddOrReplaceNrrdHeaderField<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	lhmaskfile=file.path(ObjDir,"LHMask.nrrd")
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),"space origin","(0,0,0)")
	AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),"fakefield","(1,2,3)",Force=TRUE)
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	newSpaceOrigin=c(0,0,0)
	fakefield=c(1,2,3)
	checkEqualsNumeric(h$`space origin`,newSpaceOrigin,msg="Mismatch with expected image origin (physical coords)",tol=1e-6)
	checkEqualsNumeric(h$`fakefield`,fakefield,msg="Mismatch with fake field",tol=1e-6)
	unlink(tmpdir,recursive=TRUE)
}
