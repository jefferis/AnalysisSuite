test.NrrdFlip<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	on.exit(unlink(tmpdir,recursive=TRUE))

	origlhmaskfile=file.path(ObjDir,"LHMask.nrrd")
	lhmaskfile=file.path(tmpdir,"LHMask.nrrd")
	# make a copy of the file, ensuring that it has content field
	AddOrReplaceNrrdHeaderField(origlhmaskfile,lhmaskfile,
		c(content=NrrdCrc(origlhmaskfile)),Force=T)
	
	lhm=Read3DDensityFromNrrd(lhmaskfile)
	h=ReadNrrdHeader(lhmaskfile)
	
	# flip back and forth along x axis
	NrrdFlip(lhmaskfile,axes=c(0,0))
	lhm00=Read3DDensityFromNrrd(file.path(tmpdir,"LHMask-flip00.nrrd"))
  print(attributes(lhm))
	print(attributes(lhm00))
	checkEquals(lhm,lhm00,"Check that flipping twice along an axis returns original")
	  checkEquals(NrrdCrc(file.path(tmpdir,"LHMask-flip00.nrrd")),NrrdCrc(lhmaskfile),
		"Check that flipping twice along an axis leaves Nrrd CRC unchanged")
	# flip x then z
	NrrdFlip(lhmaskfile,axes=c(0,2))
	# and z then x
	NrrdFlip(lhmaskfile,axes=c(2,0))
	checkEquals(NrrdCrc(file.path(tmpdir,"LHMask-flip02.nrrd")),
		NrrdCrc(file.path(tmpdir,"LHMask-flip20.nrrd")),
		"Check that flipping x then z and z then x give identical Nrrd CRC")
	# z then x separately
	NrrdFlip(lhmaskfile,axes=2)
	NrrdFlip(file.path(tmpdir,"LHMask-flip2.nrrd"),axes=0)
	checkEquals(NrrdCrc(file.path(tmpdir,"LHMask-flip2-flip0.nrrd")),
		NrrdCrc(file.path(tmpdir,"LHMask-flip20.nrrd")),
		"Check that flipping z then x in one step or two gives identical Nrrd CRC")
	# check content
	h20=ReadNrrdHeader(file.path(tmpdir,"LHMask-flip20.nrrd"))
	checkEquals(h20$content,
		ReadNrrdHeader(file.path(tmpdir,"LHMask-flip2-flip0.nrrd"))$content,
		"Check that content description is identical for one or two step flip")
	
	expected=paste("flip(flip(",h$content,",2),0)",sep="")
	checkEquals(h20$content,paste("flip(flip(",h$content,",2),0)",sep=""),
		"Check content field is as expected after two flips")
}
#debug(test.NrrdFlip)