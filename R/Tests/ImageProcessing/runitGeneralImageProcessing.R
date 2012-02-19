test.ReadWriteIdentityRegistration<-function(){
	ireg=structure(c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0), 
		.Dim = c(5L, 3L), .Dimnames = list(c("xlate", "rotate", "scale", "shear", "center"),
		 c("X", "Y", "Z")))
	f=WriteIdentityRegistration()
	on.exit(unlink(f,recursive=TRUE))
	x=ReadIGSRegistration(f)
	checkEqualsNumeric(ireg,do.call(rbind,x$affine_xform))
}

test.AutoCropNrrd<-function(){
	tmpdir=tempfile()
	dir.create(tmpdir)
	lhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	AutoCropNrrd(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"))
	h=ReadNrrdHeader(file.path(tmpdir,"LHMask.nrrd"))
	imageSize=c(44,42,40)
	spaceOrigin=c(4.2,7.0,5.6)
	checkEqualsNumeric(h$sizes,imageSize,msg="Mismatch with expected image dimensions in pixels")
	checkEqualsNumeric(h$`space origin`,spaceOrigin,msg="Mismatch with expected image origin (physical coords)",tol=1e-6)
	unlink(tmpdir,recursive=TRUE)	
}
