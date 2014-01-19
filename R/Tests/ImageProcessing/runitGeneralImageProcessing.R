# runTestFile(file.path(TestDir,'ImageProcessing','runitGeneralImageProcessing.R'))

require(RUnit)

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

test.CMTKStatistics<-function(){
	if(is.null(cmtk.bindir())){
		message("Can't find CMTK statistics executable, skipping test")
		return()
	}
	lhmaskfile=file.path(TestDir,"Data","LHMask.nrrd")
	statsnrrd=file.path(TestDir,"Data","dataforstats.nrrd")
	a=CMTKStatistics(lhmaskfile)
	baseline_a=structure(list(min = 0, max = 1, mean = 0.22935, sdev = 0.42042, 
	    n = 125000L, Entropy = 0.53849, sum = 28669), .Names = c("min", 
	"max", "mean", "sdev", "n", "Entropy", "sum"), class = "data.frame", row.names = c(NA, 
	-1L))
	b=CMTKStatistics(statsnrrd)
	baseline_b=structure(list(min = 0, max = 100, mean = 8e-04, sdev = 0.28284, 
	    n = 125000L, Entropy = 1e-04, sum = 100), .Names = c("min", 
	"max", "mean", "sdev", "n", "Entropy", "sum"), class = "data.frame", row.names = c(NA, 
	-1L))
	checkEqualsNumeric(a,baseline_a)
	checkEqualsNumeric(b,baseline_b)
	c=CMTKStatistics(statsnrrd,mask=lhmaskfile)
	# my hacked version of statistics provides nnz
	if('nnz'%in%names(c)){
		baseline_c=structure(list(X.M = 0:1, min = c(0, 0), max = c(0, 100), mean = c(0, 
		0.00349), sdev = c(0, 0.5906), n = c(96331L, 28669L), nnz = 0:1, 
		    Entropy = c(0, 0.00039), sum = c(0, 100)), .Names = c("X.M", 
		"min", "max", "mean", "sdev", "n", "nnz", "Entropy", "sum"), class = "data.frame", row.names = c(NA, 
		-2L))
	} else {
		baseline_c = structure(list(MaskLevel = 0:1, min = c(0, 0), max = c(0, 100), mean = c(0, 
		0.00349), sdev = c(0, 0.5906), n = c(96331L, 28669L), Entropy = c(0, 
		0.00039), sum = c(0, 100)), .Names = c("MaskLevel", "min", "max", "mean", 
		"sdev", "n", "Entropy", "sum"), class = "data.frame", row.names = c(NA, 
		-2L))
	}
	checkEqualsNumeric(c,baseline_c)
}
