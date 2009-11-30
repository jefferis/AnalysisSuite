# General Image Processing Functions combining CMTK / UNU or my R Native Code
# distinct from ImageAnalyisFunctions which are more statistical in nature
# this is really about munging files

ResampleAndFlipMasks<-function(masks,outdir,FlipBridgingReg,flipAxis=c("X","Y","Z"),targetspec,
	suffix="-resampled",gzip=TRUE){
	flipAxis=match.arg(flipAxis)
	if(!file.exists(outdir)) dir.create(outdir)
	resampledfiles=ResampleMasks(masks,outdir,FlipBridgingReg,flipAxis,targetspec,suffix=suffix,gzip)
	FlipAndORMasks(resampledfiles,outdir,FlipBridgingReg,flipAxis,gzip)
}

ResampleMasks<-function(masks,outdir,FlipBridgingReg,flipAxis=c("X","Y","Z"),targetspec,
	suffix="-resampled",gzip=TRUE){
	flipAxis=match.arg(flipAxis)
	if(!file.exists(outdir)) dir.create(outdir)
	resampledfiles=file.path(outdir,
		sub(".nrrd$",paste(suffix,".nrrd",sep=""),basename(masks)))
	for (i in seq(masks)){
		if(!exists("identityReg") || !file.exists(identityReg))
			identityReg=WriteIdentityRegistration()
		# resample - use reformatx to do this, to ensure that we get the same result
		ReformatImage(masks[i],target=targetspec, registrations=identityReg,
			filesToIgnoreModTimes=identityReg, OverWrite='update',
			output=resampledfiles[i],reformatoptions="-v --pad-out 0 --nn",dryrun=FALSE)
	}
	return(resampledfiles)
}

FlipAndORMasks<-function(masks,outdir,FlipBridgingReg,flipAxis=c("X","Y","Z"),gzip=TRUE){
	flipAxis=match.arg(flipAxis)
	if(!file.exists(outdir)) dir.create(outdir)
	for (infile in masks){
		flippedresampledfile=file.path(outdir,
			sub(".nrrd$","-flip.nrrd",basename(infile)))
		orfile=file.path(outdir,
			sub(".nrrd$","-OR.nrrd",basename(infile)))

		# make flipping registration
		if(!exists("horizontalFlipReg") || !file.exists(horizontalFlipReg))
			horizontalFlipReg=WriteFlipRegistration(infile,axis=flipAxis)
		# and flip all masks, applying the horiz bridging registration
		ReformatImage(infile,target=infile,
			registrations=c(FlipBridgingReg,horizontalFlipReg),
			filesToIgnoreModTimes=horizontalFlipReg, OverWrite='update', Verbose=T,
			output=flippedresampledfile,reformatoptions="-v --pad-out 0 --nn",dryrun=FALSE)
		# and OR ing those results
		Nrrd2op(c(infile,flippedresampledfile),orfile,'max',gzip=gzip)
	}
	unlink(horizontalFlipReg)
}

.makeReformatxTargetSpecification<-function(target){
	# a little function to make a target volume specification that 
	# CMTK's reformatx can understand
	if(is.character(target)){
		# assume this is a file
		# TODO: Intercept Amiramesh files and extract bounding box
		target=shQuote(target)
	} else if(is.vector(target)){
		# specify a target range c(Nx,Ny,Nz,dX,dY,dZ,[Ox,Oy,Oz])
		if(length(target)==9) {
			target=paste("--target-grid",
				paste(paste(target[1:3],collapse=","),paste(target[4:6],collapse=","),
				paste(target[7:9],collapse=","),sep=":"))
		} else if(length(target)==6) {
			target=paste("--target-grid",
				paste(paste(target[1:3],collapse=","),paste(target[4:6],collapse=","),sep=":"))
		} else stop("Incorrect target specification: ",target)
	} else {
		# can also give a density object
		# --target-grid
		#           Define target grid for reformating as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz]
		#           (dims:pixel:origin)
		# TODO: Double check definition of origin
		dims=dim(target)
		bb=getBoundingBox(target)
		vd=voxdim.gjdens(target)
		target=paste("--target-grid",shQuote(
			paste(paste(dims,collapse=","),paste(vd,collapse=","),paste(bb[c(1,3,5)],collapse=","),sep=":")
			))
	}
	target
}

ReformatImage<-function(floating,target,registrations,output, 
	dryrun=FALSE, Verbose=TRUE, MakeLock=TRUE, OverWrite=c("no","update","yes"),
	filesToIgnoreModTimes=NULL,
	reformatxPath="/usr/local/bin/reformatx",reformatoptions="-v --pad-out 0",...){
		# TODO improve default ouput file name
	if(missing(output)){
		output=file.path(dirname(floating),paste(basename(target),"-",basename(floating),'.nrrd',sep=""))
	} else if(isTRUE(file.info(output)$isdir)){
		output=file.path(output,paste(basename(target),"-",basename(floating),'.nrrd',sep=""))
	}
	if(is.logical(OverWrite)) OverWrite=ifelse(OverWrite,"yes","no")
	else OverWrite=match.arg(OverWrite)
	
	targetspec=.makeReformatxTargetSpecification(target)
	allinputs=c(floating,registrations)
	# if the target was a plain file add it to the inputs
	if(substring(targetspec,1,2)!="--") allinputs=c(allinputs,target)
	
	inputsExist=file.exists(allinputs)
	if(!all(inputsExist)){
		cat("Missing input files",basename(allinputs)[!inputsExist],"\n")
		return(FALSE)
	}
	if( file.exists(output) ){
		# output exists
		if(OverWrite=="no"){
			if(Verbose) cat("Output",output,"already exists; use OverWrite=\"yes\"\n")
			return(FALSE)
		} else if(OverWrite=="update"){
			# check modification times
			filesToCheck=setdiff(allinputs,filesToIgnoreModTimes)
		} else if(Verbose) cat("Overwriting",output,"because OverWrite=\"yes\"\n")
	} else OverWrite="yes" # just for the purpose of the runtime checks below 
		
	cmd=paste(shQuote(reformatxPath), reformatoptions,
		"-o",shQuote(output),"--floating",shQuote(floating),targetspec,
		paste(shQuote(registrations),collapse=" "))
	lockfile=paste(output,".lock",sep="")
	PrintCommand<-FALSE
	if(dryrun) PrintCommand<-TRUE
	if(!dryrun) {
		if(!MakeLock) system(cmd,...)
		else if(makelock(lockfile)){
			if(OverWrite=="update")
				PrintCommand<-RunCmdForNewerInput(cmd,filesToCheck,output,Verbose=Verbose,...)
			else {
				PrintCommand<-TRUE;system(cmd,...)
			}
			removelock(lockfile)
		} else if(Verbose) cat("Unable to make lockfile:",lockfile,"\n")
	}
	if(PrintCommand) cat("cmd:\n",cmd,"\n") 
	return(TRUE)
}

WriteFlipRegistration<-function(examplenrrd=nrrdfile,
	regfolder=file.path(tempdir(),"hflip.list"),axis=c("X","Y","Z"),...){
	axis=match.arg(axis)
	h=ReadNrrdHeader(examplenrrd)
	if(any(names(h)=='space directions'))
	voxdims=sqrt(rowSums(h[['space directions']]^2))
	else voxdims=h$spacings
	boundingBox=voxdims * (h$sizes - 1)
	if(length(boundingBox)!=3 || any(is.na(boundingBox)))
		stop("Unable to extract sensible bounding box")
	names(boundingBox)=c("X","Y","Z")

	hflipreg=structure(c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0), 
	.Dim = c(5L, 3L), .Dimnames = list(c("xlate", "rotate", "scale", "shear", "center"),
	 c("X", "Y", "Z")))
	hflipreg['scale',axis]=-1
	hflipreg['xlate',axis]=boundingBox[axis]

	WriteIGSRegistrationFolder(hflipreg,regfolder)
	return(regfolder)
}

WriteIdentityRegistration<-function(regfolder=file.path(tempdir(),"identityreg.list"),...){
	ireg=structure(c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0), 
		.Dim = c(5L, 3L), .Dimnames = list(c("xlate", "rotate", "scale", "shear", "center"),
		 c("X", "Y", "Z")))
	WriteIGSRegistrationFolder(ireg,regfolder)
	return(regfolder)
}

AutoCropNrrd<-function(infile, threshold=1,suffix="-acrop",
	outfile=NULL,outdir=NULL, options="")
{
	if(is.null(outfile) && is.null(outdir))
		outfile=sub("(\\.[^.]+)$",paste(suffix,"\\1",sep=""),infile)
	else if (is.null(outfile) && !is.null(outdir))
		outfile=file.path(outdir,sub("(\\.[^.]+)$",paste(suffix,"\\1",sep=""),basename(infile)))
	else if (!is.null(outfile)  && !is.null(outdir))
		outfile=file.path(outdir,outfile)
	
	# take a nrrd image and run Torsten's auto crop function
	# read in the resultant affine transformation file
	# and shift the nrrd's origin
	cropxformreg=paste(tempfile(),".list",sep="")
	options=paste('--auto-crop',threshold,'--crop-xform-out',shQuote(cropxformreg),options)
	
	tmpoutfile=paste(sep=".",outfile,"tmp.nrrd")
	cmd=paste("convert",options,shQuote(infile),shQuote(tmpoutfile))
	if(!RunCmdForNewerInput(cmd,infile,tmpoutfile)) return (FALSE)

	reg=ReadIGSRegistration(cropxformreg)
	xlate=reg$affine_xform$xlate
	unlink(cropxformreg)

	inh=ReadNrrdHeader(infile)
	outh=ReadNrrdHeader(tmpoutfile)
	originalOrigin=c(0,0,0)
	if("space origin"%in%names(inh)) originalOrigin=inh[["space origin"]]
	newOrigin=originalOrigin+xlate
	newOriginLine=paste("space origin: (",paste(newOrigin,collapse=","),")",sep="")
	
	oht=attr(outh,"headertext")
	if("space origin"%in%names(outh)){
		# replace existing space origin
		oht=sub("space origin: .*",newOriginLine,oht)
	} else {
		# just append
		oht=c(oht,newOriginLine)
	}
	# add a blank line
	oht=c(oht,"")
	tmpheader=tempfile()
	writeLines(oht,tmpheader)
	system(paste("unu data",shQuote(tmpoutfile),"| cat",tmpheader,"- >",shQuote(outfile)))
	system(cmd)
	unlink(c(tmpoutfile,tmpheader))
}