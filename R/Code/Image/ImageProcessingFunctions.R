# General Image Processing Functions combining CMTK / UNU or my R Native Code
# distinct from ImageAnalyisFunctions which are more statistical in nature
# this is really about munging files

ResampleAndFlipMasks<-function(masks,outdir,FlipBridgingReg,flipAxis=c("X","Y","Z"),targetspec,
	suffix="-resampled",gzip=TRUE){
	flipAxis=match.arg(flipAxis)
	if(!file.exists(outdir)) dir.create(outdir)
	resampledfiles=ResampleMasks(masks,outdir,targetspec,suffix=suffix)
	FlipAndORMasks(resampledfiles,outdir,FlipBridgingReg,flipAxis,gzip)
}

ResampleMasks<-function(masks,...) ResampleImages(images=masks,...,interpolation="nn")

ResampleImages<-function(images,outdir,targetspec,registrations,suffix="-resampled",
	interpolation=c("linear","nn","cubic"),TargetIsMask=FALSE,Verbose=TRUE,...){
	interpolation=match.arg(interpolation)
	if(!file.exists(outdir)) dir.create(outdir)
	resampledfiles=file.path(outdir,
		sub(".nrrd$",paste(suffix,".nrrd",sep=""),basename(images)))
	useIdentity=missing(registrations)
	if(useIdentity)
		registrations=WriteIdentityRegistration()
	# set additional reformat options
	# TargetIsMask means only calculate for pixels where mask is non-zero
	reformatoptions=paste(ifelse(TargetIsMask,"--mask","")," -v --pad-out 0 --",sep="",interpolation)
	filesToIgnoreModTimes=if(useIdentity) registrations else NULL
	for (i in seq(images)){		
		# resample - use reformatx to do this, to ensure that we get the same result
		# even if we have to provide an identity registration
		ReformatImage(images[i],target=targetspec, registrations=registrations,
			filesToIgnoreModTimes=filesToIgnoreModTimes,
			OverWrite='update',output=resampledfiles[i],
			reformatoptions=reformatoptions,dryrun=FALSE,Verbose=Verbose,...)
	}
	if(useIdentity) unlink(registrations)
	invisible(resampledfiles)
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
	reformatxPath="reformatx",reformatoptions="-v --pad-out 0",
	Push=FALSE,...){
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
		"-o",shQuote(output),ifelse(Push,"--push",""),"--floating",shQuote(floating),targetspec,
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
	outfile=NULL,outdir=NULL, options="",convertTool="convertx",Force=FALSE,UseLock=FALSE)
{
	if(is.null(outfile) && is.null(outdir))
		outfile=sub("(\\.[^.]+)$",paste(suffix,"\\1",sep=""),infile)
	else if (is.null(outfile) && !is.null(outdir))
		outfile=file.path(outdir,sub("(\\.[^.]+)$",paste(suffix,"\\1",sep=""),basename(infile)))
	else if (!is.null(outfile)  && !is.null(outdir))
		outfile=file.path(outdir,outfile)

	inh=ReadNrrdHeader(infile)
	hasOrigin=TRUE
	if(is.null(inh[["space origin"]])) hasOrigin=FALSE
	
	# Check if we actually need to run
	if(!Force && !RunCmdForNewerInput(NULL,infile,outfile)) return (FALSE)
	lockfile=paste(outfile,sep=".","lock")
	if(UseLock && !makelock(lockfile)) return (FALSE)
	
	# take a nrrd image and run Torsten's auto crop function
	if(hasOrigin){
		options=paste("--crop-by-threshold",threshold,options)
		# Torsten's tool will add to existing origin if present
		cmd=paste(convertTool,options,shQuote(infile),shQuote(outfile))
		system(cmd)
		if(UseLock) unlink(lockfile)
		return (TRUE)
	} else {
		# read in the resultant affine transformation file
		# and shift the nrrd's origin assuming that it was 0,0,0
		cropxformreg=paste(tempfile(),".list",sep="")
		options=paste("--crop-by-threshold-write-xform",threshold,options)
		tmpoutfile=paste(sep=".",outfile,"tmp.nrrd")
		# nb xform is written to stdout by default
		cmd=paste(convertTool,options,shQuote(infile),shQuote(tmpoutfile),
			">",shQuote(cropxformreg))
		system(cmd)
	}

	reg=ReadIGSRegistration(cropxformreg)
	xlate=reg$affine_xform$xlate
	unlink(cropxformreg)

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
	unlink(c(tmpoutfile,tmpheader))
	if(UseLock) unlink(lockfile)
}

#' Quantise and smooth a nrrd, by default producing (0,1) range float image
#'
#' Details
#' @param summary of each parameter
#' @return return value
#' @export
#' @seealso \code{\link{somefun}}
#' @examples
NormaliseAndSmoothNrrd<-function(infile,outfile,outdir,threshold,max,
	out_type=c("float","uint8", "uint16", "uint32"),
	centering=c("node","cell"), scalefactor="x1 x1 x1",sigma=3,kernelsigmacutoff=2.5,
	DryRun=FALSE, gzip=FALSE,UseLock=FALSE)
{
	if(missing(outfile) && !missing(outdir))
		outfile=file.path(outdir,basename(infile))
	else if(missing(outdir) && missing(outfile))
		outfile=sub("(\\.[^.]+)$",paste("-nsmooth","\\1",sep=""),infile)
	centering=match.arg(centering)
	out_type=match.arg(out_type)
	haveRun=FALSE
	cmd=paste("unu 3op clamp ",threshold,infile,max)
	cmd=paste(cmd, "| unu 2op - - ",threshold)
	cmd=paste(cmd, "| unu 2op / - ",(max-threshold),"-t float")
	kernel=paste("--kernel gauss:",sigma,",",kernelsigmacutoff,sep="")
	cmd=paste(cmd, "| unu resample --size",scalefactor,kernel)
	cmd=paste(cmd,"--center",centering)
	lockfile=paste(outfile,sep=".","lock")
	if(UseLock && !makelock(lockfile)) return (haveRun)

	# nb 9f means max compression (9), specialised for filtered data
	if(out_type!='float') {
		bits=as.integer(sub("uint","",out_type))
		cmd=paste(cmd,"|unu quantize -b",bits,'-min',0,'-max',1)
	}
	if(gzip) cmd=paste(cmd,"| unu save --format nrrd --encoding gz:9f -o",outfile)
	else cmd=paste(cmd,"-o",outfile)
	if(DryRun) print(cmd)
	else haveRun<-RunCmdForNewerInput(cmd,infile,outfile)
	if(UseLock) unlink(lockfile)
	return(haveRun)
}

SimplifyLabelFile<-function(f,omitMaterials="CellBody",includeMaterials=NULL)
{
	d=ReadAmiramesh(f)
	m=attr(d,"Materials")
	m=subset(m,!rownames(m)%in%omitMaterials & !is.na(id))
	d[!d%in%m$level]=0
	d[d%in%m$level]=1
	d
}

#' Calculate image statistics for a nrrd or other CMTK compatible file
#'
#' When given a mask returns a table that includes the number of non-zero voxels
#' in the main image for every level of label in the mask image
#' @param f Path to image file (any CMTK compatible format)
#' @param mask Optional path to a mask file
#' @param masktype Whether mask should be treated as label field or binary mask
#'        (default label)
#' @return return dataframe describing results
#' @export
#' @examples
#' #CMTKStatistics('someneuron.nrrd',mask='neuropilregionmask.nrrd')
CMTKStatistics<-function(f,mask,masktype=c("label","binary"),exe="statistics"){
	masktype=match.arg(masktype)
	if(length(f)>1) return(sapply(f,CMTKStatistics,mask=mask,masktype=masktype,exe=exe))
	args=f
	if(!missing(mask)){
		args=c(ifelse(masktype=='label','--Mask','--mask'),mask,args)
	}
	rval=system2(exe,args,stdout=TRUE)
	# there is a bug in CMTK statistics when used with a mask
	# the header says that there are two entropy columns (H1,H2)
	# but in fact there is only 1
	rval[2]=sub('H1\tH2','H',rval[2])
	rval[2]=sub('#M','MaskLevel',rval[2])
	tc=textConnection(rval)
	on.exit(close(tc))
	read.table(tc,header=TRUE,skip=1,comment.char="")
}

CMTKSimilarity<-function(floating,reference,exe="similarity",simargs=NULL,...){
	if(length(floating)>1 || length(reference)>1){
		return(t(mapply(CMTKSimilarity,floating,reference,
			MoreArgs=list(exe=exe,simargs=simargs))))
	}	
	rval=system2(exe,c(simargs,reference,floating),stdout=TRUE,...)
	scorelines=grep("^SIM",rval,val=T)
	tc=textConnection(scorelines)
	on.exit(close(tc))
	t=read.table(tc,header=TRUE)
	# data.matrix(t[-1])
	unlist(t[-1])
}
