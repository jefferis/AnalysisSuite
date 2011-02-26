# Some functions to help use Daniel Rueckert's 
# Image Registration Toolkit (IRTK)
# See: http://www.doc.ic.ac.uk/~dr/software

# See also VTKIO.R for functions to read the VTK pointset format
# used by IRTK

irtk.dof2mat<-function(doffile,matfile,Invert=FALSE,...){
	# see http://www.doc.ic.ac.uk/~dr/software/usage.html#dof2mat
	# reads in a transformation matrix from doffile and converts
	# to 4 x4 homogeneous transformation matrix
	deleteoutput=FALSE
	if(missing(matfile)) {
		matfile=tempfile()
		deleteoutput=TRUE
	}
	
	rval=.callirtk('dof2mat',
		c(shQuote(doffile),"-matout",shQuote(matfile),ifelse(Invert,character(0),"-invert")),...)
	if(rval>0) return(NULL)
	
	mat=scan(matfile,skip=1)
	if(length(mat)!=16) stop("Unable to read a 4x4 matrix. Output was:",readLines(matfile))
	if(deleteoutput) unlink(matfile)
	matrix(mat,nrow=4,byrow=T)
}

.callirtk<-function(cmd,args,irtkdir=path.expand("~/dev/registration/irtk"),Verbose=FALSE,DryRun=FALSE,...){
	cmdline=paste(file.path(irtkdir,cmd),paste(args,collapse=" "))
	if(Verbose || DryRun) print(cmdline)
	if(DryRun) return(FALSE)
	system(cmdline,...)
}

irtk.preg<-function(src, target=NULL, dofout=NULL, dofin=NULL,
	xformtype=c("rigid","affine","nonrigid"),cpspacing=10,...){
	# landmarks registration: xform maps points in target to points in src
	# 
	# this xform can later be used via "transformation" to map
	# the target image into the src space.
	# 
	# cpspacing is the control point spacing
		
	xformtype=match.arg(xformtype)
	cmd=paste("p",sep=substring(xformtype,1,1),"reg")
	landmarks=NULL
	if(is.list(src)){
		# this should be a landmark pair
		landmarks=src
		target=tempfile()
		src=tempfile()
		on.exit(unlink(c(src,target)))
		WriteVTKLandmarks(src,landmarks[[1]],"Landmark Set 1 (Source)")
		WriteVTKLandmarks(target,landmarks[[2]],"Landmark Set 2 (Target)")
	} else {
		if(is.null(target)) stop("Please specify target file")
		if(!file.exists(src)) stop("Missing src file: ",src)
		if(!file.exists(target)) stop("Missing target file: ",target)
	}

	if(is.null(dofout)) {
		if(!is.null(landmarks)){ 
			# ie we were given an R list not some files
			stop("Please supply an output file")
		} else {
			# constructing default output file based on srcfilename
			srcstem=sub("\\.[^.]+$","",basename(src))
			targetstem=sub("\\.[^.]+$","",target)
			dofout=paste(targetstem,"_",srcstem,"_",cmd,".dof",sep="")
		}
	}
		
	if(!is.null(dofin)) {
		if(!file.exists(dofin))
			stop("dofin file: ",dofin," is missing")
	}
	
	args=c(shQuote(target), shQuote(src), "-dofout",shQuote(dofout))
	if(xformtype=="nonrigid")
		args=c(args,"-ds",cpspacing)
		
	if(!is.null(dofin))
		args=c( args, paste("-dofin", shQuote(dofin)) )
	rval=.callirtk(cmd,args, ...)
	if(rval!=0) stop("error ",rval," in IRTK ",cmd)
	return(dofout)
}

irtk.transformation<-function(src, dofin, output, target,Invert=FALSE,
	datatype=c("points","image","surface"),
	interp=c("nn","linear","bspline","cspline","sinc"), moreArgs, ...){
	# See http://www.doc.ic.ac.uk/~dr/software/usage.html#transformation
	# TODO process moreArgs (direct options for image transformation)
	interp=match.arg(interp)
	datatype=match.arg(datatype)
	cmd=switch(datatype,
		surface='stransformation',points='ptransformation','transformation')
	
	if(datatype=="points" && !is.character(src)){
		# handle points provided directly as R matrix
		pointsrc=tempfile()
		WriteVTKLandmarks(pointsrc,src)
		src=pointsrc
		if(missing(output)) output=NA
		on.exit(unlink(pointsrc))
	} else if(!file.exists(src))
		stop("Cannot read source file: ",src)
	
	ReadOutputPoints=FALSE
	if(missing(output)) {
		srcstem=sub("\\.[^.]+$","",src)
		srcext=sub(".*(\\.[^.]+)$","\\1",basename(src))
		dofinstem=sub("^([^._]+)[._].*","\\1",basename(dofin))
		dofout=paste(srcstem,"-",dofinstem,srcext,sep="")
	} else if(is.na(output)){
		if(datatype=="points"){
			# assume that we want to get the points straight back into memory
			ReadOutputPoints=TRUE
			pointsout=tempfile()
			output=pointsout
			on.exit(unlink(pointsout),add=TRUE)
		} else 
			stop("Have not implemented direct reading of (s)transformation output")
	}
	if(file.access(dirname(output),2)!=0)
		stop("Cannot write to output directory:",dirname(output))

	if(!file.exists(src)) stop("Cannot read dofin file: ",dofin)

	args=c(shQuote(src),shQuote(output),"-dofin",shQuote(dofin))

	if(Invert) argcs=c(args,"-invert")

	if(datatype=="image"){
		if(interp!="nn")
			args=c(args,paste("-",interp,sep=""))
		if(!missing(target))
			args=c("-target",target)
	}

	rval=.callirtk(cmd, args, ...)
	if(rval!=0) stop("error ",rval," in IRTK ",cmd)
	if(ReadOutputPoints) ReadVTKLandmarks(output)
	else output
}
