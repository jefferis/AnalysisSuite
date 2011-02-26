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
	# landmark registration 
	# xform maps points in target to points in src
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
