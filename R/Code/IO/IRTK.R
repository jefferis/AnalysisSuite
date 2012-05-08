# Some functions to help use Daniel Rueckert's 
# Image Registration Toolkit (IRTK)
# See: http://www.doc.ic.ac.uk/~dr/software

# See also VTKIO.R for functions to read the VTK pointset format
# used by IRTK

# General notes:
# Transformation Direction
# ========================
# (from rreg)
# Estimate a rigid transformation between two images. The rigid transformation
# is represented by six degrees of freedom when registering 3-D images. Three
# degrees of freedom encode rotations about the axes and the remainder encode
# components of a translation. Three degrees of freedom can encode a rigid
# transformation between 2-D images, one for a rotation about the origin and two
# for a translation. The mandatory arguments are the file names of two images, the
# first is treated as a `target' image and the second is treated as a `source'
# image. The transformation that is estimated maps locations in the target to
# locations in the source.
# ie src <-- target

# So like CMTK, the affine registration calculated is the inverse of what one
# might naively expect.  The reason for this is that non-rigid transforms are
# more efficient if they provide a lookup from a regulat grid in target space
# back to irregular locations in the source image (which can be interpolated)

irtk.help<-function(cmd,web=FALSE){
	# simple wrapper function that returns command line usage information for 
	# irtk commands
	if(web) browseURL(paste("http://www.doc.ic.ac.uk/~dr/software/usage.html",sep="#",cmd))
	else .callirtk(cmd,"")
}

irtk.readmat<-function(matfile,endian='big'){
	# read an irtk format binary matrix file
	# this is not documented but seems fairly obvious
	# However I do wonder whether always big endian
	# OK - Not quite so obvious - the binary mat file
	# stores the inverse of what is written to stdout
	# not sure which corresponds to version discussed in documentation
	con=file(matfile,open='rb')
	on.exit(close(con))
	header=readLines(con,1) 
	# irtkMatrix 4 x 4
	if(regexpr("^irtkMatrix",header)<0) stop("This is not a valid IRTK matrix")
	dims=as.integer(unlist(strsplit(
		sub(".* (\\d+) x (\\d+)","\\1 \\2",header)," ")))
	m=readBin(con,what="double",n=prod(dims),size=8,endian=endian)
	dim(m)=dims
	m
}

irtk.dof2mat<-function(doffile,matfile,Invert=FALSE,...){
	# see http://www.doc.ic.ac.uk/~dr/software/usage.html#dof2mat
	# reads in a transformation matrix from doffile and converts
	# to 4 x4 homogeneous transformation matrix
	if(missing(matfile)) {
		matfile=tempfile()
		on.exit(unlink(matfile))
	}
	
	rval=.callirtk('dof2mat',
		c(shQuote(doffile),"-matout",shQuote(matfile),ifelse(Invert,"-invert","")),...)
	if(rval>0) return(NULL)
	
	irtk.readmat(matfile)
}

irtk.dofinvert<-function(dofin,dofout,...){
	# see http://www.doc.ic.ac.uk/~dr/software/usage.html#dofinvert
	# inverts an affine IRTK format transformation 
	if(missing(dofout)) dofout=sub("(\\.[^.]+)$","-inv\\1",dofin)
	rval=.callirtk('dofinvert',
		c(shQuote(dofin),shQuote(dofout)),...)
	if(rval>0) return(NULL)
	else return(dofout)
}

.callirtk<-function(cmd,args,irtkdir=path.expand("~/dev/registration/irtk"),Verbose=FALSE,DryRun=FALSE,...){
	cmdline=paste(file.path(irtkdir,cmd),paste(args,collapse=" "))
	if(Verbose || DryRun) print(cmdline)
	if(DryRun) return(FALSE)
	system(cmdline,...)
}

irtk.preg<-function(target, src=NULL, dofout=NULL, dofin=NULL,
	xformtype=c("rigid","affine","nonrigid"),cpspacing=10,...){
	# landmarks registration: xform maps points in src to points in target
	# ie 2->1. This is the opposite of what the IRTK docs say, but I am 
	# quite convinced. Specifically:
	# I think that affine finds a 9dof xform from 2->1
	# Then if you do dof2mat you get a binary irtk matrix file
	# containing the inverse of the 9dof xform (ie 1->2). 
	# Just to add to the confusion dof2mat will show the 9dof xform on stdout.
	# ie 2->1
	# BUT whatever is happening internally both mat and dof files specify 
	# a 1->2 transformation so the documentation is fundamentally correct
	
	# cpspacing is the control point spacing
		
	xformtype=match.arg(xformtype)
	cmd=paste("p",sep=substring(xformtype,1,1),"reg")
	landmarks=NULL
	if(is.list(target)){
		# this should be a landmark pair
		landmarks=target
		target=tempfile()
		src=tempfile()
		on.exit(unlink(c(src,target)))
		WriteVTKLandmarks(target,landmarks[[1]],"Landmark Set 1 (Target)")
		WriteVTKLandmarks(src,landmarks[[2]],"Landmark Set 2 (Source)")
	} else {
		if(is.null(src)) stop("Please specify source file")
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
			dofout=paste(targetstem,"-",srcstem,"_",cmd,".dof",sep="")
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
# Quoting from help for transformation:
# Transform one image (the source) onto the voxel lattice of a second image (the
# target) using a transformation estimate.The transformation estimate must be
# 

# A typical set of command line arguments
# might be as follows: transformation a.nii.gz out.nii.gz -dofin tr-a-b.dof
# -target b.nii.gz -linear In this example the image being transformed is
# a.nii.gz, the voxel lattice of the resulting image (out.nii.gz) will match that
# of the target image (b.nii.gz). The transformation given, tr-a-b.dof, is used to
# transform the source intensities to the target but it should be noted that this
# transformation maps locations in the target image to locations in the source
# image. The intensity at each voxel of out.nii.gz is `pulled-back' from the
# corresponding source location. If no file is given for the transformation, a
# default identity transformation is used. I.e. it is assumed that the target and
# source images share the same world coordinate system (although they need not
# share the same voxel lattice). An interpolation method (`-linear' in the above
# example) is usually required because the source locations corresponding to
# target voxels may not coincide with the source voxel lattice. If no
# interpolation scheme is specified, nearest neighbour interpolation is used by
# default. If no target image is specified, the voxel lattice of the source is
# used as a default.
	
	
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
		output=paste(srcstem,"-",dofinstem,srcext,sep="")
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

	if(Invert) args=c(args,"-invert")

	if(datatype=="image"){
		if(interp!="nn")
			args=c(args,paste("-",interp,sep=""))
		if(!missing(target))
			args=c(args,"-target",target)
	}

	rval=.callirtk(cmd, args, ...)
	if(rval!=0) stop("error ",rval," in IRTK ",cmd)
	if(ReadOutputPoints) ReadVTKLandmarks(output)
	else output
}

irtk.transformedPoints<-function(xyzs=NULL,dofins,direction=c("inverse","forward"),
	transforms=c("warp","affine"),...){
	# provide a wrapper for irtk.transformation suitable for calling
	# from TransformNeuron
	# I don't know how to access the affine transform preceding a non-rigid transform
	# so right now the transforms switch is ignored
	
	transforms=match.arg(transforms,several.ok=TRUE)

	direction=match.arg(direction) #nb inverse implies from sample to ref

	# massage xyzs input to a 3 col matrix
	if(is.data.frame(xyzs)) xyzs=data.matrix(xyzs)
		if(ncol(xyzs)>3){
		if(all(c("X","Y","Z")) %in% colnames(xyzs))
		xyzs=xyzs[,c("X","Y","Z")]
		else xyzs=xyzs[,1:3]
	}

	l=list(pre=xyzs)
	if(length(transforms)!=length(dofins)) 
		stop("Must provide the same number of transform descriptors as dof files")
	names(dofins)=transforms
	for (t in transforms){
		l[[t]]=irtk.transformation(xyzs,dofins[t],output=NA,
			Invert=ifelse(direction=="forward",TRUE,FALSE))
	}
	l
}

irtk.TransformNeuron<-function(neuron,dofin=NULL,transform=c("warp","affine"),...){
	# FIXME - transform specifiers (warp, affine etc) are irrelevant for IRTK 
	transform=match.arg(transform,several.ok=FALSE)
	xyzmatrix(neuron)<-irtk.transformedPoints(xyzmatrix(neuron),
		dofins=dofin,transforms=transform,...)[[transform]]
	neuron
}
