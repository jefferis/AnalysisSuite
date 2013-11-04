is.neuron<-function(n,Strict=FALSE) {
	# If Strict is FALSE will also return TRUE
	# if n is a list which looks like a neuron
	inherits(n,"neuron") ||
		(!Strict && is.list(n) && !is.null(n$SegList))
}

#' Arithmetic for neuron coordinates
#'
#' If x is one number or 4-vector, multiply xyz and diameter by that
#' If x is a 3-vector, multiply xyz only
#' TODO Figure out how to document arithemtic functions in one go
#' @param n a neuron
#' @param x (a numeric vector to multiply neuron coords in neuron)
#' @return modified neuron
#' @export
#' @examples
#' n1<-MyNeurons[[1]]*2
#' n2<-MyNeurons[[1]]*c(2,2,2,2)
#' stopifnot(all.equal(n1,n2))
#' n3<-MyNeurons[[1]]*c(2,2,4)
`*.neuron` <- function(n,x) {
	# TODO look into S3 generics for this functionality
	
	nd=n$d[,c("X","Y","Z","W")]
	stopifnot(is.numeric(x))
	lx=length(x)
	if(lx==1) nd[,-4]=nd[,-4]*x
	else if(lx%in%c(3,4)) nd[,1:lx]=t(t(nd[,1:lx])*x)
	else stop("expects a numeric vector of length 1, 3 or 4")
	n$d[,colnames(nd)]=nd
	n
}

`+.neuron` <- function(n,x) {
	if(!is.numeric(x))
		stop("expects a numeric vector")
	nd=n$d[,c("X","Y","Z","W")]
	lx=length(x)
	if(lx==1) nd[,-4]=nd[,-4]+x
	else if(lx%in%c(3,4)) nd[,1:lx]=t(t(nd[,1:lx])+x)
	else stop("expects a numeric vector of length 1, 3 or 4")
	n$d[,colnames(nd)]=nd
	n
}

`-.neuron` <- function(n,x) n+(-x)
`/.neuron` <- function(n,x) n*(1/x)

#' Divide neuron coords by a factor (and optionally center)
#'
#' Note that if scale=TRUE, the neuron will be rescaled to unit sd in each axis
#' likewise if center=TRUE, the neuron will be centred around the axis means
#' @param scale 3-vector used to divide x,y,z coords
#' @param center 3-vector to subtract from x,y,z coords
#' @return neuron with scaled coordinates
#' @export
#' @seealso \code{\link{scale.default}}
#' @examples
#' n1.scaledown=scale(MyNeurons[[1]],c(2,2,3))
#' n1.scaleup=scale(MyNeurons[[1]],1/c(2,2,3))
scale.neuron<-function(n,scale,center=F){
	d=xyzmatrix(n)
	ds=scale(d,scale=scale,center=center)
	n$d[,colnames(d)]=ds
	n
}

all.equal.neuron<-function(target,current,tolerance=1e-6,check.attributes=FALSE,
	fieldsToCheck=c("NeuronName", "NumPoints", "StartPoint", "BranchPoints",
		"EndPoints", "NumSegs", "SegList", "d"), fieldsToCheckIfPresent="nTrees",
	CheckSharedFieldsOnly=FALSE, ...){
	if(length(fieldsToCheck)==1 && is.na(fieldsToCheck))
		fieldsToCheck=names(current)
		
	if(!is.neuron(target) || !is.neuron(current))
		return ("target and current must both be neurons")
	fieldsInCommon=intersect(names(target),names(current))
	# figure out which of the optional fields to check are present
	fieldsToCheckIfPresent=intersect(fieldsInCommon,fieldsToCheckIfPresent)
	# and add those to the fields to check 
	fieldsToCheck=unique(c(fieldsToCheck,fieldsToCheckIfPresent))
	if(CheckSharedFieldsOnly){
		fieldsToCheck=intersect(fieldsInCommon,fieldsToCheck)
	} else{
		# check all core fields
		missingfields=setdiff(fieldsToCheck,names(current))
		if(length(missingfields)>0)
			return(paste("Current missing fields: ",missingfields))
		missingfields=setdiff(fieldsToCheck,names(target))
		if(length(missingfields)>0)
			return(paste("Target missing fields: ",missingfields))		
	}
	all.equal(target[fieldsToCheck],current[fieldsToCheck],
		tolerance=tolerance, check.attributes=check.attributes, ...)
}

as.neuron<-function(n){
	if(is.null(n)) return (NULL)
	if(!is.neuron(n,Strict=TRUE)) class(n)=c("neuron",class(n))
	n
}

# so you can just do:
# plot(ANeuron)
plot.neuron<-function(...) plotneuron2d(...)

read.neuron<-function(f, ...){
	# generic function to read in neuron from any kind of file we know about
	# should return exactly one neuron on success
	if(!file.exists(f)) stop("Unable to read file: ",f)
	ext=tolower(sub(".*\\.([^.]+$)","\\1",basename(f)))

	if(ext=="asc")
		n=ReadNeuronFromAsc(f, ...)
	else if(ext=="swc")
		n=ReadNeuronFromSWC(f, ...)
	else if(ext=="rds")
		n=readRDS(f)
	else if(ext=="rda"){
		objname=load(f,envir=environment())
		if(length(objname)>1) stop("More than 1 object in file:",f)
		n=get(objname,envir=environment())
	} else {
    h=readLines(f,1) # nb readLines can cope with gzipped data
    
		if(regexpr("amira",h,ignore.case=TRUE)>0){
      # check to see what kind of amiramesh neuron we have
      ftype=try(ReadAmiramesh.Header(f)$Parameters$ContentType[1])
      if(inherits(ftype,'try-error') || is.null(ftype))
        stop("Unable to indentify amiramesh neuron")
      else if(ftype=='HxLineSet')
        n=ReadNeuronFromAM(f, ...)
      else 
        n=ReadNeuronFromAM3D(f, ...)
    } else if(regexpr("xml",h,ignore.case=TRUE)>0)
			n=ReadNeuronsFromLongairTraces(f, MergePaths=TRUE, ...)
		else if(regexpr("^;",h)>0)
			n=ReadNeuronFromAsc(f, ...)
		else stop("Unable to identify file format of file: ",f)
	}
	# we can normally rely on dotprops objects to have the correct class
	if(is.neuron(n,Strict=FALSE) && !is.dotprops(n)) as.neuron(n)
	else n
}

#' Write out a neuron in any of the file formats we know about
#'
#' If filename is not specified the neuron's InputFileName field will be checked.
#' If this is missing there will be an error.
#' If dir is specified it will be combined with basename(filename).
#' If filename is specified but ftype is not, it will be inferred from filename's
#'   extension.
#' @param n A neuron
#' @param filename Path to output file
#' @param dir Path to directory (this will replace dirname(filename) if specified)
#' @param ftype File type (a unique abbreviation of 
#'   'swc','lineset.am','skeletonize.am','neurolucida.asc','borst','rds')
#' @param suffix Will replace the default suffix for the filetype and should
#'   include the period eg suffix='.amiramesh' or suffix='_reg.swc'
#' @param ... Additional arguments passed to selected  WriteNeuron function
#' @return return value
#' @export
#' @seealso \code{\link{WriteSWCFile, WriteNeuronToAM, WriteNeuronToAM3D, 
#'   WriteAscFromNeuron, WriteBorstFile,saveRDS}}
write.neuron<-function(n,filename=NULL,dir=NULL,ftype=c('swc','lineset.am',
    'skeletonize.am','neurolucida.asc','borst','rds'),suffix=NULL,...){
  if(is.dotprops(n)){
    # we only know how to save dotprops objects in R's internal format
    ftype='rds'
  }
  if(is.null(filename)){
    # no filename was specified - use the one embedded in neuron
    filename=n$InputFileName
    if(is.null(filename)) stop("No filename specified and neuron does not have an InputFileName")
    if(!is.null(suffix)){
      # we specified an explicit suffix - use this to identify file type
      ftype_from_ext=switch(tolower(suffix),.swc='swc',.asc='neurolucida.asc',
        .am='lineset.am',.amiramesh='lineset.am',.borst='borst',.rds='rds',NA)
      if(is.na(ftype_from_ext) && length(ftype!=1)){
        stop("file suffix: ",suffix,
          " does not uniquely identify filetype nor has this been specified directly")
      } else if(length(ftype_from_ext)>1)
        ftype=ftype_from_ext
      else
        ftype=match.arg(ftype)
    } else {
      # use the file type to specify the suffix
      ftype=match.arg(ftype)
      suffix=sub(".*?([^.]+)$",".\\1",ftype)
    }
  } else {
    ext=sub(".*(\\.[^.]+)$","\\1",filename)
    ftype_from_ext=switch(tolower(ext),.swc='swc',.asc='neurolucida.asc',
      .am='lineset.am',.amiramesh='lineset.am',.borst='borst',.rds='rds',NA)
    if(!is.na(ftype_from_ext) && length(ftype)!=1)
      ftype=ftype_from_ext
    else ftype=match.arg(ftype)
  }
  # replace with an explicit suffix if desired
  if(!is.null(suffix)) filename=sub("(\\.[^.]+)$",suffix,filename)
  if(!is.null(dir)){
    filename=file.path(dir,basename(filename))
  }
  if(ftype=='rds'){
    saveRDS(n,file=filename,...)
  } else if(ftype=='lineset.am'){
    WriteNeuronToAM(n,filename,...)
  } else if(ftype=='skeletonize.am'){
    WriteNeuronToAM3D(n,filename,...)
  } else if(ftype=='neurolucida.asc'){
    WriteAscFromNeuron(n,filename,...)
  } else if(ftype=='swc'){
    WriteSWCFile(n,filename,...)
  } else if(ftype=='borst'){
    WriteBorstFile(n,filename,...)
  } else {
    stop("Unimplemented file type ",ftype)
  }
}
#' transform 3d points using a function, an affine matrix, or CMTK registration
#'
#' If a CMTK registration file is supplied the warping registration will be
#' applied.
#' @param x A neuron
#' @param reg matrix, path or function describing registration to apply
#' @param na.action What to do if a point can't be transformed
#' @param FallBackToAffine Whether to try calculating affine for a CMTK 
#'   registration if the warp fails.
#' @return transformed neuron
#' @export
#' @seealso \code{\link{transform.dotprops}}
transform.points3d<-function(x,reg,na.action=c('warn','drop','error'),FallBackToAffine=FALSE,...)
{
  na.action=match.arg(na.action)
  xyz=xyzmatrix(x)
  if(is.function(reg)){
  	# we've been given a function - apply this to points
  	pointst=reg(xyz)
  } else if(is.matrix(reg)) {
  	# we've been given an affine transformation matrix
  	pointst=TransformPoints(xyz,reg)
  } else {
  	# assume we've been given a CMTK registration file
  	pointst=transformedPoints(xyzs=xyz,warpfile=reg,transforms='warp',...)[['warp']]
  	naPoints=is.na(pointst[,1])
  	if(any(naPoints) && FallBackToAffine){
  		affpoints=transformedPoints(xyzs=xyz[naPoints,],warpfile=reg,
  			transforms='affine',...)[['affine']]
  		pointst[naPoints,]=affpoints
  	}
  }
  naPoints=is.na(pointst[,1])
  if(any(naPoints)){
    if(na.action=='drop') pointst=pointst[!naPoints,]
    else if(na.action=="warn")
      warning("There were ",length(naPoints),' that could not be transformed')
    else if(na.action=="error")
      stop("There were ",length(naPoints),' that could not be transformed')
  }
  xyzmatrix(x)<-pointst
  x
}

#' transform a neuron using a function, an affine matrix, or CMTK registration
#'
#' If a CMTK registration file is supplied the warping registration will be
#' applied.
#' @param x A neuron
#' @param reg matrix, path or function describing registration to apply
#' @param na.action What to do if a point can't be transformed
#' @param FallBackToAffine Whether to try calculating affine for a CMTK 
#'   registration if the warp fails.
#' @return transformed neuron
#' @export
#' @seealso \code{\link{transform.dotprops}}
transform.neuron<-function(x,reg,na.action='error',FallBackToAffine=TRUE,...) {
  transform.points3d(x=x,reg=reg,na.action=na.action,FallBackToAffine=FallBackToAffine,...)
}
