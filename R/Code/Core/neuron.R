# so you can just do:
# plot(ANeuron)
plot.neuron<-function(...) plotneuron2d(...)

# This function can be called directly or using the registration above
# it will be called as a fall-through by nat::read.neuron
read.neuron.extra<-function(f, ...){
	# generic function to read in neuron from any kind of file we know about
	# should return exactly one neuron on success
	if(!file.exists(f)) stop("Unable to read file: ",f)
	ext=tolower(sub(".*\\.([^.]+$)","\\1",basename(f)))

	if(ext=="asc") {
		n=ReadNeuronFromAsc(f, ...)
	} else {
		h=readLines(f,1) # nb readLines can cope with gzipped data
		
		if(regexpr("^;",h)>0)
			n=ReadNeuronFromAsc(f, ...)
		else stop("Unable to identify file format of file: ",f)
	}
	# we can normally rely on dotprops objects to have the correct class
	if(is.neuron(n,Strict=FALSE) && !is.dotprops(n)) as.neuron(n)
	else n
}

# register a catch all function to read any file types that aren't handled by
# other readers.
# NB1 that calling the format zzz ensures that it will be the last format to be
# called (since they are processed in alphabetical order). 
# NB2 note that setting magic to a function that swallows all arguments and
# always returns TRUE captures all remaining file types
nat::registerformat("zzz", read="read.neuron.extra", magic=function(...) TRUE,
  class='neuron')

#' [Deprecated] Write out a neuron in any of the file formats we know about 
#'
#' This function is being retained only for the ability to write neurolucida and
#' Borst format neurons.
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
write.neuron.extra<-function(n,filename=NULL,dir=NULL,ftype=c('neurolucida.asc','borst'),suffix=NULL,...){
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
  if(ftype=='neurolucida.asc'){
    WriteAscFromNeuron(n,filename,...)
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
#' @param x An object
#' @param reg matrix, path or function describing registration to apply
#' @param na.action What to do if a point can't be transformed
#' @param FallBackToAffine Whether to try calculating affine for a CMTK 
#'   registration if the warp fails.
#' @return transformed neuron
#' @export
#' @seealso \code{\link{transform.dotprops}}
transform.points3d<-function(x,reg,na.action=c('warn','drop','error'),FallBackToAffine=FALSE,...)
{
	.Deprecated('nat::xform')
	xform(x,reg,na.action=na.action,FallBackToAffine=FallBackToAffine,...)
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
	.Deprecated('nat::xform')
	xform(x,reg,na.action=na.action,FallBackToAffine=FallBackToAffine,...)
}
