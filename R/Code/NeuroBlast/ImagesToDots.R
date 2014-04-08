# Functions to process images of neurons and generate more compact
# dot-based representations
# Based on Nick Masse's matlab code (esp extract_properties.m)

require(RANN)

#' Find the nearest point(s) in a neuron/dotprops to a 3d point
#'
#' For example see how far a marked cell body location is from a cell
#' @param dp the dotprops or neuron
#' @param point the point 
#' @param k, the number of nearest neighbours (default 1)
#' @param ReturnDistance, return distance as 4th column
#' @return the nearest point(s)
#' @export
#' @examples
#' nearestpoint(dps[['ChaMARCM-F000673_seg001']],c(255.883,141.301,90.7188))
nearestpoint<-function(dp,point,k=1){
	xyz=xyzmatrix(dp)
	if(is.null(dim(point))){
		if(length(point)==3) point=matrix(point,ncol=3)
		else stop("cannot understand passed point: ",point)
	}
	nnl=nn2(xyz,point,k=k)
	xyz[nnl$nn.idx,]
}

length.dotprops<-function(dp) nrow(dp$points)

# redefining length upsets str
str.dotprops<-function(dp,...) {class(dp)<-"list";str(dp,...)}

#' Deprecated: Carry out in memory digest of a dotprops object
#'
#' see nat::ndigest
#'
#' This takes care to remove non-essential attributes
#' @param dp Dotprops object
#' @param ... additional parameters passed to digest function
#' @return character vector of digest value
#' @export
#' @import digest
#' @seealso \code{\link{digest}}
#' @examples
#' digest(dps[[1]])
digest.dotprops<-function(dp,...){
  .Deprecated("nat::ndigest")
  nat::ndigest(as.dotprops(dp), ...)
}

# Set up a private function that directly accesses internal LAPACK routine of eigen
if(R.version$major<3){
  .internal_lars<-function(x) .Call("La_rs", x, only.values=FALSE, PACKAGE = "base")
} else {
  .internal_lars<-function(x) .Internal(La_rs(x, only.values=FALSE))
}

#' Deprecated method to calculate dotprops
#' see nat::dotprops
DotProperties<-function(points,k=20,UseLabels=TRUE,na.rm=FALSE){
	.Deprecated('dotprops','nat')
	nat::dotprops(x=points,k=k,Labels=UseLabels,na.rm=TRUE)
}

DotPropertiesFromFile<-function(f, xformfun=NULL, ...){
	ext=sub(".*\\.([^.]+$)","\\1",basename(f))
	l=list()
	if(ext=="nrrd"){
		x=Read3DDensityFromNrrd(f)
		l$points=ind2coord(x)
	} else if(ext=='csv') {
		l$points=read.csv(f,header=FALSE)
		colnames(l$points)=c("X","Y","Z")
	} else {
		# maybe some kind of neuron
		t=try(read.neuron(f))
		if(inherits(t,'try-error'))
			stop("Cannot extract dots from file type: ",ext)
		l$points=t$d[,c("X","Y","Z")]
	}
	if(!is.null(xformfun))
		l$points=xformfun(l$points)
	l=dotprops(l$points, ...)
	attr(l,'file')=f
	fi=file.info(f)
	attr(l,'mtime')=fi$mtime
	attr(l,'size')=fi$size
	l
}

#' Deprecated Transform dot property object using specified registration
#'
#' see nat::xform
#'
#' @param dp dotprops object to transform
#' @param reg Path to CMTK registration file OR function to transform points
#' @param k Number of neighbour points to use when recalculating dot properties
#' @param RecalculateDotProps Whether to recalculate tangent vector etc after 
#'   applying transformation
#' @param ... additional arguments passed to transformedPoints or reg function
#' @return return points
#' @export
#' @seealso \code{\link{transformedPoints}}
transform.dotprops<-function(dp,reg,k, RecalculateDotProps=T,na.action=c('warn','drop','error'),...) {
	na.action=match.arg(na.action)
	if(!RecalculateDotProps) stop("RecalculateDotProps must always be TRUE!")
	.Deprecated('xform','nat')
	nat::xform(dp,reg,k=k)
}

#' Compute point & tangent vector similarity score between two dotprops objects
#'
#' UseAlpha determines whether the alpha values (eig1-eig2)/sum(eig1:3)
#' are passed on to WeightedNNBasedLinesetMatching. These will be used to scale
#' the dot products of the direction vectors for nearest neighbour pairs.
#' @param dp1,dp2 dotprops objects
#' @param UseAlpha Whether to scale dot product of tangent vectors (default=F)
#' @return Return value of NNDistFun passd to WeightedNNBasedLinesetMatching
#' @export
#' @seealso \code{\link{DotProperties}},\code{\link{WeightedNNBasedLinesetMatching}}
#' @examples
WeightedNNBasedLinesetMatching.dotprops<-function(dp1,dp2,UseAlpha=FALSE,...){
	if(UseAlpha)
		WeightedNNBasedLinesetMatching(dp1$points,dp2$points,dvs1=dp1$vect,dvs2=dp2$vect,
			alphas1=dp1$alpha,alphas2=dp2$alpha,...)
	else 
		WeightedNNBasedLinesetMatching(dp1$points,dp2$points,dvs1=dp1$vect,dvs2=dp2$vect,...)
}


# stop()
# 
# points=ReadAmiramesh("/GD/projects/Nick/FruCloneClustering/data/SAKW13-1_gj_dimension_reduced.am")
