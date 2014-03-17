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

ind2coord<-function(inds, ...) UseMethod("ind2coord")

ind2coord.array<-function(inds, voxdims, origin, ...){
	dims=dim(inds)
	# can supply a gjdens object as 2nd param - which will then
	# provide origin, voxdims
	if(missing(voxdims)){
		if(is.gjdens(inds))
			voxdims=voxdim.gjdens(inds)
		else
			stop("no voxdims supplied and inds has no physical dimension attributes")
	} else if(is.gjdens(voxdims)){
		if(missing(origin))
			origin=matrix(getBoundingBox(voxdims),nrow=2)[1,]
		voxdims=voxdim.gjdens(voxdims)
	}

	if(missing(origin)){
		if(is.gjdens(inds))
			origin=matrix(getBoundingBox(inds),nrow=2)[1,]
		else
			origin=rep(0,length(dims))
	}

	ind2coord.default(inds, dims=dims, voxdims=voxdims, origin=origin, ...)
}

ind2coord.default<-function(inds, dims, voxdims, origin, axperm=NULL){
	# ind2coord find XYZ coords corresponding to 1D indices into a 3D image
	# 
	# Usage: coords = ind2coord(dims, inds, voxdims, [axperm])
	# 
	# Input:
	# inds    - indices into an image array 
	#           (either 1d when dims must be present or a logical array)
	# dims    - dimensions of 3d image array
	# voxdims - vector of 3 voxel dimensions (width, height, depth, dx,dy,dz)
	# axperm  - 3-vector containing reordering of axes, negatives imply flip
	# 
	# coords  - 3xN XYZ triples 
	# 
	# Permutations and Flips:
	# -----------------------
	# axperm = [1 2 3]  do nothing
	# axperm = [-1 2 3] flip the points along the array's first axis
	# axperm = [-2 1 3] flip 1st axis and then swap 1st and 2nd
	# 
	# axperm = [1 -2 3] fixes coords from images loaded by readpic
	# axperm = [2 1 3] fixes coords from images loaded by imread (eg tif)
	# 
	# See also coord2ind, ind2sub
	
	# FIXME - see xyzpos.gjdens for details of handling voxels 
	# amira or imagej style (cell vs node in nrrd terminology)

	if(length(dims) != 3 )
		stop('coords2ind only handles 3d data')
	if(is.matrix(voxdims))
		voxdims=as.numeric(voxdims)
	if(length(voxdims)!=length(dims))
		stop('number of voxel dimensions must match dimensionality of data')

	if(is.array(inds)){
		if(is.logical(inds))
			pixcoords = which(inds, arr.ind=TRUE)
		else if(is.integer(inds) || is.raw(inds))
			pixcoords = which(inds>0, arr.ind=TRUE)
		else stop("cannot handle numeric arrays - make a logical")
	} else if(is.logical(inds)){
		# 1d logical
		pixcoords=arrayInd(which(inds),.dim=dims)
	} else {
		# numbers 
		pixcoords=arrayInd(inds,.dim=dims)
	}
	
	if(nrow(pixcoords)==0) return(NULL)
	
	if(!is.null(axperm)){
		stop("axis permutation is not yet implemented")
		# flip and swap axes if required
		# for i=1:3
		# 	# NB (pix)coords are 0 indexed whereas subscripts are 1 indexed,
		# 	# so we subtract 1 implicitly or explixicitly below
		# 	if(axperm(i)<0)
		# 		# flip axis (NB siz(i)-1-pixcoords would give 1-indexed flip)
		# 		pixcoords(i,:)=siz(i)-pixcoords(i,:);
		# 	else
		# 		pixcoords(i,:)=pixcoords(i,:)-1;
		# 	end
		# end
		# pixcoords=pixcoords(abs(axperm),:);
	}

	# then convert from pixel coords to physical coords
	# transpose to allow multiplication, then back again to give 3 cols
	# note that we must subtract 1 from 1-indexed pixcoords
	rval = if(missing(origin))
		t(t(pixcoords-1)*voxdims)
	else
		t(t(pixcoords-1)*voxdims+origin)
	colnames(rval)=c("X","Y","Z")
	rval
}

coord2ind<-function(coords,imdims,voxdims,origin,aperm,Clamp=FALSE,CheckRanges=!Clamp){
	# finds 1d indices into 3d image array
	# coords  - N x 3 XYZ triples
	# imdims - dimensions of 3d img array (or the array itself)
	# voxdims - vector of 3 voxel dimensions (width, height, depth, dx,dy,dz)
	# aperm   - permutation order for axes
	
	if(is.array(imdims)){
		if(missing(voxdims))
			voxdims=as.numeric(voxdim.gjdens(imdims))
		if(missing(origin))
			origin=getBoundingBox(imdims)[c(1,3,5)]
		imdims=dim(imdims)
	}
	
	if(length(imdims) != 3)
		stop('coord2ind only handles 3d data')

	if(!is.matrix(coords))
		coords=matrix(coords,byrow=TRUE,ncol=length(coords))
	if(!missing(origin))
		coords=t(t(coords)-origin)

	# first convert from physical coords to pixel coords
	# FIXME surely coords are 0 indexed
	pixcoords=t(round(t(coords)/voxdims))+1

	# make sure no points are out of range
	if(Clamp){
		pixcoords[,1]=pmin(imdims[1],pmax(1,pixcoords[,1]))
		pixcoords[,2]=pmin(imdims[2],pmax(1,pixcoords[,2]))
		pixcoords[,3]=pmin(imdims[3],pmax(1,pixcoords[,3]))
	} else if(CheckRanges){
		ranges=apply(pixcoords,2,range)
		if(any(ranges[2,]>imdims) || any(ranges[1,]<1))
			stop("pixcoords out of range")
	}

	# convert to 1d indices
	if (!missing(aperm))
		imdims=imdims[aperm]
	sub2ind(imdims,pixcoords)
}

sub2ind<-function(dims,coords){
	# emulate matlab's sub2ind command

	# convert vector containing 1 coordinate into matrix
	if(!is.matrix(coords))
		coords=matrix(coords,byrow=TRUE,ncol=length(coords))
	if(length(dims)!=ncol(coords)){
		stop("coords must have the same number of columns as dimensions in dims")
	}
	k=cumprod(c(1,dims[-length(dims)]))
	ndx=1
	for(i in 1:length(dims)){
		v=coords[,i]
		if(any(v<1) || any(v>dims[i]))
			stop("index out of range")
		ndx=ndx+(v-1)*k[i]
	}
	ndx
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
