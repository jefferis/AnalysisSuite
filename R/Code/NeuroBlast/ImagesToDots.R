# Functions to process images of neurons and generate more compact
# dot-based representations
# Based on Nick Masse's matlab code (esp extract_properties.m)

require(RANN)

is.dotprops<-function(dp) inherits(dp,"dotprops")

as.dotprops<-function(dp){
	if(is.null(dp)) return (NULL)
	if(!is.dotprops(dp)) class(dp)=c("dotprops",class(dp))
	if(is.null(colnames(dp$points))) colnames(dp$points) <-c("X","Y","Z") 
	dp
}

plot3d.dotprops<-function(dp, scalevecs=1.0, alpharange=NULL,
	PlotPoints=FALSE, PlotVectors=TRUE, UseAlpha=FALSE,...){
	# rgl's generic plot3d will dispatch on this
	if (!is.null(alpharange))
		dp=subset(dp,dp$alpha<=alpharange[2] & dp$alpha>=alpharange[1])
	rlist=list()
	if(PlotPoints)
		rlist$points=points3d(dp$points,...)
	if(PlotVectors){
		halfvect=dp$vect/2*scalevecs
		if(UseAlpha) halfvect=halfvect*dp$alpha
		starts=dp$points-halfvect
		stops=dp$points+halfvect
		interleaved=matrix(t(cbind(starts,stops)),ncol=3,byrow=T)
		rlist$segments=segments3d(interleaved,...)
	}
	invisible(rlist)
}

xyzmatrix<-function(x,y=NULL,z=NULL,Transpose=FALSE) {
	# quick function that gives a generic way to extract coords from 
	# classes that we care about and returns a matrix
	# nb unlike xyz.coords this returns a matrix (not a list)
	x=if(is.neuron(x)) x$d[,c("X","Y","Z")]
	else if(is.dotprops(x)) x$points
	else if(!is.null(z)){
		cbind(x,y,z)
	} else x
	mx=data.matrix(x)
	if(Transpose) t(mx) else mx
}

#' Assign xyz elements of neuron or dotprops object
#'
#' Can also handle matrix like objects with cols called X,Y,Z
#' @param n dotprops/neuron/data.frame/named matrix
#' @param value Nx3 matrix specifying xyz coords
#' @return Original object with modified coords
#' @export
#' @seealso \code{\link{xyzmatrix}}
#' @examples
#' n=MyNeurons[[1]]
#' xyzmatrix(n)<-xyzmatrix(n)
#' stopifnot(isTRUE(
#'   all.equal(xyzmatrix(n),xyzmatrix(MyNeurons[[1]]))
#' ))
`xyzmatrix<-`<-function(n,value){
  if(is.neuron(n)) n$d[,c("X","Y","Z")]=value
  else if(is.dotprops(n)) n$points[,c("X","Y","Z")]=value
  else if(c("X","Y","Z") %in% colnames(n)) n[,c("X","Y","Z")]=value
  else stop("Not a neuron or dotprops object or a matrix-like object with XYZ volnames")
  n
}

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

#' Subset points in dotprops object that match given conditions
#'
#' subset defines either logical or numeric indices, in which case these are
#' simply applied to the matrices that define the points, vect etc
#' OR a function (which is called with the 3d points array and returns T/F vector)
#' OR another dotprops in which case prune.dotprops is called
#' @param dp A dotprops object
#' @param subset indices, function or another dotprops (see Details)
#' @param ... Additional parameters passed to prune.dotprops (see Details)
#' @return subsetted dotprops object
#' @export
#' @seealso \code{\link{prune.dotprops}}
#' @examples
#' \dontrun{
#' s3d=select3d()
#' dp1=subset(dp,s3d(points))
#' # special case of previous version
#' dp2=subset(dp,s3d)
#' stopifnot(all.equal(dp1,dp2))
#' dp2=subset(dp,alpha>0.5 & s3d(pointd))
#' dp3=subset(dp,1:10)
#' }
subset.dotprops<-function(dp,subset,...){
	e <- substitute(subset)
	r <- eval(e, dp, parent.frame())
	if (!is.logical(r) && !is.numeric(r)) {
		# a function that tells us whether a point is in or out
		if(is.function(r)) r=subset(dp$points)
		else if(is.dotprops(r)) return(prune.dotprops(dp,subset,...))
		else stop("Cannot evaluate subset")
	}
	if(is.logical(r)) r <- r & !is.na(r)
	else if(!is.numeric(r)) stop("Subset must evaluate to a logical or numeric index")
	
	dp$points=dp$points[r,,drop=F]
	dp$alpha=dp$alpha[r]
	dp$vect=dp$vect[r,,drop=F]
	if(!is.null(dp$labels)) dp$labels=dp$labels[r]
	dp
}

length.dotprops<-function(dp) nrow(dp$points)

# redefining length upsets str
str.dotprops<-function(dp,...) {class(dp)<-"list";str(dp,...)}

#' Carry out in memory digest of a dotprops object
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
  # remove mtime and file attributes
  atts=attributes(dp)
  mostattributes(dp)<-atts[setdiff(names(atts),c("mtime",'file'))]
  digest(dp,...)
}

# Set up a private function that directly accesses internal LAPACK routine of eigen
if(R.version$major<3){
  .internal_lars<-function(x) .Call("La_rs", x, only.values=FALSE, PACKAGE = "base")
} else {
  .internal_lars<-function(x) .Internal(La_rs(x, only.values=FALSE))
}

DotProperties<-function(points,k=20,UseLabels=TRUE){
	# store labels from SWC format data if this is a neuron
	Labels=if(UseLabels && is.neuron(points)) points$d$Label else NULL
	points=xyzmatrix(points)
	npoints=nrow(points)
	if(npoints<k) stop("Too few points to calculate properties")
	if(ncol(points)!=3) stop("points must be a N x 3 matrix")
	
	alpha=rep(0,npoints)
	vect=matrix(0,ncol=3,nrow=npoints)

	nns=nn2(points,points,k=k)
	# transpose points to 3xN because 
	# R arithemtic of matric / vector operates column-wise
	pointst=t(points)
	for(i in 1:npoints){
		indNN=nns$nn.idx[i,]
		
		pt=pointst[,indNN]
		cpt=pt-rowMeans(pt)
		
		inertia=matrix(0,ncol=3,nrow=3)
		diag(inertia)=rowSums(cpt^2)
		inertia[1,2]<-inertia[2,1]<-sum(cpt[1,]*cpt[2,])
		inertia[1,3]<-inertia[3,1]<-sum(cpt[1,]*cpt[3,])
		inertia[2,3]<-inertia[3,2]<-sum(cpt[2,]*cpt[3,])
		
		# call internal LAPACK routine of eigen directly
		z<-.internal_lars(inertia)
		ord <- rev(seq_along(z$values))
		v1d1=list(values = z$values[ord], vectors = z$vectors[,ord, drop = FALSE])

		alpha[i]=(v1d1$values[1]-v1d1$values[2])/sum(v1d1$values)
		vect[i,]=v1d1$vectors[,1]
	}
	rlist=list(points=points,alpha=alpha,vect=vect)
	rlist$labels=Labels
	attr(rlist,'k')=k
	return(as.dotprops(rlist))
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
	} else if(ext=="swc") {
		l$points=ReadSWCFile(f)[,c("X","Y","Z")]
	} else if(ext=='csv') {
		l$points=read.csv(f,header=FALSE)
		colnames(l$points)=c("X","Y","Z")
	} else 
		stop("Cannot extract dots from file type: ",ext)
	if(!is.null(xformfun))
		l$points=xformfun(l$points)
	l=DotProperties(l$points,...)
	attr(l,'file')=f
	fi=file.info(f)
	attr(l,'mtime')=fi$mtime
	attr(l,'size')=fi$size
	l
}

#' Transform dot property object using specified registration
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
	if(is.function(reg)){
	  # we've been given a function - apply this to points
	  pointst=reg(dp$points,...)
	} else {
	  # we've been given a CMTK registration file
	  pointst=transformedPoints(xyzs=dp$points,warpfile=reg,transforms='warp',...)[['warp']]
	}
	
	naPoints=is.na(pointst[,1])
	if(any(naPoints)){
		if(na.action=='warn')
			warning("Dropping ",sum(naPoints),' points')
		else if (na.action=='error')
			stop("Error: Failed to transform ",sum(naPoints),' points')
		
		pointst=pointst[!naPoints,]
	}
	dpn=DotProperties(pointst,k)
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
