# Functions to process images of neurons and generate more compact
# dot-based representations
# Based on Nick Masse's matlab code (esp extract_properties.m)

require(RANN)

DotProperties<-function(points,k=20){
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
		z<-.Call("La_rs", inertia, only.values=FALSE, PACKAGE = "base")
		ord <- rev(seq_along(z$values))
		v1d1=list(values = z$values[ord], vectors = z$vectors[,ord, drop = FALSE])

		alpha[i]=(v1d1$values[1]-v1d1$values[2])/sum(v1d1$values)
		vect[i,]=v1d1$vectors[,1]
	}
	return(list(points=points,alpha=alpha,vect=vect))
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
	
	rval = if(missing(origin))
		t(t(pixcoords)*voxdims)
	else
		t(t(pixcoords)*voxdims+origin)
	colnames(rval)=c("X","Y","Z")
	rval
}

# stop()
# 
# points=ReadAmiramesh("/GD/projects/Nick/FruCloneClustering/data/SAKW13-1_gj_dimension_reduced.am")
