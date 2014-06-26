# Affine.R
# ########
# Functions to generate 3D affine transformation matrices
# for points in homogeneous co-ordinates
###########
# 040219
# 040422 - added functions to compose affine matrices according
# to details supplied by Torsten Rohlfing

#RELEASE
#BEGINCOPYRIGHT
###############
# R Source Code to accompany the manuscript
#
# "Comprehensive Maps of Drosophila Higher Olfactory Centers: 
# Spatially Segregated Fruit and Pheromone Representation"
# Cell (2007), doi:10.1016/j.cell.2007.01.040
# by Gregory S.X.E. Jefferis*, Christopher J. Potter*
# Alexander M. Chan, Elizabeth C. Marin
# Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo
#
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# See flybrain.stanford.edu for further details
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################
#ENDMAINCOPYRIGHT
# Copyright Addendum:
# DecomposeAffineToIGSParams and ComposeAffineFromIGSParams
# are translations to R of code originally copyright Torsten Rohlfing 
# 

matI=matrix(c(1,0,0,0,
              0,1,0,0,
              0,0,1,0,
              0,0,0,1), ncol=4, byrow=TRUE)

AffineRotation<-function(rx=0,ry=0,rz=0,Degrees=F){
	# expects rotations around the x, y and z axes
	# in radians unless Degrees == T
	if(Degrees){
		rx=rx*pi/180
		ry=ry*pi/180
		rz=rz*pi/180	
	}
	
	zrot=as.matrix( rbind(
		c(cos(rz),-sin(rz),0,0),
		c(sin(rz),cos(rz),0,0),
		c(0,0,1,0),
		c(0,0,0,1)  ))
	yrot=as.matrix( rbind(
		c(cos(ry),0,-sin(ry),0),
		c(0,1,0,0),
		c(sin(ry),0,cos(ry),0),
		c(0,0,0,1)  ))
	xrot=as.matrix( rbind(
		c(1,0,0,0),
		c(0,cos(rx),-sin(rx),0),
		c(0,sin(rx),cos(rx),0),
		c(0,0,0,1)  ))
	# return the product of these three transformations
	xrot%*%yrot%*%zrot
}

AffineInvRotation<-function(rx=0,ry=0,rz=0,Degrees=F){
	# expects rotations around the x, y and z axes
	# in radians unless Degrees == T
	if(Degrees){
		rx=rx*pi/180
		ry=ry*pi/180
		rz=rz*pi/180	
	}
	
	zrot=as.matrix( rbind(
		c(cos(rz),-sin(rz),0,0),
		c(sin(rz),cos(rz),0,0),
		c(0,0,1,0),
		c(0,0,0,1)  ))
	yrot=as.matrix( rbind(
		c(cos(ry),0,-sin(ry),0),
		c(0,1,0,0),
		c(sin(ry),0,cos(ry),0),
		c(0,0,0,1)  ))
	xrot=as.matrix( rbind(
		c(1,0,0,0),
		c(0,cos(rx),-sin(rx),0),
		c(0,sin(rx),cos(rx),0),
		c(0,0,0,1)  ))
	# return the product of these three transformations
	zrot%*%yrot%*%xrot
}


AffineScale<-function(sx=1,sy=1,sz=1){
	as.matrix( rbind(
			c(sx,0,0,0),
			c(0,sy,0,0),
			c(0,0,sz,0),
			c(0,0,0,1) ))
}

AffineTranslate<-function(tx=0,ty=0,tz=0,WidthPixel=1,HeightPixel=1,DepthPixel=1){
	as.matrix( rbind(
			c(1,0,0,tx/WidthPixel),
			c(0,1,0,ty/HeightPixel),
			c(0,0,1,tz/DepthPixel),
			c(0,0,0,1) ))
}

AffineTranslateRotateScale<-function(tx=0,ty=0,tz=0,rx=0,ry=0,rz=0,
  sx=1,sy=1,sz=1,Degrees=F,...){
	AffineScale(sx,sy,sz)%*%AffineRotation(rx,ry,rz,Degrees=Degrees)%*%AffineTranslate(tx,ty,tz,...)
}

#' Compose homogeneous affine matrix from CMTK registration parameters
#' 
#' Deprecated
#' @details If the \code{legacy} parameter is not set explicitly, then it will 
#'   be set to \code{TRUE} if params has a version attribute <2.4 or FALSE 
#'   otherwise.
#' @param params 5x3 matrix of CMTK registration parameters or list of length 5.
#' @param legacy Whether to assume that parameters are in the format uses by 
#'   CMTK <=2.4.0 (default FALSE, see details).
#' @param tx,ty,tz Translation along x, y and z axes (default 0)
#' @param rx,ry,rz Rotation about x, y and z axes (in degrees, default 0)
#' @param sx,sy,sz Scale for x, y and z axes (default 1)
#' @param shx, shy, shz Shear for x,y,z axes (default 0)
#' @param cx,cy,cz Centre for rotation
#' @return 4x4 homogeneous affine transformation matrix
#' @details translation and centre components are assumed to be in physical 
#'   coordinates.
#' @export
ComposeAffineFromIGSParams<-function(params=NULL, tx=0, ty=0, tz=0, rx=0, ry=0, 
  rz=0, sx=1, sy=1, sz=1, shx=0, shy=0, shz=0, cx=0, cy=0, cz=0, legacy=NA){

	.Deprecated('nat::cmtkparams2affmat')
	cmtkparams2affmat(params=params, tx=tx, ty=ty, tz=tz, rx=rx, ry=ry, 
		rz=rz, sx=sx, sy=sy, sz=sz, shx=shx, shy=shy, shz=shz, cx=cx, cy=cy, cz=cz, 
		legacy=NA)
}

MakeAffineFromReg<-function(fileName,Inverse=FALSE){
	# Return the fully composed homogeneous affine transformation matrix
	# based on the affine registration file
	# Can return the inverse of the affine if Inverse=TRUE
	reg=ReadIGSRegistration(fileName)
	if(!is.null(reg$spline_warp)){
		AffInstructions=reg$spline_warp$affine_xform
		# a warp transformationm so get the second affine
	} else {
		AffInstructions=reg$affine_xform
	}
	AffInstructions=do.call(rbind,AffInstructions)
	#return(AffInstructions)
	AffMat=ComposeAffineFromIGSParams(AffInstructions)
	if(Inverse) AffMat=solve(AffMat)
	return(AffMat)
}

TransformPoints<-function(Points,AffMat,dim=3){
	# function to Transform a set of points in which each point
	# occupies 1 row (and dim cols)
	# using an affine transformation matrix

	# convert points to homogeneous co-ords by appending 1 to all of them
	# note that it is necessary to transpose to
	# have points as column vectors
	if(is.vector(Points)) Points=matrix(Points,nrow=1)
	MyPoints=t(cbind(Points[,1:dim,drop=F],1))
	# apply transformation
	# transpose back to get a list of points and
	# drop the last column 1 to get dim numbers again for each dim-d point
	rval=t(AffMat%*%MyPoints)[,1:dim,drop=F]
	# add back the other cols if there were more than dim cols in
	# the supplied points (eg width)
	if(ncol(Points)>3) rval=cbind(rval,Points[,-(1:dim)])

	return(rval)
}

CalculateAffineFromLandmarkPairs<-function(landmarks,StartMat,
	SwapTransform=FALSE,Verbose=FALSE,costfn=function(x) sum(abs(x)),method="BFGS",...){

	# Get landmark data
	if(!is.list(landmarks)) d=ReadAmiraLandmarks(landmarks)
	else d=landmarks
	if(!is.list(d)) stop("Unable to read landmarks file")
	if(SwapTransform) d=d[2:1]
	
	# Choose a sensible starting matrix
	if(missing(StartMat)) {
		StartMat=as.matrix( rbind(
			c(1,0,0,0),
			c(0,1,0,0),
			c(0,0,1,0),
			c(0,0,0,1)  ) )
		# set a default translation
		StartMat[1:3,4]=colMeans(d[[1]]-d[[2]])
	}
	if(Verbose) cat("Starting Matrix is:",StartMat,"\n")
	# Use optim
	# need a little cost function
	# use a least squares min?
	
	SquaredDiffUnderAffine<-function(affmatcomps){
		affmat=matrix(affmatcomps,ncol=4,byrow=T)
		affmat=rbind(affmat,c(0,0,0,1))
		diffmat=TransformPoints(d[[1]],affmat)-d[[2]]
		rval=sum(diffmat^2)
		if(Verbose) cat(rval,",")
		rval
	}
	SquaredDiffUnderAffine<-function(affmatcomps){
		affmat=matrix(affmatcomps,ncol=4,byrow=T)
		affmat=rbind(affmat,c(0,0,0,1))
		diffmat=TransformPoints(d[[1]],affmat)-d[[2]]
		rval=costfn(diffmat)
		if(Verbose) cat("\n[",affmat,sep=",","]")
		if(Verbose) cat(rval,",")
		rval
	}
	
	initaffmatcomps=t(StartMat)[1:12]
	optresults=optim(initaffmatcomps,SquaredDiffUnderAffine,method=method,...)
	if(optresults$convergence!=0)
		cat("Warning: optimisation returned convergence error code",optresults$convergence,"\n")
	finalaffmatcomps=optresults$par
	cat("Final score = ",optresults$value,"\n")

	finalaffmat=matrix(finalaffmatcomps,ncol=4,byrow=T)
	finalaffmat=rbind(finalaffmat,c(0,0,0,1))
	finalaffmat
}

CalculateIGSParamsFromLandmarkPairs<-function(landmarks,dofs=c(3,6,9,12),
	StartParams=c(0,0,0,0,0,0,1,1,1,0,0,0,84.39,84.39,43.5),
	SwapTransform=FALSE,Verbose=FALSE,costfn=function(x) sum(abs(x)),method="BFGS",...){
	dofs=dofs[1]
	# Get landmark data
	if(!is.list(landmarks)) d=ReadAmiraLandmarks(landmarks)
	else d=landmarks
	if(!is.list(d)) stop("Unable to read landmarks file")
	if(SwapTransform) d=d[2:1]

	# Check StartParams
	if(is.matrix(StartParams))  StartParams=as.vector(t(StartParams))
	# Use optim
	# need a little cost function
	# use a least squares min?
	
	CostFNUnderTParams<-function(params){
		fullparams=StartParams; fullparams[1:dofs]=params
		affmat=ComposeAffineFromIGSParams(fullparams)
		diffmat=TransformPoints(d[[1]],affmat)-d[[2]]
		rval=costfn(diffmat)
		if(Verbose) cat(rval,",")
		rval
	}

	optresults=optim(StartParams[1:dofs],CostFNUnderTParams,method=method,...)
	finalparams=optresults$par
	finalfullparams=StartParams; finalfullparams[1:dofs]=finalparams
	if(optresults$convergence!=0)
		cat("Warning: optimisation returned convergence error code",optresults$convergence,"\n")
	cat("Final score = ",optresults$value,"\n")
	
	m=matrix(finalfullparams,ncol=3,byrow=TRUE)
	attr(m,"optresults")=optresults
	m
}

#' Decompose homogeneous affine matrix to CMTK registration parameters
#'
#' @details deprecated
#' @param matrix 4x4 homogeneous affine matrix
#' @param centre Rotation centre
#' @return 5x3 matrix of CMTK registration parameters
#' @export
#' @seealso \code{\link{ComposeAffineFromIGSParams}}
DecomposeAffineToIGSParams<-function(matrix,centre=c(0,0,0)){
	.Deprecated('nat::affmat2cmtkparams')
	nat::affmat2cmtkparams(matrix,centre=centre)
}
