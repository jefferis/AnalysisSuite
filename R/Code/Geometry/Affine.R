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
# DecomposeAffineToIGSParams and ComposeAffineFromIGSParams.named
# are translations to R of code originally copyright Torsten Rohlfing 

# 
# source(file.path(CodeDir,"Affine.R"))
# 
# Identity matrix:
# [1, 0, 0, 0]
# [0, 1, 0, 0]
# [0, 0, 1, 0]
# [0, 0, 0, 1]

matI=as.matrix( rbind(
		c(1,0,0,0),
		c(0,1,0,0),
		c(0,0,1,0),
		c(0,0,0,1)  ) )


# Z rotation:
# [cos t, 0, -sin t, 0]
# [sin t, 0, cos t, 0]
# [0, 0, 1, 0]
# [0, 0, 0, 1]

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

AffineShear<-function(hxy,hxz,hyz){
}

AffineTranslateRotateScale<-function(tx=0,ty=0,tz=0,rx=0,ry=0,rz=0,
  sx=1,sy=1,sz=1,Degrees=F,...){
	AffineScale(sx,sy,sz)%*%AffineRotation(rx,ry,rz,Degrees=Degrees)%*%AffineTranslate(tx,ty,tz,...)
}

ComposeAffineFromIGSParams<-function(params, legacy=FALSE){
	# Expects a 5 x 3 matrix of paramaters
	# in the same order as the affine transformation parameters
	# in the affine registration files
	
	if(is.list(params)) params=unlist(params)
	if(is.vector(params)) params=matrix(params,ncol=3,byrow=T)
	if(identical(dim(params), as.integer(c(3,5)))) params = t(params) else {
		if(!identical(dim(params),as.integer(c(5,3)))) stop("unable to parse params")
	}
	
	return(ComposeAffineFromIGSParams.named( tx=params[1,1],ty=params[1,2],tz=params[1,3],
	rx=params[2,1],ry=params[2,2],rz=params[2,3],
	sx=params[3,1],sy=params[3,2],sz=params[3,3],
	shx=params[4,1],shy=params[4,2],shz=params[4,3],
	cx=params[5,1],cy=params[5,2],cz=params[5,3], legacy=legacy))
}

ComposeAffineFromIGSParams.named<-function(tx=0,ty=0,tz=0,rx=0,ry=0,rz=0,sx=1,sy=1,
  sz=1,shx=0,shy=0,shz=0,cx=0,cy=0,cz=0, legacy=FALSE){
	# Compose an affine transformation matrix in an identical fashion to 
	# IGS's affine matrix
	#- params[0..2] tx,ty,tz
	# params[3..5] rx,ry,rz (in degrees)
	# params[6..8] sx,sy,sz
	# params[9..11] shx,shy,shz
	# params[12..14] cx,cy,cz
	DegToRad=function(theta) theta/360*2*pi
	
	alpha = DegToRad(rx)
	theta = DegToRad(ry)
	  phi = DegToRad(rz)
	 
	cos0 = cos(alpha)
	sin0 = sin(alpha)
	cos1 = cos(theta)
	sin1 = sin(theta)
	cos2 = cos(  phi)
	sin2 = sin(  phi)

	sin0xsin1 = sin0 * sin1
	cos0xsin1 = cos0 * sin1

	rval=matrix(0,4,4)
	diag(rval)<-1
	# nb in R matrices are indexed m[row,col]
	# whereas in C looks like T indexed them
	# m[col-1][row-1]
	# Regexps to transform these forms:
	# \[\d\]\[\d\]
	# (\2+1,\1+1)
	# 
	rval[0+1,0+1] =  cos1*cos2
	rval[1+1,0+1] = -cos1*sin2
	rval[2+1,0+1] = -sin1
	rval[0+1,1+1] =  (sin0xsin1*cos2 + cos0*sin2)
	rval[1+1,1+1] = (-sin0xsin1*sin2 + cos0*cos2)
	rval[2+1,1+1] =  sin0*cos1
	rval[0+1,2+1] =  (cos0xsin1*cos2 - sin0*sin2)
	rval[1+1,2+1] = (-cos0xsin1*sin2 - sin0*cos2)
	rval[2+1,2+1] =  cos0*cos1

  if(legacy){
    rval[1:3,1]=rval[1:3,1]*sx
    rval[1:3,2]=rval[1:3,2]*sy
    rval[1:3,3]=rval[1:3,3]*sz
    shears=c(shx,shy,shz)
    if(TRUE){
      # generate shears in broken CMTK <2.4.0 style
      for (i in 3:1 ) {
        shear=matrix(0,4,4)
        diag(shear)<-1
        # i/2 {0,0,1} for i={0,1,2}
        # (i/2)+(i%2)+1 {1,2,2} for i={0,1,2}
        # shear[i/2][(i/2)+(i%2)+1] = dofs[9+i];
        shear[c(2,3,3)[i],c(1,1,2)[i]]=shears[i]
        rval = shear%*%rval
      }
    } else {
      # Generate shears in broken early IGS form
      for (i in 3:1 ) {
        for (j in 1:3) {
          rval[j,i] =rval[j,i] +shears[i] * rval[j,i%%3+1]
        }
      }
    }
  } else {
    # generate scales and shears according to CMTK >=v.2.4.0 / svn r5050
    scaleShear=matrix(0,4,4)
    diag(scaleShear)=c(sx,sy,sz,1)
    scaleShear[0+1,1+1]=shx
    scaleShear[0+1,2+1]=shy
    scaleShear[1+1,2+1]=shz
    
    # NB matrix multiplication must be in opposite order from C original
    rval = rval%*%scaleShear
  }
	
	# transform rotation center
	cM = c(  cx*rval[0+1,0+1] + cy*rval[0+1,1+1] + cz*rval[0+1,2+1],
	  cx*rval[1+1,0+1] + cy*rval[1+1,1+1] + cz*rval[1+1,2+1],
	  cx*rval[2+1,0+1] + cy*rval[2+1,1+1] + cz*rval[2+1,2+1]  )

	# set translations
	rval[0+1,3+1] = tx - cM[1] + cx
	rval[1+1,3+1] = ty - cM[2] + cy
	rval[2+1,3+1] = tz - cM[3] + cz
	
	return(rval)
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

	
AffineTranslateRotateScale(
-6.400000095, 4.300000064, 12.00000005,
-3.800542432, 1.256377663, -0.9860968818,
1.054021854, 0.9631939015, 1.039477509 ,
Degrees=T,Width=.33,Height=.33,Depth=1)

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


FindTorstenFromAffine<-function(affmat,centre=c(84.39,84.39,43.5)){
# INCOMPLETE
	TorstenGuess=rbind(affmat[1:3,4],rep(0,3),rep(1,3),rep(0,3),centre)
	TorstenGuessComps=t(TorstenGuess)[1:12]
	TorstenGuess=optim()
}

#' Decompose homogeneous affine matrix to CMTK registration parameters
#'
#' @param matrix 4x4 homogeneous affine matrix
#' @param centre Rotation centre
#' @return 5x3 matrix of CMTK registration parameters
#' @export
#' @seealso \code{\link{ComposeAffineFromIGSParams}}
DecomposeAffineToIGSParams<-function(matrix,centre=c(0,0,0)){
# C matrices indexed [C][R] (vs [R,C] in R)
# so it seems easiest to transpose for use in R
	matrix=t(matrix)
	# params will contain Torsten's transformation parameters
	params=numeric(15)
	
	# translation entries
	params[1:3] = matrix[4,1:3]

	cM=c(
		sum(centre[1:3]*matrix[1:3,1]),
		sum(centre[1:3]*matrix[1:3,2]),
		sum(centre[1:3]*matrix[1:3,3]) )

	params[1:3] = params[1:3] + cM[1:3] - centre[1:3]
	params[13:15]=centre

	# QR decomposition
	matrix2d=t(matrix[1:3,1:3])
	qr.res=qr(matrix2d)
	Q=qr.Q(qr.res)
	R=qr.R(qr.res)
	R[lower.tri(R)]=0
	
	for (k in 1:3) {
		# if scale is negative, make positive and correct Q and R accordingly (we will figure out later if the overall transformation is a true rotation or has a negative determinant)
		if ( R[k,k] < 0 ) {
			R[k,1:3] = -R[k,1:3]
			Q[1:3,k] = -Q[1:3,k]
		}
		
		# scale
		params[6 + k] = R[k,k]
		
		# report error on singular matrices.
		if ( params[6+k]	< .Machine$double.eps ) stop("singular matrix")
		
		# shear: i,j index the upper triangle of aMat, which is R from QR
		i = (c(0, 0, 1)+1)[k]  # i.e. i := { 0, 0, 1 }
		j = (c(1, 2, 2)+1)[k]  # i.e. j := { 1, 2, 2 }
		params[9+k] = R[i,j]
	}
	
	# =========================================================================
	# 
	# THE FOLLOWING CODE WAS ADOPTED AND MODIFIED FROM VTK, The Visualization
	# Toolkit.
	# 
	#		Program:	 Visualization Toolkit
	#		Language:	 C++
	#		Thanks:		 Thanks to David G. Gobbi who developed this class.
	# 
	# Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
	# All rights reserved.
	# 
	# Redistribution and use in source and binary forms, with or without
	# modification, are permitted provided that the following conditions are met:
	# 
	#	 * Redistributions of source code must retain the above copyright notice,
	#		 this list of conditions and the following disclaimer.
	# 
	#	 * Redistributions in binary form must reproduce the above copyright notice,
	#		 this list of conditions and the following disclaimer in the documentation
	#		 and/or other materials provided with the distribution.
	# 
	#	 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
	#		 of any contributors may be used to endorse or promote products derived
	#		 from this software without specific prior written permission.
	# 
	#	 * Modified source versions must be plainly marked as such, and must not be
	#		 misrepresented as being the original software.
	# 
	# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
	# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
	# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
	# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	# 
	# =========================================================================
	
	if (det(matrix[1:3,1:3])<0){
		# negative determinant, not a true rotation
		# negate x scale
		params[6+1] = -params[6+1];
		# also negative shears related to x
		params[9:10+1] = -params[9:10+1];
	}
	
	VTK.AXIS.EPSILON=1E-8

	# FIXED FROM HERE
	# rotation
	# first rotate about y axis
	x2 = Q[2,1] / params[7]
	y2 = Q[3,1] / params[7]
	z2 = Q[1,1] / params[7]
		
	x3 = Q[2,3] / params[9]
	y3 = Q[3,3] / params[9]
	z3 = Q[1,3] / params[9]
		
	dot = x2 * x2 + z2 * z2
	d1 = sqrt (dot)
		
	if (d1 < VTK.AXIS.EPSILON) {
		cosTheta = 1.0
		sinTheta = 0.0
	} else {
		cosTheta = z2 / d1
		sinTheta = x2 / d1
	}
	rad2deg=function(theta) theta/(2*pi)*360

	params[6] = rad2deg(-atan2 (sinTheta, cosTheta)) # theta
		
		# now rotate about x axis
	dot = x2 * x2 + y2 * y2 + z2 * z2
	d = sqrt (dot)
		
	if (d < VTK.AXIS.EPSILON) {
		sinPhi = 0.0
		cosPhi = 1.0
	} else 
		if (d1 < VTK.AXIS.EPSILON) {
		sinPhi = y2 / d
		cosPhi = z2 / d
		} else {
		sinPhi = y2 / d
		cosPhi = ( x2 * x2 + z2 * z2) / (d1 * d)
		}
		
	params[5] = rad2deg(-atan2 (sinPhi, cosPhi)) # phi 
		
	# finally, rotate about z
	x3p = x3 * cosTheta - z3 * sinTheta
	y3p = - sinPhi * sinTheta * x3 + cosPhi * y3 - sinPhi * cosTheta * z3
	dot = x3p * x3p + y3p * y3p
	d2 = sqrt (dot)
	if (d2 < VTK.AXIS.EPSILON) {
		cosAlpha = 1.0
		sinAlpha = 0.0
	} else {
		cosAlpha = y3p / d2
		sinAlpha = x3p / d2
	}
		
	params[4] = rad2deg(-atan2 (sinAlpha, cosAlpha)) # alpha

	# /** END OF ADOPTED VTK CODE **/
	return (matrix(params,ncol=3,byrow=TRUE))
}
