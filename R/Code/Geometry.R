# Geometry.R
# 2005-02-23
# some routines to allow geometric calculations that I will
# need for analysis of tracing data after transformation according
# to warps calculated from Torsten's software.

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

# source(file.path(CodeDir,"Geometry.R"))
FindPlaneFrom3Points<-function(ps){
	ps=data.frame(ps)
	colnames(ps)=c("X","Y","Z")
	# given three points, in which each point is a row 
	# and three cols represent X, Y and Z coords,
	# find the equation of the plane containing them
	# in the form
	# Ax + By + Cz + D = 0
	# copied straight from:
	# http://astronomy.swin.edu.au/~pbourke/geometry/planeeq/
	
	A = ps$Y[1] * (ps$Z[2] - ps$Z[3]) + ps$Y[2] * (ps$Z[3] - ps$Z[1]) + ps$Y[3] * (ps$Z[1] - ps$Z[2]) 
	B = ps$Z[1] * (ps$X[2] - ps$X[3]) + ps$Z[2] * (ps$X[3] - ps$X[1]) + ps$Z[3] * (ps$X[1] - ps$X[2]) 
	C = ps$X[1] * (ps$Y[2] - ps$Y[3]) + ps$X[2] * (ps$Y[3] - ps$Y[1]) + ps$X[3] * (ps$Y[1] - ps$Y[2]) 
	D = -1*(ps$X[1] * (ps$Y[2] * ps$Z[3] - ps$Y[3] * ps$Z[2]) + ps$X[2] * (ps$Y[3] * ps$Z[1] - ps$Y[1] * ps$Z[3]) + ps$X[3] * (ps$Y[1] * ps$Z[2] - ps$Y[2] * ps$Z[1]))
	
	return(c(A,B,C,D))
}

FindPlaneFromPointAndNormal<-function(p,n){
	# p is point, n is normal vector
	if(is.null(dim(p))) p=matrix(p,ncol=3)
	if(is.null(dim(n))) n=matrix(n,ncol=3)
	k=apply(n*p,1,sum)
	return(cbind(n,-k))
}



SideOfPlane<-function(ps,plane){
	if(is.null(dim(ps))) ps=matrix(ps,ncol=3)
	# The sign of s = Ax + By + Cz + D 
	# determines which side the point (x,y,z) lies with respect to the plane. 
	#If s > 0 then the point lies on the same side as the normal (A,B,C). 
	#If s < 0 then it lies on the opposite side, 
	#if s = 0 then the point (x,y,z) lies on the plane.
	plane[1]*ps[,1]+plane[2]*ps[,2]+plane[3]*ps[,3]+plane[4]
}

IntersectionLineAndPlane<-function(Line,Plane){
	# the line is defined by two points in a 2X3 where cols are X,Y,Z
	# and rows are the different points
	# the plane is defined the coefficients A,B,C,D
	# for  Ax + By + Cz + D = 0
	
	u=(Plane[1:3]*Line[1,] +Plane[4])/(Plane[1:3]*(Line[1,]-Line[2,]))
	return(Line[1,]+u*(Line[2,]-Line[1,]))
}

IntersectionLineSegmentAndPlane<-function(LineSegs,Plane){
	# the line seg is defined by two points in a 2x3nN where cols are X,Y,Z
	# and rows are the different points and 3rd dim is different line segs
	# the plane is defined the coefficients A,B,C,D
	# for  Ax + By + Cz + D = 0
	# This function will return NA if the plane does not lie between 
	# the 2 points
	
	if(identical(dim(LineSegs),as.integer(c(2,3)))){
		# only one point in 2X3 
		u=(sum(Plane[1:3]*LineSegs[1,]) +Plane[4])/(sum(Plane[1:3]*(LineSegs[1,]-LineSegs[2,])))
		#cat("u =",u,"\n")
		if(u>1 || u<0) return(c(NA,NA,NA))
		return(LineSegs[1,]+u*(LineSegs[2,]-LineSegs[1,]))
	} else {
		# 2x3xN
		# output will have cols for different line segs, rows for each of
		# three points
		u=(apply(Plane[1:3]*LineSegs[1,,],2,sum) +Plane[4])/apply(Plane[1:3]*(LineSegs[1,,]-LineSegs[2,,]),2,sum)
		#rval=LineSegs[1,,]+u*(LineSegs[2,,]-LineSegs[1,,])
		# nb t so that rows are points and cols are x,y,z
		rval=t(LineSegs[1,,])+apply(LineSegs[2,,]-LineSegs[1,,],1,"*",u)
		#cat("u =",u,"\n")
		rval[u>1 | u<0,]=NA
		rval
	}	
}
	
RegularTetrahedron<-function(cent=rep(0,3),a=1){
	# Returns a matrix giving the vertices of a regular tetrahedron
	# centred on c 
	coords=rbind(c(0,0,0),c(1,0,0),c(1/2,sqrt(3)/2,0),c(1/2,sqrt(3)/6,sqrt(3)/2))
	scale(coords,cent=-cent,scale=rep(1/a,3))
}

TetrahedronVolume<-function(v){
	# Expects a matrix v with 4 rows and 3 columns
	# Each row gives the coords of one of the vertices
	
	# nb regular tetrahedron has V=1/3*A*h
	# = 1/3*sqrt(3)/2*a * a/2 * a*sqrt(3)/2
	# where a is side length.
	
	# this regular tetrahedron has coords
	# c=rbind(c(0,0,0),c(1,0,0),c(1/2,sqrt(3)/2,0),c(1/2,sqrt(3)/6,sqrt(3)/2))
	if(isTRUE(all.equal(dim(c),c(4,3)))){
		abs(1/6*det(cbind(v,rep(1,4))))
	} else {
		return(NULL)
	}
}

	
DefineLHEntryPlane<-function(){
	# There should also be a little amira script to define this plane
	# v
	xorig=77
	dx=45
	dz=80
	g=dz/dx
	PointsOnLHEntryPlane=rbind(c(xorig,0,0) ,c(xorig+dx,0,dz),c(xorig+dx/2,50,dz/2))
	FindPlaneFrom3Points(PointsOnLHEntryPlane)
}

LHEntryPlane=DefineLHEntryPlane()
#In amira 
#u = 1...100
#v = 1...50
#x=v+75
#y=u*1.5
#z=v*110/70
MBEntryPlane=FindPlaneFromPointAndNormal(c(52,0,61),c(70,0,27))
#In amira 
#u = 1...100
#v = 1...50
#x=v+23
#y=u*1.5
#z=v*(-80/26)+150

IntersectNeuronWithPlane<-function(ANeuron,APlane){
	# iterate over segments
	# for each segment, find closest point
	# then find neighbour on opposite side
	# compute intersection
	
	DistancesToPlane=list()
	for(i in 1:ANeuron$NumSegs){
		DistancesToPlane=c(DistancesToPlane,list(cbind(ANeuron$d[ANeuron$SegList[[i]],"PointNo"],SideOfPlane(ANeuron$d[ANeuron$SegList[[i]],c("X","Y","Z")],APlane))))
	}
	IntersectionPoints=matrix(NA,nrow=0,ncol=3)
	for(i in 1:ANeuron$NumSegs){
		ZeroPoints=DistancesToPlane[[i]][,2]==0
		if(any(ZeroPoints)) IntersectionPoints=rbind(IntersectionPoints,
			ANeuron$d[DistancesToPlane[[i]][which(ZeroPoints)[1],1],c("X","Y","Z")])
		else { 
			if(any(DistancesToPlane[[i]][,2]>0) && any(DistancesToPlane[[i]][,2]<0)){
				# there are both pos and neg distances.
			}
		}
	}
	
	
	return(DistancesToPlane)
}

IntersectNeuronWithPlane<-function(ANeuron,APlane){
	# consider every line in neuron (ie each point and its parent)
	# make a 3d matrix R,C,i Rows are point 1 and 2, cols are XYZ and i are
	# multiple lines

	if(nrow(ANeuron$d)<2) stop("Need at least 2 points to define a line in Neuron")
	d=ANeuron$d
	points=unique(unlist(ANeuron$SegList))
	
	PointMatrix=array(0,dim=c(2,3,length(points)-1))
	# set up the ends of the lines
	PointMatrix[2,,]=t(as.matrix(d[points[-1],c("X","Y","Z")]))
	# set up the starts of the lines
	PointMatrix[1,,]=t(as.matrix(d[d$Parent[points[-1]],c("X","Y","Z")]))
	rval=IntersectionLineSegmentAndPlane(PointMatrix,APlane)
	# return any non NA rows
	rval[!is.na(rval[,1])]
}

FindPlaneOrigin<-function(Plane){
	# Not sure this is possible!
}
CrossProduct<-function(a,b){
	c(a[2]*b[3]-a[3]*b[2],
		a[3]*b[1]-a[1]*b[3],
		a[1]*b[2]-a[2]*b[1])
}
DotProduct<-function(a,b){
	sum(a*b)
}


FindXYPosOnPlane<-function(Points,Plane,Origin=FindPlaneOrigin(Plane)){
	# find the XY coords of a set of points on the plane
	# using the definition of the origin in the plane
	# 
	# UNFINISHED!
	u=c(0,1,0) # y vector
	n=Plane[1:3] # normal vector of plane
	v=CrossProduct(u,n)
	#normalise
	v=v/sqrt(sum(v^2))
	
	cbind(apply(Points,1,DotProduct,u),apply(Points,1,DotProduct,v))
}

Make3DBins<-function(d,n=10,dx=signif(diff(range(z$Y))/n,2),dy=dx,dz=dx){
	#Takes a set of 3D data points and finds a set of 3D bins that
	#encompass the range of points
	# returns a 3D array
}
