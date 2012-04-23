# source(file.path(CodeDir,"MyDensity.R"))
# Code that can be used for smoothing in place of matlab kde or sm
# 2006-02-23 Presently using this for Mushroom body data

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

bin3d<-function(d,nbins=rep(50,3),bounds){
	#function to bin a set of data in 3d
#	d=subset(d,X>=bounds[1] & X<=bounds[2] & Y>=bounds[3] & Y<=bounds[4] & Z>=bounds[5] & Z<=bounds[6])
	d=d[(d[,1]>=bounds[1] & d[,1]<=bounds[2] & d[,2]>=bounds[3] & d[,2]<=bounds[4] & d[,3]>=bounds[5] & d[,3]<=bounds[6]),]
	xbins=seq(bounds[1],bounds[2],len=nbins[1]+1)
	ybins=seq(bounds[3],bounds[4],len=nbins[2]+1)
	zbins=seq(bounds[5],bounds[6],len=nbins[3]+1)
	xidxs=cut(d[,1],xbins,right=F)
	yidxs=cut(d[,2],ybins,right=F)
	zidxs=cut(d[,3],zbins,right=F)
	table(xidxs,yidxs,zidxs)
}

smooth2dbysm=function(display='none',...){
		# calls sm.density and returns an object in something close to my own
		# density format.
		r=sm.density(display=display,...)
		rval=r$estimate
		attr(rval,'x')=r$eval.points[,1]
		attr(rval,'y')=r$eval.points[,2]
		rval
}

smooth3d=function(d,sigma2=diag(2.5,3),
	bounds=apply(d[,1:3],2,range),
	nbins=rep(50,3),
	ksize,ksizesds=2,Verbose=TRUE,...){
	
	require(AnalyzeFMRI)
	
	# note that sigma2 is the co-variance matrix so is
	# the sd^2 for a standard gaussian
	
	# ksizesds gives the number of standard deviations to use when
	# calculating the kernel size
	
	bounds=as.vector(bounds)
	# first bin data
	if(length(nbins)==1) nbins=rep(nbins,3)
	if(length(sigma2)==1) sigma2=diag(sigma2,3)
	b=bin3d(d,nbins,bounds)

	voxdim=as.vector(diff(matrix(bounds,nrow=2))/nbins)
	if(Verbose) print(voxdim)
	
	# Calculate minimum size of kernel - note fudge factor of 0.6
	# OLD STYLE - too small for small kernels
	# too large for large ones
	kthresh=ceiling(max(diag(sigma2)/voxdim/.6))

	# Calculate minimum size of kernel
	# defaults to +/- 2 SDs ie 4SDs
	kthresh=ceiling(max( sqrt(diag(sigma2))/voxdim ))*ksizesds*2
	# nb must be odd
	if (kthresh%%2==0) kthresh=kthresh+1

	if(missing(ksize)) {
			ksize=kthresh
	} else {
			if(ksize<kthresh) warning("Kernel is below recommended size")
	}	
	if(Verbose) cat("kernel size is:",ksize,"\n")
	rval=GaussSmoothArray(b,sigma=sigma2,voxdim=voxdim,ksize=ksize,...)
	
	BoundingBox=as.vector(rbind(voxdim/2,-voxdim/2))+bounds
	
	attr(rval,'x')=seq(BoundingBox[1],BoundingBox[2],length=nbins[1])
	attr(rval,'y')=seq(BoundingBox[3],BoundingBox[4],length=nbins[2])
	attr(rval,'z')=seq(BoundingBox[5],BoundingBox[6],length=nbins[3])
	attr(rval,'BoundingBox')=BoundingBox
	attr(rval,'bounds')=bounds
	attr(rval,'kernel')=list(voxdim=voxdim,ksize=ksize,sigma2=sigma2)	
	rval
}


kde3d<-function (x, y=NULL, z=NULL,f, h, n = 25, lims = c(range(x), range(y),range(z))) 
{
		if(is.null(y)){
				y=x[,2]
				z=x[,3]
				x=x[,1]
		}
		nx <- length(x); 		
		if (length(y) != nx || length(y)!=length(z)) 
				stop("data vectors must be the same length")

		if (!missing(f) && length(f) != nx )
				stop("data vectors and weighting factors (f) must be the same length")

		if(length(h)==1) h=rep(h,3)
				
		# the grid over which the densities will be evaluated
		# note this means that we are not really using the mid-point of each
		# bin (rather the left most val for the left most bin and the right
		# most for the right most bin)

		gx <- seq(lims[1], lims[2], length = n)
		gy <- seq(lims[3], lims[4], length = n)
		gz <- seq(lims[5], lims[6], length = n)

		# 
		# for each point convert the position of the bin to so many standard
		# deviations away from the point
		ax <- outer(gx, x, "-")/h[1]
		ay <- outer(gy, y, "-")/h[2]
		az <- outer(gz, z, "-")/h[3]
# 		z <- matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n, 
# 				nx))/(nx * h[1] * h[2])
# 		return(list(x = gx, y = gy, z = z))

		prodx=matrix(dnorm(ax), nrow=n)
		prody=matrix(dnorm(ay), nrow=n)
		prodz=matrix(dnorm(az), nrow=n)
		
		# f (where present) is the weighting for each point
		if(!missing(f)) prodx=scale(prodx,cent=F,scale=f)
		
		prodxyz=array(0,dim=rep(n,3))
		for(i in 1:n){
				for(j in 1:n){
						for(k in 1:n){
								prodxyz[i,j,k]=sum(prodx[i,]*prody[j,]*prodz[k,])
								#prodxyz[i,j,]=outer(prodx[i,],prody[j,])*prodz[k,])
						}
				}
		}
		
		# normalise so that the resultant density 
		# sums to 1 over range.
		prodxyz=prodxyz/sum(prodxyz)
		
		# I may wish to go on to normalise by the total arbour length
		# involved in the calculation / the volume in question to give
		# arbour density in microns / um3 = um^-2
		
		return(list(x = gx, y = gy, z = gz, dens=prodxyz))
}




smooth3d.old<-function (d, h, n = 25,oversample=4,nbins=n*oversample,lims = c(range(x), range(y),range(z))) 
{
		# prebin the data before doing a 3d smoothing
		
		r=apply(d,2,range)
		w=apply(r,2,diff)
		m=apply(r,2,mean)
		mind=apply(r,2,min)
		
		scaled_d=round(scale(d,cent=mind,scale=w/nbins))
		
		prebin=integer(nbins^3)
		dim(prebin)<-rep(nbins,3)
		
		for(i in 1:nrow(d)){
				x=scaled_d[i,]
				prebin[x[1],x[2],x[3]]=prebin[x[1],x[2],x[3]]+1
		}
		
		#k=2*2*2*oversample+1 #(kernel width)
		smoothed_scaled_d=GaussSmoothArray(prebin,voxdim=w/nbins,ksize=(round(2*2*2*oversample/(w[1]/nbins))%/%2)*2+1,sigma=diag(h,3))
		unscaled_smoothed_scaled_d=scale(d,cent=-mind,scale=nbins/w)

}

# remove unused functions
rm(smooth3d.old,kde3d)
