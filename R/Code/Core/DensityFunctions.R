#@+leo-ver=4-thin
#@+node:jefferis.20051014173041.13:@thin R/DensityFunctions.R
#@@language r
# DensityFunctions.R
# 
# source(file.path(CodeDir,"DensityFunctions.R"))

# Helper functions for handling densities - started after the move
# to matlab based 3D density estimation

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

#@+others
#@+node:jefferis.20060522195714.1:unmask
unmask=function(d,mask,default=NA,attributes.=attributes(mask),copyAttributes=TRUE){
	rval=mask
	if(copyAttributes) attributes(rval)=attributes.
	rval[mask!=1]=default
	rval[mask==1]=d
	rval
}
#@nonl
#@-node:jefferis.20060522195714.1:unmask
#@+node:jefferis.20060522182249.2:DensityCumSumFracThreshold
DensityCumSumFracThreshold<-function(d,fractions,decreasing=TRUE,returnValue=TRUE,...){
# This function will calculate the threshold value for which the cumulative sum
# is greater than a given fraction of the total sum.
    sortedVals=sort(d,decreasing=decreasing,...)
    cs=cumsum(sortedVals)
    totalSum=sum(sortedVals)
    # fix so that can cope with more than one fraction
    #thresholdIdx=min(which(cs>=fractions*totalSum))
    thresholdIdxs=sapply(fractions,function (f) min(which(cs>=f*totalSum)))
    if(returnValue) return(sortedVals[thresholdIdxs])
    else return(thresholdIdxs)
}
#@nonl
#@-node:jefferis.20060522182249.2:DensityCumSumFracThreshold
#@+node:jefferis.20060314194700:expand.grid.gjdens
expand.grid.gjdens<-function(d){
    # takes the x,y and z attributes of d and
    # makes an n x 3 matrix containing the grid points
    
    # this is basically the guts of expand.grid
    # but returns a nice clean matrix
    
    dims=dim(d)
    orep <- prod(dims)
    nargs=3
	
	if(all(c("x","y","z") %in% names(attributes(d)))){
		args=attributes(d)[c("x","y","z")]
	} else {
		args=list()
		boundingBox=matrix(getBoundingBox(d),nrow=2)
		for(i in seq(dims)){
			args[[i]]=seq(from=boundingBox[1,i],to=boundingBox[2,i],length=dims[i])
		}
	}
    rep.fac <- 1
    rval=matrix(nrow=orep,ncol=length(dims))
    for (i in 1:nargs) {
        x <- args[[i]]
        nx <- length(x)
        orep <- orep/nx
        x <- x[rep.int(rep.int(seq(length = nx), rep.int(rep.fac, 
            nx)), orep)]
        if (!is.factor(x) && is.character(x)){
            cat("converting to factor")
            x <- factor(x, levels = unique(x))
        }
        rval[,i]=x
        rep.fac <- rep.fac * nx
    }
    return(rval)
}
#@nonl
#@-node:jefferis.20060314194700:expand.grid.gjdens
#@+node:jefferis.20060314184616.2:Interpolate3D
Interpolate3D<-function(d,newx,newy,newz,method=c("constant","linear")){
    method=match.arg(method)
    if(missing(newy) & missing(newz)){
        # assume that we have a gj density array
        if(is.array(newx)) newx=attributes(newx)
        # process if we have a list
        if(is.list(newx)) {
            newz=newx$z; newy=newx$y; newx=newx$x
        }
    }
    
    require(e1071) # library with primitive 3d interpolation function
    
    # make the grid for interpolation
    new.grid=data.matrix(expand.grid(x=newx,y=newy,z=newz))
    
    # NB adims is the grid positions of the original data, d
    newd=interpolate(new.grid,d,adims=attributes(d)[c("x","y","z")],method=method)
    
    # find the new bounding box
    BoundingBox=as.vector(sapply(list(newx,newy,newz),range))
    
    attributes(newd)<-c(attributes(newd),list(x=newx,y=newy,z=newz,BoundingBox=BoundingBox))
    dim(newd)=sapply(list(newx,newy,newz),length)
    return(newd)
}
#@nonl
#@-node:jefferis.20060314184616.2:Interpolate3D
#@+node:jefferis.20060314184616.1:NearestNeighbourInterpolate
NearestNeighbourInterpolate<-function(d,newx,newy,newz,method=c("constant","linear")){
    method=match.arg(method)
    if(missing(newy) & missing(newz)){
        # assume that we have a gj density array
        if(is.array(newx)) x=attributes(newx)
        # process if we have a list
        if(is.list(newx)) {
            newz=newx$z; newy=newx$y; newx=newx$x
        }
    }
    
    # now do the interpolation
    # for each intersection of the grid defined by x, y, z
    # find the corresponding nearest neighbour in d
    
    # INCOMPLETE - see Interpoalate3D     
}
#@nonl
#@-node:jefferis.20060314184616.1:NearestNeighbourInterpolate
#@+node:jefferis.20060305164352:IntegrateDensity
IntegrateDensity<-function(d){
    # first handle the case of a list of density objects
    if(is.list(d)) return(sapply(d,IntegrateDensity))
    # Now assume that we have a single density object
    ndims=length(dim(d))
   
	voxelDimensions=voxdim.gjdens(d)
	
	voxelVolume=prod(voxelDimensions)
	if(ndims==3) return( voxelVolume*sum(d) )
	if(ndims==2) {
		projdimlength=length(attributes(d)[[attr(d,"ProjDim")]])
		if(projdimlength>0){
			return(projdimlength * voxelVolume*sum(d) )
		}
	}
	else stop("IntegrateDensity only handles 2d or 3d data")
}
#@nonl
#@-node:jefferis.20060305164352:IntegrateDensity
#@+node:jefferis.20060302181445:makeScaleBar
makeScaleBar<-function(levels,col,nlevels=NULL,zlim=NULL,horizontal=TRUE,lab="Density",
    mar=c(4,2,2,2)+0.1,border=NULL, ...){
    
	# allow 
	if(!is.null(zlim) ){
		nc <- length(col)
		if ( (any(!is.finite(zlim)) || diff(zlim) < 0)) 
			stop("invalid z limits")
		if (diff(zlim) == 0) 
			zlim <- if (zlim[1] == 0) 
				c(-1, 1)
			else zlim[1] + c(-0.4, 0.4) * abs(zlim[1])
		#z <- (z - zlim[1])/diff(zlim)
		#zi <- floor((nc - 1e-05) * z + 1e-07)
		#zi[zi < 0 | zi >= nc] <- NA
		levels=seq(from=zlim[1],to=zlim[2],len=nc+1)
	}
    # allow scaleinfo objects to be passed directly
    if(missing(col) && is.list(levels)){
        col=levels$col
        levels=levels$levels
	}
	if(horizontal){
        par(mar=mar)
        plot(range(levels), c(0,1), type="n",
            xaxs="i", yaxs="i", xlab=lab, ylab="", yaxt="n", ...)
        rect(levels[-length(levels)], 0, border=border,
             levels[-1], col = col  , 1)
    } else {
        par(mar=mar[c(2,1,3,4)])
        plot(c(0,1), range(levels), type="n",
            xaxs="i", yaxs="i", xlab="", ylab=lab, xaxt="n", ...)
        rect(0, levels[-length(levels)], border=border,
             1, levels[-1], col = col)
    }
}
#@-node:jefferis.20060302181445:makeScaleBar
#@+node:jefferis.20051015010751:densityArrayFromList
densityArrayFromList<-function(d, Simplify=FALSE){
	if(is.list(d)){
		# needs to be turned into an array with dims
		#  x y z "t"
		subdims=dim(d[[1]])
		dd=do.call(cbind,d)
        # Simplify => return a matrix rather than array
		if(Simplify) dim(dd)<-c(prod(subdims),length(d)) 
            else dim(dd)<-c(subdims,length(d))
		a=attributes(d[[1]])[c("BoundingBox","x","y","z")]
		attributes(dd)=c(attributes(dd),a)
		ndims=length(dim(dd))
		dimnames(dd)[[ndims]]<-names(d)
        class(dd)<-c("gjdensityArray",class(dd))
		return (dd)
	}
	if(!is.array(d)) stop ("must pass a density list or array")
	if(!inherits(d,"gjdensityArray")) class(d)<-c("gjdensityArray",class(d))
    return (d)
}

densityListFromArray<-function(d, attributes.d=attributes(d)){
	if(is.array(d)){
		# needs to be turned into a list 
		#  x y z "t"
		dims=dim(d); ndims=length(dims)
		subdims=dims[-ndims]; lastdim=dims[ndims]
		lastdimnames=dimnames(d)[ndims]
		dim(d)=c(prod(subdims),lastdim)
		l=list()
		for(g in seq(lastdim)){
			l[[g]]=d[,g]
			if(!is.null(attributes.d)) mostattributes(l[[g]])=attributes.d
#			class(dd)<-c("gjdensityArray",class(dd))
			dim(l[[g]])=subdims
		}
		if(!is.null(lastdimnames)) names(l)=lastdimnames
		return(l)
	}
	if(!list(d)) stop ("must pass a density list or array")
	else d
}


#@-node:jefferis.20051015010751:densityArrayFromList
#@+node:jefferis.20051014173041.14:rgl.gjdens
rgl.gjdens<-function(d,x=attr(d,"x"),y=attr(d,"y"),z=attr(d,"z"),scaleFact=1,
    UseSprites=FALSE,...){
    # attempt to do some kind of 3D plot of image density
    # using jet colours?
    # plot dots everywhere above some threshold in white, with fixed size but variable alpha
    vd=as.vector(d)
    zmax=max(vd)
    vd=vd/zmax*scaleFact
    thresh=quantile(vd,0.9)

    # make a set of points corresponding to the 3d grid
    # in an array, left most index moves fastest
    nx=length(x)
    ny=length(y)
    nz=length(z)
    
    xyzpoints=data.matrix(expand.grid(z,y,x))
#    xyzpoints=data.matrix(expand.grid(attr(d,"x"),attr(d,"y"),attr(d,"z")))
    xyzpoints=xyzpoints[vd>thresh,]
    if(!UseSprites) rgl.points(xyzpoints[,1],xyzpoints[,2],xyzpoints[,3],alpha=vd[vd>thresh],...)		
    else rgl.sprites(xyzpoints[,1],xyzpoints[,2],xyzpoints[,3],alpha=vd[vd>thresh],
            textype='alpha',texture=system.file("textures/particle.png",package="rgl"),...)		
}
#@nonl
#@-node:jefferis.20051014173041.14:rgl.gjdens
#@+node:jefferis.20051014173041.16:Masks (TODO)
MaskDensity<-function(d,mask){
		# function to make take eg 3D LH density and mask it
		# so that all areas outside LH are set to 0
}

MakeMaskToFitBounds<-function(mask,bounds,n){
		# take a mask and interpolate it to fit 
		# 
}
is.gjdens<-function(d){
	if(inherits(d,"gjdens")) {
		return(TRUE)
	} else {
		if( ( inherits(d,"array") || inherits(d,"array") ) && !is.null(attr(d,"BoundingBox"))){
			return(TRUE)
		}		
	}
	return(FALSE)
}

	
voxdim.gjdens<-function(d){
	if(all(c("x","y","z") %in% names(attributes(d)))){
		originaldims=sapply(attributes(d)[c("x","y","z")],length)
	} else {
		originaldims=dim(d)
	}
	if (!is.null(attr(d,"bounds")))
		# bounds = outer limit of voxels
		return(diff(matrix(attr(d,"bounds"),nrow=2))/originaldims)
	else if (!  is.null(attr(d,"BoundingBox"))) {
		# BoundingBox = CENTRES of outer voxels (like Amira)
		# therefore includes 1 fewer voxel in each dimension
		return(diff(matrix(attr(d,"BoundingBox"),nrow=2))/(originaldims-1))
	} 
	#warning("Cannot find bounds or BoundingBox attribute")
	return(NULL)
}

getBounds<-function(b){
	# nb Bounds = outer edges of outer voxels 
	# (ie probably what most people expect & more like ImageJ than Amira)
	if(is.character(b)) {
		if(toupper(b)=='LH') return(c(95,165,60,130,0,70))
		if(toupper(b)=='MB') return(c(45,115,60,130,60,95))
		if(file.exists(b)){
			b=read.im3d(b,ReadData=FALSE)
		} else stop(paste("Don't know how to get bounds for",b))
	}
	if(!is.null(attr(b,"bounds"))) return(attr(b,"bounds"))
	else if(!is.null(attr(b,"BoundingBox"))){
		bounds<-matrix(attr(b,"BoundingBox"),nrow=2)
		# nb if any of the dimensions are 1 then the voxel dimension
		# cannot be calculated
		dims=dim(b)
		if(is.null(dims))
			dims=sapply(c('x','y','z'),function(d) length(attr(b,d)),USE.NAMES=F)
		halfVoxelDims=apply(matrix(bounds,nrow=2),2,diff)/(dims-1)/2
		bounds[1,]=bounds[1,]-halfVoxelDims
		bounds[2,]=bounds[2,]+halfVoxelDims
		# zap small gets rid of FP rounding errors
		return(zapsmall(bounds))
	}
	else return(NULL)
}

getBoundingBox<-function(b,bounds=attr(b,"bounds"),voxdim=voxdim.gjdens(b)){
	# nb BoundingBox = CENTRES of outer voxels (like Amira)
	if(!is.null(attr(b,"BoundingBox"))) return(attr(b,"BoundingBox"))
	else if(!is.null(bounds) && !is.null(voxdim)){
		
		if(is.vector(bounds)) bounds<-matrix(bounds,nrow=2)
		halfVoxelDims=voxdim/2
		bounds[1,]=bounds[1,]+halfVoxelDims
		bounds[2,]=bounds[2,]-halfVoxelDims
		# zap small gets rid of FP rounding errors
		return(zapsmall(bounds))
	}
	else if(is.character(b) && file.exists(b)){
		return(attr(read.im3d(b,ReadData=FALSE),'BoundingBox'))
	}
	return(NULL)
}

#' read 3d image using one of the available readers returning in gjdens format
read.im3d<-function(f,ReadData=TRUE,...){
	if(is.nrrd(f)) Read3DDensityFromNrrd(f,ReadData=ReadData,...)
	else stop("Unable to read data from: ",f)
}

#@nonl
#@-node:jefferis.20051014173041.16:Masks (TODO)
#@+node:jefferis.20051014173041.15:image.gjdens (TODO)
contour.gjdens<-function(x=NULL,y=NULL,
	z, xlim=range(x,finite=TRUE), ylim=range(y,finite=TRUE),
    plotdims=NULL,flipdims='y', axes=FALSE, drawlabels=FALSE, ...){
	
	if (missing(z)) {
		if (!is.null(x)) {
			z <- x
			attributes(z)<-attributes(x)
			x <- NULL
		}
		else stop("no 'z' matrix specified")
	}

	if(!is.null(plotdims)){
		plotdims=tolower(plotdims); plotdims=unlist(strsplit(plotdims,split=""))
		if(is.null(x)) x=attr(z,plotdims[1])
		if(is.null(y)) y=attr(z,plotdims[2])
	} else if(!is.null(attr(z,"ProjDim"))){
		# if this is a projection, then choose correct axes to display
		plotdims=setdiff(c("x","y","z"),attr(z,"ProjDim"))
		if(is.null(x)) x=attr(z,plotdims[1])
		if(is.null(y)) y=attr(z,plotdims[2])
	} else if (all( c("x","y")%in%names(attributes(z)) )){
		if(is.null(x)) x=attr(z,"x")
		if(is.null(y)) y=attr(z,"y")
	}
	# If we still haven't set anything, then use default
	if(is.null(x)) x=seq(0,1,len=nrow(z))
	if(is.null(y)) y=seq(0,1,len=ncol(z))
	if(is.null(plotdims)) plotdims=c("x","y")

	# what about transposing matrix if axes have been swapped?
	numbers=1:3;names(numbers)=letters[24:26]
	if(numbers[plotdims[1]]>numbers[plotdims[2]]){
		z=t(z)
	}
	if(is.numeric(flipdims)) flipdims=names(numbers(flipdims))
	if(is.character(flipdims)) {
		flipdims=unlist(strsplit(tolower(flipdims),split=""))
		# nb at this stage we assume z looks like our axes, so
		# we don't try to match axis names etc
		if(plotdims[1]%in%flipdims) z<-flip.array(z,1)
		if(plotdims[2]%in%flipdims) z<-flip.array(z,2)
	}
	
	# now reverse axes if required
	if(plotdims[1]%in%flipdims) x=-rev(x)
	if(plotdims[2]%in%flipdims) y=-rev(y)
	
	contour(x=x,y=y,z=z,xlim=xlim,ylim=ylim,axes=FALSE,
		xlab=plotdims[1],ylab=plotdims[2],drawlabels=drawlabels, ...)

	if(axes){
		axis(2,pretty(par("usr")[3:4]),abs(pretty(par("usr")[3:4])))
		axis(1,pretty(par("usr")[1:2]),abs(pretty(par("usr")[1:2])))
	}
}


image.gjdens<-function(x=NULL,y=NULL,
	z, zlim=NULL, xlim=range(x,finite=TRUE), ylim=range(y,finite=TRUE),
    plotdims=NULL,flipdims='y',filledContour=FALSE,asp=NA, axes=FALSE,
    xlab=NULL,ylab=NULL,
    nlevels=20,levels = pretty(zlim, nlevels+1),
    color.palette=jet.colors,col = color.palette(length(levels) - 1),...){
		
    # # function which will extend R's image function
    # by 
    # 1. allowing images to have inverted x or ylims.
    #    - What would be the best here - I thought of just checking whether
    #    xlim or ylim were reversed, but I think that an actual flip command
    #    is preferable.
    # 2. Defaulting to a suitable colour ramp
    #@    << handle z as first argument >>
    #@+node:jefferis.20051016235253.1:<< handle z as first argument >>
    if (missing(z)) {
    	if (!is.null(x)) {
    		z <- x
    		attributes(z)<-attributes(x)
    		x <- NULL
    	}
    	else stop("no 'z' matrix specified")
    }
    #@-node:jefferis.20051016235253.1:<< handle z as first argument >>
    #@nl
    #@    << set x y z >>
    #@+node:jefferis.20051024111842:<< set x y z >>
    # don't try and plot anything if we have malformed zlims
    if(is.null(zlim)) zlim=range(z,finite=TRUE)
    if(!all(is.finite(zlim))){
    	warning(paste("supplied zlim is not finite:",zlim))
    	zlim=c(0,0)
    }
    
    if(!is.null(plotdims)){
        plotdims=tolower(plotdims); plotdims=unlist(strsplit(plotdims,split=""))
        if(is.null(x)) x=attr(z,plotdims[1])
        if(is.null(y)) y=attr(z,plotdims[2])
    } else if(!is.null(attr(z,"ProjDim"))){
        # if this is a projection, then choose correct axes to display
        plotdims=setdiff(c("x","y","z"),attr(z,"ProjDim"))
        if(is.null(x)) x=attr(z,plotdims[1])
        if(is.null(y)) y=attr(z,plotdims[2])
    } else if (all( c("x","y")%in%names(attributes(z)) )){
        if(is.null(x)) x=attr(z,"x")
        if(is.null(y)) y=attr(z,"y")
    }
    # If we still haven't set anything, then use default
    if(is.null(x)) x=seq(0,1,len=nrow(z))
    if(is.null(y)) y=seq(0,1,len=ncol(z))
    if(is.null(plotdims)) plotdims=c("x","y")
    #@-node:jefferis.20051024111842:<< set x y z >>
    #@nl
    #@    << handle flipdims >>
    #@+node:jefferis.20051016235253:<< handle flipdims >>
    numbers=1:3;names(numbers)=letters[24:26]
    
    # what about transposing matrix if axes have been swapped?
    if(numbers[plotdims[1]]>numbers[plotdims[2]]){
        z=t(z)
    }
    
    if(is.numeric(flipdims)) flipdims=names(numbers(flipdims))
    if(is.character(flipdims)) {
        flipdims=unlist(strsplit(tolower(flipdims),split=""))
        # nb at this stage we assume z looks like our axes, so
        # we don't try to match axis names etc
        if(plotdims[1]%in%flipdims) z<-flip.array(z,1)
        if(plotdims[2]%in%flipdims) z<-flip.array(z,2)
    }
    
    # now reverse axes if required
    if(plotdims[1]%in%flipdims) x=-rev(x)
    if(plotdims[2]%in%flipdims) y=-rev(y)
    #@nonl
    #@-node:jefferis.20051016235253:<< handle flipdims >>
    #@nl
    
    
    if(filledContour){
	    plot(0,0,xlim,ylim,type='n',asp=asp,axes=FALSE,
        xlab=plotdims[1],ylab=plotdims[2],xaxs="i",yaxs="i",...)
        if (!is.double(z)) storage.mode(z) <- "double"
        .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
        col = col))
    } else {
        image(x=x,y=y,z=z,zlim=zlim,xlim=xlim,ylim=ylim,col=col,asp=asp,axes=FALSE,
        xlab=plotdims[1],ylab=plotdims[2],...)
    }
    if(axes){
        axis(2,pretty(par("usr")[3:4]),abs(pretty(par("usr")[3:4])))
        axis(1,pretty(par("usr")[1:2]),abs(pretty(par("usr")[1:2])))
	}
    # Return info that will be useful for creating scalebars
	invisible(list(zlim=zlim,nlevels.actual=length(levels),nlevels.orig=nlevels,
        levels=levels,colors=col))
}

# function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)), 
	# z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
	# zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), 
	# nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) - 
		# 1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
	# xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, 
	# ...) 
#@nonl
#@-node:jefferis.20051014173041.15:image.gjdens (TODO)
#@+node:jefferis.20051014173041.17:Masks
MakeMaskFromDensity<-function(d,bounds=attr(d,"bounds"),BoundingBox=attr(d,"BoundingBox"),threshold=-1){
		# function that makes a mask object of the same size as
		# a particular density and sets values to 0 or 1 depending
		# on whether they exceed a threshold
		
		# Threshold could perhaps also be a plugin function
		# to allow fancier setting of levels

		m=integer(length(d))
		dim(m)<-dim(d)
		if(missing(BoundingBox))
		attr(m,"BoundingBox")=BoundingBox
		attr(m,"BoundingBox")=bounds
		m[d>threshold]=1
		m
}
#@-node:jefferis.20051014173041.17:Masks
#@+node:jefferis.20051014173041.18:projection
projection<-function(a,projdim='z',projfun=c('integrate','mean','sum'),warn=F,na.rm=T,mask=NULL,...){
    ndims=length(dim(a))
    if(is.character(projdim)){
        projdim=tolower(projdim)
        projdim=which(letters==projdim)-which(letters=="x")+1
    }
	if(ndims<3) {
		if(warn) warning("3D arrays only in zproj - no z axis in array")
		return (a)
	}
	if(!is.null(mask)) a[mask==0]=NA

	if(is.character(projfun) && projdim==ndims) {
		# This will do fast sums over the last dimension for 3 funcs
		projfun=match.arg(projfun)
		# These functions are handled specially
		if( projfun=="sum" || projfun=="integrate" ){
			rval=rowSums(a,dims=ndims-1,na.rm=na.rm)
		} else if(projfun=="mean"){ 		
			rval=rowMeans(a,dims=ndims-1,na.rm=na.rm)
		}
	} else {
		# This is the more general routine
		if(is.character(projfun) && projfun=="integrate") projfun.fun=sum
		else projfun.fun=match.fun(projfun)
		# Now do the projection
		margins=setdiff(1:ndims,projdim)
		rval=apply(a,margins,projfun.fun,na.rm=na.rm,...)
	}
	if(is.character(projfun) && projfun=="integrate") {
		dx=voxdim.gjdens(a)[projdim]
		rval=rval*dx
	}

	# copy over attributes
	attributeNamesToCopy=setdiff(names(attributes(a)),names(attributes(rval)))
    attributes(rval)=c(attributes(rval),attributes(a)[attributeNamesToCopy])
	# ... and set the ProjDim to the correct letter
	projDimChar=letters[23+projdim]
    attr(rval,'ProjDim')=if(!is.na(projDimChar)) projDimChar else projdim
    rval
}
#@nonl
#@-node:jefferis.20051014173041.18:projection
#@+node:jefferis.20051014173041.19:flip.array
flip.array=function(a,flipdim='X'){
	ndims=length(dim(a))
	if(ndims>3) stop("Can't handle more than 3D arrays")
	
	if(is.character(flipdim)){
		flipdim=tolower(flipdim)
		# fixed a bug which prevented y axis flips
		flipdim=which(letters==flipdim)-which(letters=="x")+1
	}
	if(!flipdim%in%seq(len=length(dim(a)))){
		stop("Can't match dimension to flip")
	}
	
	revidxs=dim(a)[flipdim]:1
	
	if(ndims==3){
		if(flipdim==1) rval=a[revidxs,,]
		if(flipdim==2) rval=a[,revidxs,]
		if(flipdim==3) rval=a[,,revidxs]
	}
	else if(ndims==2){
		if(flipdim==1) rval=a[revidxs,]
		if(flipdim==2) rval=a[,revidxs]
	}
	else if (ndims==1){
		return(flip.vector(a))
	}
        attributes(rval)=attributes(a)
	return (rval)
}


flip.vector=function(x) rev(x)

flip.matrix=function(x,...) {
		flip.array(x,...)
}
#@-node:jefferis.20051014173041.19:flip.array
#@+node:jefferis.20051015023212:blendDensities.gjdensityArray
blendDensities.gjdensityArray<-function(dens, mat, scaleFact,
    RemoveMissingRows=TRUE,RemoveMissingCols=TRUE ) {
    # function to take an array of named densities and a
    # matrix of named cols (corresponding to the named densities)
    # and rows corresponding to the named output groupings
    # each entry in the matrix says how much of density j will
    # end up in ouput density i
    # 
    # NB the output will be an array with the dimensions 1:n-1
    # identical to the original density array and last dimension
    # equal to the number of cols in mat
    
    if(RemoveMissingCols){
        # remove any empty cols from mixmat
        emptyCols=colSums(mat,na.rm=T)==0
        mat=mat[,!emptyCols]
    }
    if(RemoveMissingRows){
        # remove any empty rows rom mixmat
        emptyRows=rowSums(mat,na.rm=T)==0
        mat=mat[!emptyRows,]
    }
    if(nrow(mat)==0) return(NULL)

    if (missing(scaleFact)){
        # How about changing the scaleFact so that 
        scaleFact=rep(1,ncol(mat))
        names(scaleFact)=colnames(mat)
    }		
    if (!is.matrix(mat)) mat=as.matrix(mat)
    
		ndims=length(dim(dens))
		densColNames=dimnames(dens)[[ndims]]
    colsToMix=intersect(colnames(mat),densColNames)    
    
    if( !all(colsToMix%in%names(scaleFact)) ){
        stop("scaleFact must contain named scale factors for each entry in denslist that will be mixed")
    }
    #print(scaleFact)
    
    # reduce the number of dimensions to allow multiplication
    subdims=dim(dens)[-ndims]
    dim(dens)=c(prod(subdims),dim(dens)[ndims])
    # put back the names of last dimension ie colnames
    colnames(dens)<-densColNames
    
    densarrayout=dens[,colsToMix]%*%t(mat[,colsToMix])
    mostattributes(densarrayout)<-attributes(dens)
    dim(densarrayout)<-c(subdims,nrow(mat))
    class(densarrayout)<-c("gjdensityArray",class(densarrayout))
	attr(densarrayout,"mixmat")=mat
	attr(densarrayout,"mixedCols")=colsToMix
    densarrayout
}
#@nonl
#@-node:jefferis.20051015023212:blendDensities.gjdensityArray
#@+node:jefferis.20051014173041.20:blendDensities
blendDensities<-function(denslist, mat, scaleFact,
    RemoveMissingRows=TRUE,RemoveMissingCols=TRUE ) {
    # function to take a list of named densities and a
    # matrix of named cols (corresponding to the named densities)
    # and rows corresponding to the named output groupings
    # each entry in the matrix says how much of density j will
    # end up in ouput density i
    # 
    # NB the output will be a list of identically sized densities
    # the list will be as long as the nrows of the matrix
    
    if(RemoveMissingCols){
        # remove any empty cols from mixmat
        emptyCols=colSums(mat,na.rm=T)==0
        mat=mat[,!emptyCols]
    }
    if(RemoveMissingRows){
        # remove any empty rows rom mixmat
        emptyRows=rowSums(mat,na.rm=T)==0
        mat=mat[!emptyRows,]
    }
    if(nrow(mat)==0) return(NULL)

    if (missing(scaleFact)){
        # How about changing the scaleFact so that 
        scaleFact=rep(1,ncol(mat))
        names(scaleFact)=colnames(mat)
    }		
    if (!is.matrix(mat)) mat=as.matrix(mat)
    
    colsToMix=intersect(colnames(mat),names(denslist))
    
    
    if( !all(colsToMix%in%names(scaleFact)) ){
        stop("scaleFact must contain named scale factors for each entry in denslist that will be mixed")
    }
    
    denslistout=list()		
    for (od in rownames(mat)){
        for( g in colsToMix){
                #cat("g=",g,"od=",od,"\n")
            if(length(denslistout[[od]])==0){
                denslistout[[od]]=mat[od,g]*denslist[[g]]*scaleFact[g]
            } else {
                denslistout[[od]]=denslistout[[od]]+mat[od,g]*denslist[[g]]*scaleFact[g]
            }
            
            if(any(is.na(denslistout[[od]]))){
                #cat("Gone NA!\n")
            }
            
        }
	}
	attr(denslistout,"mixmat")=mat
	attr(denslistout,"mixedCols")=colsToMix
    denslistout
}
#@nonl
#@-node:jefferis.20051014173041.20:blendDensities
#@+node:jefferis.20051031112814:getMBTips
# function to find the tips of the mushroom body collaterals of a neuron

getMBTips<-function(ANeuron,verbose=TRUE){
	# Just get all the points of all the segments
	# and see which ones are endpoints
	allMBPoints=unique(unlist(ANeuron$SegList[unlist(ANeuron$MBSegNos)]))
	MBTips=intersect(ANeuron$EndPoints,allMBPoints)
	if(verbose && length(MBTips)<1) warning("This neuron has no MB tips")
	return(MBTips)
}
#@-node:jefferis.20051031112814:getMBTips
#@-others
#@+at
# 
#@-at
#@-node:jefferis.20051014173041.13:@thin R/DensityFunctions.R
#@-leo

clampmax<-function(xmin,xmax) {
	# this fn returns a new function that will find the maximum of its inputs
	# and then clamp the return value between xmin and xmax
	# +/- Inf are converted to NA
	if(missing(xmax)) {
		xmax=xmin[2]
		xmin=xmin[1]
	}
	# Example: image.gjdens(projection(d,projfun=clampmax(0,15)))
	function(x,...){
		r=max(x,...)
		if(r==Inf || r==-Inf) 
			NA
		else if(r<xmin)
			xmin 
		else if(r>xmax)
			xmax
		else r
	}
}

xyzpos.gjdens<-function(d,ijk)
{
	# return the xyz position for a pixel location (i,j,k)
	# This will be the pixel centre based on the bounding box
	# Note that ijk will be 1-indexed according to R's convention

	# transpose if we have received a matrix (with 3 cols i,j,k) so that
	# multiplication below doesn not need to be changed
	if(is.matrix(ijk)) ijk=t(ijk)
	if(any(ijk<1)) warning("expects 1-indexed pixel coordinate so pixels <1 make little sense")
	dxyz=as.vector(voxdim.gjdens(d))
	origin=getBoundingBox(d)[c(1,3,5)]
	xyz=(ijk-1)*dxyz+origin
	if(is.matrix(xyz)) t(xyz) else xyz
}

ijkpos.gjdens<-function(d,xyz,roundToNearestPixel=TRUE)
{
	# return the ijk position for a physical location (x,y,z)
	# This will be the pixel centre based on the bounding box
	# Note that ijk will be 1-indexed according to R's convention

	# transpose if we have received a matrix (with 3 cols x,y,z) so that
	# multiplication below doesn not need to be changed
	if(is.matrix(xyz)) xyz=t(xyz)

	dxyz=as.vector(voxdim.gjdens(d))
	BB=getBoundingBox(d)
	origin=BB[c(1,3,5)]
	farcorner=BB[c(2,4,6)]

	ijk=(xyz-origin)/dxyz+1
	if(roundToNearestPixel) {
		ijk=round(ijk)
		if(any(ijk<1) || any(ijk>dim(d))) warning("pixel coordinates outside image data")
	}
	if(is.matrix(ijk)) t(ijk) else ijk
}

MinBoundingBox<-function(d,threshold=0,aspixels=TRUE)
{
	# Returns the bounding box that contains all voxels above threshold
	# find voxels above threshold
	z=d>threshold
	# Make a list containing projections down onto each axis
	l=lapply(seq(length(dim(d))), function(x) apply(z,x,sum,na.rm=TRUE))
	newBBAsPixels=sapply(l,function(x) range(which(x>0)))
	if(aspixels) newBBAsPixels
	else xyzpos.gjdens(d,newBBAsPixels)
}

CropToBoundingBox<-function(d,bb)
{
	# Crops a density to the nearest pixels mataching a specified bounding box 
	# Can get the bounding box from another density
	if(is.gjdens(bb)) bb=getBoundingBox(bb)
	if(!is.matrix(bb)) matrix(bb,ncol=3)
	if(!all(dim(bb)==c(2,3)))
		stop("Incorrect bounding box specification (x0,x1,y0,y1,z0,z1)")
	ijks=ijkpos.gjdens(d,bb)
	newd=d[ijks[1]:ijks[2],ijks[3]:ijks[4],ijks[5]:ijks[6]]
	dims=dim(newd)
	# This seems to unset dims - don't know why
	mostattributes(newd)=attributes(d)
	dim(newd)<-dims
	actualbb=xyzpos.gjdens(d,ijks)
	attr(newd,"BoundingBox")=actualbb
	newd
}
