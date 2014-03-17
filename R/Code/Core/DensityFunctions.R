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
  .Deprecated('nat::imexpand.grid')
  if(!inherits(d,'im3d')){
    # not an im3d, so make a fake im3d
    # ensuring that we get dims properly set (as well as bounding box)
    dims=dim(d)
    if(is.null(dims)){
      dims=sapply(attributes(d)[c('x','y','z')], length)
    }
    d=im3d(dims=dims,BoundingBox=getBoundingBox(d))
  }
  nat::imexpand.grid(d)
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
	.Deprecated("nat::imscalebar")
	nat::imscalebar(levels=levels,col=col,nlevels=nlevels,zlim=zlim,
		horizontal=horizontal,lab=lab,mar=mar, border=NULL, ...)
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
  .Deprecated('nat::voxdims')
  if(!inherits(d,'im3d')){
    # not an im3d, so copy attributes to make a fake object
    d0=numeric()
    mostattributes(d0)=attributes(d)
    class(d0)=c('im3d',class(d0))
    d=d0
  }
  nat::voxdims(d)
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
		input
		halfVoxelDims=voxdim/2
		bounds[1,]=bounds[1,]+halfVoxelDims
		bounds[2,]=bounds[2,]-halfVoxelDims
		# zap small gets rid of FP rounding errors
		return(zapsmall(bounds))
	}else if(is.character(b)){
		return(nat::boundingbox(b))
	}
	return(NULL)
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
	z, zlim=range(z, finite = TRUE), xlim=NULL, ylim=NULL,
    plotdims=NULL,flipdims='y',filledContour=FALSE,asp=NA, axes=FALSE,
    xlab=NULL,ylab=NULL,
    nlevels=20,levels = pretty(zlim, nlevels+1),
    color.palette=jet.colors, col = color.palette(length(levels) - 1), ...){
	

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
	if(!inherits(z,'im3d')) class(z)=c('im3d',class(z))
	if(is.null(x)) x=attr(z,'x')
	if(is.null(y)) x=attr(z,'y')
	image(z, xlim = xlim, ylim = ylim, zlim = zlim, plotdims = plotdims, 
		flipdims = flipdims, filled.contour = filledContour, asp = asp, axes = axes, 
		xlab = xlab, ylab = ylab, nlevels = nlevels, levels = levels, 
		color.palette = color.palette, col = col, ...)
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
# MakeMaskFromDensity-> nat::threshold
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

xyzpos.gjdens<-function(d,ijk)
{
  .Deprecated('nat::xyzpos')
  if(!inherits(d,'im3d')){
    # not an im3d, so copy attributes to make a fake object
    d0=numeric()
    mostattributes(d0)=attributes(d)
    class(d0)=c('im3d',class(d0))
    d=d0
  }
  nat::xyzpos(d, ijk)
}

ijkpos.gjdens<-function(d,xyz,roundToNearestPixel=TRUE)
{
  .Deprecated('nat::ijkpos')
  if(!inherits(d,'im3d')){
    # not an im3d, so copy attributes to make a fake object
    d0=numeric()
    mostattributes(d0)=attributes(d)
    class(d0)=c('im3d',class(d0))
    d=d0
  }
  nat::ijkpos(d, xyz, roundToNearestPixel=roundToNearestPixel)
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
