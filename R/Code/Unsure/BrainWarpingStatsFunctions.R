# BrainWarpingStatsFunctions.R
# Functiosn for Displacement and Volume Comparisons
# including across the sexes for the amount of warping 
# for each brain

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

# Need to write functions to
# calculate deformations at given points - probably voxel centres for
# MB/LH densities

# source(file.path(CodeDir,"BrainWarpingStatsFunctions.R"))

transformedPoints=function(Brain=NULL,xyzs=NULL,
	warpfile=NULL,AllRegDir=get('AllRegDir',envir=.GlobalEnv),
	gregxform=file.path(cmtk.bindir(check=TRUE),"gregxform"),direction=c("inverse","forward"),
	transforms=c("warp","affine"),gregxoptions="-b"){
	if(file.access(gregxform,mode=1)<0 || file.info(gregxform)$isdir)
		stop(paste("Unable to access gregxform executable at:",gregxform))
	
	if(any(grep("\\-b",gregxoptions))) binary=TRUE else binary=FALSE
	gregxform=paste(gregxform,gregxoptions)
	
	transforms=match.arg(transforms,several.ok=TRUE)
	direction=match.arg(direction) #nb inverse implies from sample to ref

	# massage xyzs input to a 3 col matrix
	if(is.data.frame(xyzs)) xyzs=data.matrix(xyzs)
		if(ncol(xyzs)>3){
		if(all(c("X","Y","Z")) %in% colnames(xyzs))
		xyzs=xyzs[,c("X","Y","Z")]
		else xyzs=xyzs[,1:3]
	}

	tmpfile=tempfile(); tmpfile2=tempfile()
	
	if(binary) writeBin(as.vector(t(xyzs)),tmpfile,size=4) 
	  else write.table(xyzs,col.names=F,row.names=F,file=tmpfile)

	#gregxform=file.path("/Users/jefferis/bin/gregxform")
	if(direction=="forward") gregxform=paste(gregxform,"--forward")

	if(is.null(warpfile)){
		# First check if we have an entry in MNInfo
		if(is.null(warpfile<-MNInfo$WarpFile[MNInfo$Brain==Brain])){
			# if not use traceinfo
			if(!exists("TraceInfo")) {
				warning("You must supply a warpfile, or make sure that there is an entry in the MNInfo$WarpFile column")
				return(NULL)
			}			
			warpfile=file.path(RootDir,"allreg",TraceInfo$StudyList[TraceInfo$Brain==Brain])
		}
	} 
	# check the warpfile name
	if(isTRUE(charmatch("registration",basename(warpfile))==1))
		warpfile=dirname(warpfile)
	if(!file.exists(warpfile)){
		# try and see if we can find the warpfile in another directory
		warpfile=file.path(AllRegDir,basename(warpfile))
		if(!file.exists(warpfile)) {
			warning(paste("unable to locate warpfile:",warpfile))
			return(NULL)
		}
	}

	l=list(pre=xyzs)
	if("affine"%in%transforms){
		# afflistfile=sub(file.path("^(.*)","warp","(.*)_warp_.*\\.list"),file.path("\\1","affine","\\2_9dof.list"),warpfile)
		# now use new affine switch of gregxform
		system(paste(gregxform,"--affine",warpfile,"< ",tmpfile,">",tmpfile2),
			intern=FALSE,ignore.stderr=FALSE)
		if(binary){
			# cat("Reading Binary data\n")
			l$affine=matrix(readBin(tmpfile2,"numeric",size=4,n=nrow(xyzs)*3),ncol=3,byrow=T)
			attributes(l$affine)=attributes(xyzs)
		} else {
			l$affine=read.table(tmpfile2,col.names=c("X","Y","Z"))
		}

	}
	
	if("warp"%in%transforms){
		system(paste(gregxform,warpfile,"< ",tmpfile,">",tmpfile2),
			intern=FALSE,ignore.stderr=TRUE)
		if(binary){
			# cat("Reading Binary data\n")
			l$warp=matrix(readBin(tmpfile2,"numeric",size=4,n=nrow(xyzs)*3),ncol=3,byrow=T)
			attributes(l$warp)=attributes(xyzs)
		} else {
			l$warp=read.table(tmpfile2,na.strings="ERR",col.names=c("X","Y","Z"))
		}
	}
	unlink(c(tmpfile,tmpfile2))
	return(l)
}

TransformNeuron<-function(neuron,warpfile=NULL,transform=c("warp","affine"),...){
	transform=match.arg(transform,several.ok=FALSE)
	tps=transformedPoints(neuron$NeuronName,xyzmatrix(neuron),
	  warpfile=warpfile,transforms=transform,...)[[transform]]
	if(!is.null(tps)) xyzmatrix(neuron)<-tps
	else {
		warning("Failed to transform neuron")
		return(NULL)
	}
	return(neuron)
}

TransformNeuronSimple<-function(neuron,transform=c("original","affine")){
	# assumes that we are given a transformed neuron and want the original 
	# back again, or the affine applied to the original
	transform=match.arg(transform,several.ok=FALSE)
	if(!is.neuron(neuron)) neuron=GetNeuron(neuron)
	x<-TransformNeuron(neuron,direction='forward',transform='warp')
	if(transform=="original") return(x)
	if(transform=="affine") {
		return(TransformNeuron(x,direction='inverse',transform='affine'))
	}
}

TransformSurfFile<-function(surffile,outfile,warpfile=NULL,transform=c("warp","affine"),
	na.action=c("ignore","affine","warp-nocheck"),...){
	transform=match.arg(transform)
	na.action=match.arg(na.action)
	surfheaderlines<-readLines(surffile,n=500)
	vertexdef=grep("^Vertices",surfheaderlines)
	nVertices=as.integer(sub("Vertices\\s+","",surfheaderlines[vertexdef],perl=T))
	if(is.na(nVertices) || nVertices<1) stop("could not parse file")
	xyz=read.table(surffile,skip=vertexdef,nrows=nVertices)
	txyz=transformedPoints(surffile,xyz,warpfile=warpfile,transforms=transform,...)[[transform]]
	if(transform=="warp" && na.action!="ignore" && any(is.na(txyz))) {
		narows=apply(txyz,1,function(x) any(is.na(x)))
		if(na.action=='affine'){
			# use affine transformed points to replace NAs
			txyz[narows,]=transformedPoints(surffile,xyz[narows,],warpfile=warpfile,
				transforms='affine',...)[['affine']]
			warning("Falling back to affine transformation for ",sum(narows),"/",length(narows),
			" points.")
		}
		else{
			# use warp without gregxform checking if inversion was successful
			txyz[narows,]=transformedPoints(surffile,xyz[narows,],warpfile=warpfile,
				transforms='warp',gregxoptions="--binary --no-check",...)[['warp']]
			warning("Falling back to warp-nocheck transformation for ",	sum(narows),"/",length(narows),
				" points.")
		}
	}
	# now splice file back together again
	writeLines(surfheaderlines[1:vertexdef],outfile)
	write.table(txyz,outfile,append=TRUE,col.names=FALSE,row.names=FALSE)
	# now we need to append the rest of the file
	cmd=paste("tail -n +",vertexdef+nVertices," ",shQuote(surffile)," >> ",shQuote(outfile),sep="")
	system(cmd)
}

#' Mirror neuron across midplane of axis, optionally correcting for asymmetry
#'
#' @param n A neuron
#' @param mirrorAxisSize The size of the axis to mirror across. If vector of 
#'  length 2 then assume this represents min and max bouding box for this axis.
#' @param mirrorAxis One of X,Y or Z. Can also be a gjdens object or 3d image file
#' @warpfunction A function that corrects for any asymmetry after s simple
#'  mirroring.
#' @details If mirrorAxis specifies an image file or a gjdens object 
#'  (ie loaded 3d image data) then \code{getBoundingBox} will be used to find
#'  the appropriate min and max bounding box values for the 
#' @details If there is a non-zero origin, we want to sum min and max bounds 
#'  to get the equivalent of midline pos x 2 for the mirroring since
#'  newpos = midline_pos * 2 - oldpos == max_bound + min_bound - oldpos
#' @seealso \code{\link{transform.neuron},\link{TransformNeuron},
#'  \link{getBoundingBox}}
#' @export
MirrorNeuron<-function(n,mirrorAxisSize,mirrorAxis=c("X","Y","Z"),warpfunction,...){
	mirrorAxis=match.arg(mirrorAxis)
	if(!is.numeric(mirrorAxisSize)){
		# assume that this represents an object or file that we want to use to find
		# bounding box
		bb=matrix(getBoundingBox(mirrorAxisSize),
			ncol=3,dimnames=list(NULL,c("X","Y","Z")))
		mirrorAxisSize=bb[,mirrorAxis]
	} 
	# if there is a non-zero origin, we want to sum min and max bounds to get
	# the equivalent of midline pos x 2 for the mirroring since
	# newpos = midline_pos * 2 - oldpos
	if(length(mirrorAxisSize)==2) mirrorAxisSize=sum(mirrorAxisSize)
	
	n$d[,mirrorAxis]=mirrorAxisSize-1*n$d[,mirrorAxis]
	if(!missing(warpfunction))
	warpfunction(n,...)
	else n
}

CompareNeuronTransformations<-function(n,transforms=c("original","affine","warp"),...){
	transforms=match.arg(transforms,several.ok=TRUE)
	if(any(transforms=="original")) plotneuron3d(TransformNeuronSimple(n,"original"), Col='red',Clear=F)
	if(any(transforms=="affine")) plotneuron3d(TransformNeuronSimple(n,"affine"), Col='green',Clear=F)
	# warped version (ie that in MyNeurons)
	if(any(transforms=="warp")) plotneuron3d(n,Col='blue',Clear=F)
}


calcOriginalPos<-function(Brain=NULL,xyzs=NULL,...){
	# returns a list containing the position of a grid in reference space
	# after transformation to the position only post affine and after warp
	OriginalPos=transformedPoints(Brain=Brain,xyzs=xyzs,direction="inverse",transform="warp",...)$warp
	# remove any NAs - gregxform seems to choke on NAs for some reason
	NotNARows=!is.na(OriginalPos[,1])
	OriginalPos2=OriginalPos[NotNARows,]
	Warp.result=(transformedPoints(Brain=Brain,OriginalPos2,transform="warp",direction="forward"))$warp
	Affine.result=(transformedPoints(Brain=Brain,OriginalPos2,transform="affine",direction="forward"))$affine
	Warp<-(Affine<-OriginalPos)
	Warp[NotNARows,]=Warp.result
	Affine[NotNARows,]=Affine.result
	DeltaWarpAffine=Warp-Affine
# 	
# 	DeltaWarpAffine.result=(transformedPoints(Brain=Brain,OriginalPos2,transform="warp",direction="forward"))$warp-
# 		(transformedPoints(Brain=Brain,OriginalPos2,transform="affine",direction="forward"))$affine
# 	# This little shuffle is to get NAs in empty slots
# 	DeltaWarpAffine=OriginalPos 

	return(list(pre=OriginalPos,post=xyzs,warp=Warp,affine=Affine,delta=DeltaWarpAffine))
}


findJacobian<-function(Brain=NULL,xyzs=NULL,warplistfile=NULL,
	gregxform=file.path(cmtk.bindir(check=TRUE),"gregxform"),gregxoptions="-j -b"){
	
	gregxform=paste(gregxform,gregxoptions)
	if(any(grep("\\-b",gregxoptions))) binary=TRUE else binary=FALSE

	tmpfile=tempfile(); tmpfile2=tempfile()
	
	if(binary) writeBin(as.vector(t(xyzs)),tmpfile,size=4) else write.table(xyzs,col.names=F,row.names=F,file=tmpfile)
	
	if(is.null(warplistfile)){
		# we haven't been supplied with a warp file try and find one
		warplistfile=findWarpFileFromBrainName(Brain)
	}
	# bail out if we can't find the registration file
	if(!isTRUE(file.exists(warplistfile))) {
		warning(paste("Unable to find a find a registration file for Brain =",Brain))
		return(NULL)
	}
	
	system(paste(gregxform,warplistfile,"< ",tmpfile,">",tmpfile2),
				intern=FALSE,ignore.stderr=FALSE)
	if(binary){
		# cat("Reading Binary data\n")
		jacobians=readBin(tmpfile2,"numeric",size=4,n=nrow(xyzs))
	} else {
		jacobians=scan(tmpfile2,na.strings="ERR",quiet=TRUE)
	}
	unlink(c(tmpfile,tmpfile2))
	return(jacobians)
}

findJacobianVolume<-function(...){
	jj=findJacobian(...)
	sideLength=round(length(jj)^(1/3))
	if( (sideLength*sideLength*sideLength)!=length(jj)){
		stop("Cannot make an equal sided cubic volume")
	}
	if(!is.null(jj)) dim(jj)<-rep(sideLength,3)
	jj
}

findGlobalScaling<-function(Brain=NULL,warplistfile=NULL,
	gregxform=file.path(cmtk.bindir(check=TRUE),"gregxform"),gregxoptions="-g"){
	
	if(is.null(warplistfile)){
		# we haven't been supplied with a warp file try and find one
		warplistfile=findWarpFileFromBrainName(Brain)
	}
	# bail out if we can't find the registration file
	if(!isTRUE(file.exists(warplistfile))) {
		warning(paste("Unable to find a find a registration file for Brain =",Brain))
		return(NULL)
	}

	globalScaling=as.numeric(system(paste(gregxform,gregxoptions,warplistfile),intern=TRUE,ignore.stderr=FALSE))
			
	return(globalScaling)
}		

findWarpFileFromBrainName<-function(Brain){
	foundBrain = pmatch(Brain,as.character(MNInfo$Brain))
	
	if(length(foundBrain)==1){
		warpfile<-MNInfo$WarpFile[foundBrain]
	} else {
		# if not use traceinfo
		if(!exists("TraceInfo")) {
			warning("You must make sure that there is an entry in the MNInfo$WarpFile column")
			return(NULL)
		}			
		warpfile=file.path(AllRegDir,TraceInfo$StudyList[TraceInfo$Brain==Brain])
	}

	# check the warpfile name
	if(isTRUE(charmatch("registration",basename(warpfile))==1))	warpfile=dirname(warpfile)
	if(!file.exists(warpfile)) {
		warning(paste("unable to locate warpfile:",warpfile))
		return(NULL)
	}
	return(warpfile)
}

findWarpFileFromBrainNames<-function(Brain){
	# this depends on the TraceInfo dataframe, which I have so far not
	# released.
	rownames(TraceInfo)=as.character(TraceInfo$Brain)
	filenames=TraceInfo[Brain,"Path.reg"]
	filenames[!is.na(filenames)]=file.path(AllRegDir,filenames[!is.na(filenames)])
	return (filenames)
}

findGlobalScaling.file<-function(filename){
	if(file.info(filename)$isdir){
		# this is a directory, so see if we can find the registration
		dirname=filename
		filename=dir(dirname,patt="^registration(\\.gz){0,1}",full.names=T)[1]
		if(is.na(filename)) 
			stop(paste("Unable to read registration file in",dirname))
	}

	if(isTRUE(grep("\\.gz$",filename)==1)){
		warpfile=gzfile(filename)
	} else warpfile=file(filename)
	first20lines=readLines(warpfile,20)
	close(warpfile)
	scalelines=grep('^[\\t\\s]*scale',first20lines,perl=T)
	if(length(scalelines)>1) scaleline=17
	else if(length(scalelines)==1) scaleline=scalelines[1]
	else return (NA)
	scaleterms=unlist(strsplit(first20lines[scaleline]," "))[2:4]
	prod(as.numeric(scaleterms))
}

VolumeChange<-function(xyzs,...){
	# find the vertices of tetrahedra centered on grid points
	newxyzs=t(apply(xyzs,1,RegularTetrahedron))
	dim(newxyzs)<-c(nrow(xyzs)*4,3)
	# 
	calcOriginalPos(xyzs=newxyzs,...)
}

makeGrid<-function(bounds="LH",type=c('centres','margins'),spacing=1.4) {
	if(is.character(bounds)) bounds=getBounds(bounds)
	else if(inherits(bounds,"array")){
		if(all(c("x","y","z") %in% names(attributes(bounds)))){
			a=attributes(bounds)
			grid=expand.grid(a$x,a$y,a$z)
			colnames(grid)=c("X","Y","Z")
			return(grid)			
		}
		else bounds=getBounds(bounds)
	}
	if(length(spacing)==1 && ncol(bounds)>1) spacing=rep(spacing,ncol(bounds))
	corrn=spacing/2
	type=match.arg(type)
	if(type=="margins"){
		grid=expand.grid(
			seq(bounds[1],bounds[2],by=spacing[1]),
			seq(bounds[3],bounds[4],by=spacing[2]),
			seq(bounds[5],bounds[6],by=spacing[3]) )
	} else {
		grid=expand.grid(
			seq(bounds[1]+corrn[1],bounds[2]-corrn[1],by=spacing[1]),
			seq(bounds[3]+corrn[2],bounds[4]-corrn[2],by=spacing[2]),
			seq(bounds[5]+corrn[3],bounds[6]-corrn[3],by=spacing[3]) )
	}
	grid=do.call(cbind,grid)
	colnames(grid)=c("X","Y","Z")
	return(grid)
}


makeXYZCube=function(d){
	dims=dim(d)
	ndims=length(dims)
	if(ndims>=3) return(d)
	# are there a cubic number of rows
	sideLength=dims[1]^(1/3)
	sideLength.int=round(sideLength)
	cat(sideLength)
	if(isTRUE(all.equal(sideLength,sideLength.int))){
		if(is.data.frame(d)) d=data.matrix(d)
		if(ndims==1){
			dim(d)<-rep(sideLength.int,3)
			return(d)
		}
		if(ndims==2){
			dim(d)<-c(rep(sideLength.int,3),dims[2])
			return(d)
		} else {
			stop("Can't handle dims>2")
		}
	} else {
		warning("Unable to make cube")
		return(d)
	}
}

simple.ttest<-function(x,y){
# 	float df = nValuesX + nValuesY - 2;
# 	float svar=((nValuesX-1)*varianceX + (nValuesY-1)*varianceY) / df;
# 	t = (avgX - avgY) / sqrt( svar * (1.0/nValuesX + 1.0/nValuesY));
	nx=length(x);	ny=length(y)
	df=nx+ny-2
	svar=((nx-1)*var(x) + (ny-1)*var(y)) / df
	(mean(x) - mean(y)) / sqrt( svar * (1.0/nx + 1.0/ny))
}

fast.ttest<-function(x,y){
	# where x and y are matrices
	# each row is an observation point and each column an individual
	# nrow(x) and y must be equal
	if( (nrows<-nrow(x)) != nrow(y)) stop ("Number of rows in x and y must be equal")
	nx=ncol(x);	 ny=ncol(y)
	df=nx+ny-2
	varsx=rowMeans(x*x) - rowMeans(x)^2;  varsy=rowMeans(y*y) - rowMeans(y)^2
	svars=(nx*varsx + ny*varsy) / df
	rval=(rowMeans(x) - rowMeans(y)) / sqrt( svars * (1.0/nx + 1.0/ny))
	return(rval)
}

fast.sd<-function(x,biased=FALSE) {
	# use built-in sd except for matrices
	if(!is.matrix(x)) return(sd(x))
	# for matrices where 
	# each row is an observation point and each column an individual
	n=nrow(x); oneovernminusone=1/(n-1)
	xbar=colMeans(x)
	sqrt(  ( colSums(x*x) - n*(xbar*xbar) )*(oneovernminusone)  )
}

permtvals.simple=function(jacdata,males,females){
	both=c(males,females)
	males.perm=sample(both,length(males))
	females.perm=setdiff(both,males.perm)
	tvals=fast.ttest(jacdata[,males.perm],jacdata[,females.perm])
	attr(tvals,"perm")<-list(males=males.perm,females=females.perm)
	return(tvals)
}	

parallelpermtvals.simple=function(jacdata,males,females,n){
	# Expects just one vector each of male and female data
	# and a number of permutations to run
	permdata=t(replicate(n,sample(jacdata)))
	fast.ttest(permdata[,males],permdata[,females])
#	attr(tvals,"perm")<-list(males=males.perm,females=females.perm)
#	return(tvals)
}
permtvals=function(jacdata,males,females,mask,standardAttributes){
	both=c(males,females)
	males.perm=sample(both,length(males))
	females.perm=setdiff(both,males.perm)
	
	
	if(missing(mask) && missing(standardAttributes)){
		standardAttributes=structure(list(dim = c(50, 50, 50), 
			BoundingBox = c(95.7, 164.3, 60.7, 129.3, 0.7, 69.3),
			x = seq(from=95.7,to=164.3,by=1.4),
			y = seq(from=60.7,to=129.3,by=1.4),
			z = seq(from=0.7,to=69.3,by=1.4)),
		.Names = c("dim", "BoundingBox", "x", "y", "z"))
	} else if( !missing(mask) && missing(standardAttributes) ) {
		standardAttributes=attributes(mask)
	} else if( missing(mask)  && !missing(standardAttributes) ){
		mask=array(1,dim=standardAttributes$dim)
	}
	tvals=array(NA,dim=standardAttributes$dim)
	tvals[mask==1]=fast.ttest(jacdata[,males.perm],jacdata[,females.perm])
	
	attributes(tvals)=standardAttributes
	attr(tvals,"perm")<-list(males=males.perm,females=females.perm)
	return(tvals)
}

ExtremeTvals<-function(jacdata,males,females,perms,quantiles=c(0,0.001,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.999,1)){
	# Try to make null distributions of test statistics for
	# Deformation-Based Morphometry data
	# perms is either the number of permutations
	# or a matrix of permutations with nMales+nFemale cols
	# and nPerms rows
	findPerm<-function(males,females){
		both=c(males,females)
		males.perm=sample(both,length(males))
		females.perm=setdiff(both,males.perm)
		c(males.perm,females.perm)
	}
	
	if(length(perms)==1) perms=t(replicate(perms,findPerm(males,females)))
	n=nrow(perms)
	maleIdxs=seq(males);femaleIdxs=seq(females)+length(males)
	
	evdtvals=t(apply(perms,1,function(x) 
		quantile(fast.ttest(jacdata[,x[maleIdxs]],jacdata[,x[femaleIdxs]]),quantiles,na.rm=TRUE) ))
	attr(evdtvals,'perms')=perms
	return(evdtvals)
}

# Read in the sex information from the PN2 FMPro database
# and merge it with info from spreadsheet
MergeSexInfo<-function(x){
	d=read.csv(file.path(NotesDir,"PN2Slides.txt"))
	d$SlideKey2 =factor(gsub("^A(.*)","\\1",as.character(d$SlideKey)))
	x$SlideKey=gsub("([A-Z]{2,3})[0-9]{1,3}[LRTB]","\\1",as.character(x$Brain))
	y=merge(x,d[,c("SlideKey2","Sex")],by.x="SlideKey",by.y="SlideKey2",all.x=TRUE,suffixes=c("",".db"))
	# Clean up the factor levels (as.char allows easy addition of a new level)
	y$Sex.db=as.character(y$Sex.db)
	y$Sex.db[is.na(y$Sex.db)]=""
	# Integrate the Sex data from FMPro database with that from spreadsheet
	unknownSexInSpreadsheet=!y$Sex%in%c("M","F")
	y$Sex=as.character(y$Sex)
	y$Sex[unknownSexInSpreadsheet]=y$Sex.db[unknownSexInSpreadsheet]
	y$Sex=factor(y$Sex)
	if(is.factor(x$ID) && length(x$ID)==nlevels(x$ID)){
		rownames(y)=as.character(y$ID)
		y=y[x$ID,]
	}
	y[,setdiff(colnames(y),"Sex.db")]
}

deformationField=function(warplistfile,xyzs,xs=seq(0,168.78,len=16),ys=seq(0,168.78,len=16),zs=c(40,77),
	gregxform=file.path(cmtk.bindir(check=TRUE),"gregxform"),direction=c("inverse","forward"), request=c("warp","affine"),opts=NULL){
	
	direction=match.arg(direction)
	#nb inverse implies from sample to ref

	tmpfile=tempfile()
	tmpfile2=tempfile()
	if(missing(xyzs))
		xyzs=expand.grid(X=xs,Y=ys,Z=zs)
	write.table(xyzs,col.names=F,row.names=F,file=tmpfile)

	if(direction=="forward") gregxform=paste(gregxform,"--forward")
	if(!is.null(opts)) gregxform=paste(gregxform,opts)
	
	l=list(pre=xyzs)
	if("affine"%in%request){
		system(paste(gregxform,"--affine",warplistfile,"< ",tmpfile,">",tmpfile2),
			intern=FALSE,ignore.stderr=FALSE)
		l$aff=read.table(tmpfile2,col.names=c("X","Y","Z"))
	}
	
	if("warp"%in%request){
		system(paste(gregxform,warplistfile,"< ",tmpfile,">",tmpfile2),
			intern=FALSE,ignore.stderr=TRUE)
		l$warp=read.table(tmpfile2,na.strings="ERR",col.names=c("X","Y","Z"))
	}
	unlink(c(tmpfile,tmpfile2))
	return(l)
}

plotFieldComparison<-function(l,compare=c("aff","warp"),
	scl=c(1,1,1),ylim=c(168.78,0),xlim=c(0,168.78),withPlot=TRUE,axes=F,
	withPoints=FALSE,withArrows=TRUE,withGrid=FALSE,
	zfunc='all',margins=FALSE,pcol='blue',gcol='green',acol='red'){

	if(is.function(zfunc)){
		pre=subset(l[[compare[1]]],l$pre$Z==zfunc(l$pre$Z))
		post=subset(l[[compare[2]]],l$pre$Z==zfunc(l$pre$Z))		
	} else{
		pre=subset(l[[compare[1]]])
		post=subset(l[[compare[2]]])		
	}

	#cat(length(scl))
	if(length(scl)==1) scl=rep(scl,2)
	if(length(scl)==2) scl=c(scl,1)
	ylim=ylim*scl[2]
	pre=as.data.frame(t(t(pre)*scl))
	post=as.data.frame(t(t(post)*scl))

	if(!withPlot){
		pre$Y=ylim[1]-pre$Y
		post$Y=ylim[1]-post$Y
	}

	# plot points at initial pos by default
	if(is.character(withPoints) && withPoints=="post") ps=post else ps=pre
	if(is.character(withPoints)) withPoints=TRUE

	if(!margins) par(mar=rep(0,4))
	if(withPlot && withPoints) plot(ps$X,ps$Y,pch=20,cex=.75,
		col=pcol,ylim=ylim,xlim=xlim,xaxs='i',yaxs='i',axes=axes)
	if(withPlot && !withPoints) plot(ps$X,ps$Y,type='n',
		ylim=ylim,xlim=xlim,xaxs='i',yaxs='i',axes=axes)

	if(!withPlot && withPoints) points(ps$X,ps$Y,pch=20,cex=.75,col=pcol)
	if(withArrows) arrows(pre$X,pre$Y,post$X,post$Y,col=acol,len=0.05)

	box()
	#abline(h=ylim[1]);abline(h=ylim[0])
	#abline(v=xlim[1]);abline(v=xlim[0])

	if(is.character(withGrid)){
		if(withGrid=="pre") pgs=pre
		if(withGrid=="post") pgs=post
		withGrid=TRUE
	} else pgs=post

	if(withGrid==TRUE){
		z=by(pgs,l$pre$Y[l$pre$Z==zfunc(l$pre$Z)],function(x) lines(x$X,x$Y,col=gcol))
		z=by(pgs,l$pre$X[l$pre$Z==zfunc(l$pre$Z)],function(x) lines(x$X,x$Y,col=gcol))
	}
}

# Function to calculate the "average" warping across the whole brain
# This will use Chris' RegistrationSurfaceInside-labels.am

AverageJacobianAcrossMaskedBrain<-function(brain,mask,grid){
	# TODO - INCOMPLETE
	
	if(is.character(mask)) mask=Read3DDensityFromAmiraLattice(mask)
	pmin(attr(lhdxyz[[1]],"BoundingBox"),attr(i$S,"BoundingBox"))[1:3]
	#if(missing(grid)) 
}

