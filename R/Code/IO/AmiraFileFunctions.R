# AmiraFileFunctions.R

# 2005-02-03
# Functions to parse AmiraMesh 3D format - the native
# ouput of the skeletonize plugin and to read and write the density
# data in Amira file formats.
# At the moment depends on SWCFunctions.R since the amiramesh 
# data format can be easily converted to SWC.  However there
# is some redundancy in this approach

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

# source(file.path(CodeDir,"AmiraFileFunctions.R"))

require(tools) # for md5sum

ReadAmiramesh<-function(filename,DataSectionsToRead=NULL,Verbose=FALSE,AttachFullHeader=FALSE,Simplify=TRUE,endian=NULL){
  .Deprecated('nat::read.amiramesh')
  nat::read.amiramesh(file=filename,sections=DataSectionsToRead,
    header=AttachFullHeader,simplify=Simplify,Verbose=Verbose,endian=endian)
}

ReadAmiramesh.Header<-function(con,Verbose=TRUE,CloseConnection=NULL){
  .Deprecated('nat::read.amiramesh.header')
  
  if( is.logical(CloseConnection) && xor(is.character(con), CloseConnection) )
    stop("If con is a character vector, CloseConnection must now be TRUE.",
        " If not CloseConnection must now be FALSE")
  nat::read.amiramesh.header(file=con,Verbose=Verbose)
}

WriteNeuronToAM<-function(ANeuron,AMFile=NULL,
	suffix="am",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE,
  WriteRadius=TRUE){
	
	file=AMFile; dir=NULL
	if(!is.null(AMFile) && isTRUE(file.info(AMFile)$isdir)){
		# we've been given a directory
		# we want to write a file into this directory with same name as original
		file=NULL
		dir=AMFile
	}
	.Deprecated('nat::write.neuron')
	nat::write.neuron(ANeuron, file=file, dir=dir, format='hxlineset', 
		Force=Force, MakeDir=MakeDir, WriteAllSubTrees=WriteAllSubTrees,
		ScaleSubTreeNumsTo1=ScaleSubTreeNumsTo1, WriteRadius=WriteRadius)
}

WritePointsToAM<-function(d,AMFile=NULL,suffix="am",Force=F,MakeDir=T){
	# write out a set of points in the  basic AmiraMesh format

	if(is.neuron(d)){
		if(is.null(AMFile)) AMFile=paste(sub("(.*)\\.[^.]*$","\\1",d$InputFileName),sep=".",suffix)
		d=d$d
	} else if (is.null(AMFile)){
		warning("No file name specified in WritePointsToAM")
		return()
	}
	
	if(!Force && file.exists(AMFile) ){
		warning(paste(AMFile,"already exists; use Force=T to overwrite"))
		return()
	}
	if(!file.exists(dirname(AMFile))){
		# either bail
		if(!MakeDir){
			warning(paste(dirname(AMFile),"does not exist; use MakeDir=T to overwrite"))
			return()
		} else {
			# or try to make a directory
			if(!dir.create(dirname(AMFile))){
				warning(paste("Unable to create",dirname(AMFile)))
			}
		}
	}
	if(!file.create(AMFile)){
		warning(paste("Unable to write to file",AMFile))
		return()
	}
  
  # assign column names 
  if(ncol(d)==3){
    colnames(d)=c("X","Y","Z")
  } else if(is.null(colnames(d))){
		stop("Unable to identify X,Y,Z coordinates")
	} else {
      d=d[,c("X","Y","Z")]
	}
	
	cat("Writing to",AMFile,"\n")
	# Write the header
	cat("# AmiraMesh ASCII 1.0\n",file=AMFile)
	fc=file(AMFile,open="at") # ie append, text mode

	cat("# Created by WritePointsToAM - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	nVertices=nrow(d)
	cat("define Markers",nVertices,"Parameters {\nContentType \"LandmarkSet\",nSets 1\n}\n",file=fc)
# 	cat("Parameters {\n",file=fc)
# 	cat("    ContentType \"HxLineSet\"\n",file=fc)
# 	cat("}\n\n",file=fc)

	cat("Markers { float[3] Coordinates } = @1\n",file=fc)
	cat("\n",file=fc)
	
	# Write the 3D coords
	cat("@1 # ",nVertices,"xyz coordinates\n",file=fc)
	#write(t(ANeuron$d[,c("X","Y","Z")]),ncolumns=3,file=fc)
	write.table(d[,c("X","Y","Z")],col.names=F,row.names=F,file=fc)
	
	close(fc)
}

WriteNeuronToAM3D<-function(ANeuron,AMFile=NULL,
	suffix="am3",Force=F,MakeDir=T,WriteAllSubTrees=TRUE,ScaleSubTreeNumsTo1=TRUE,sep=NULL){
	# write out a neuron in the specialised skeletonize AM3D format 
	# (as opposed to the basic AmiraMesh format which is the native format
	# of amira for linesets)
	# WriteAllSubTrees will write out all the stores subtrees in a neuron 
	# which has multiple subtrees (which is often true of ill-formed 
	# skeletonize neurons).  It will also add a data field that can be used
	# to visualised different subtrees eg by colouring
	
	file=AMFile; dir=NULL
	if(!is.null(AMFile) && isTRUE(file.info(AMFile)$isdir)){
		# we've been given a directory
		# we want to write a file into this directory with same name as original
		file=NULL
		dir=AMFile
	}
	.Deprecated('nat::write.neuron')
	nat::write.neuron(ANeuron, file=file, dir=dir, format='hxskel', Force=Force,
		MakeDir=MakeDir,WriteAllSubTrees=WriteAllSubTrees,
		ScaleSubTreeNumsTo1=ScaleSubTreeNumsTo1,sep=sep)
}

ReadMBLHFromAMSurf<-function(AMSurfFile,Components=c("MB","LH"),WithConvexHull=F,...){
	# Function to read in a neuron from an AmiraMesh 3D ascii file
	# this is the file format of the Skeletonize plugin
	# I have adapted this from ReadNeuronFromAsc
	# However not all the options are functional now since they are not
	# really appropriate to the new format
	# In particular there is no contour data associated with this file
	# format.
	# nb if WithConvHull ==T then use convex hull to reduce number of data 
	# points in LH and MB contours 

	if(length(intersect(Components,c("MB","LH")))>0 ){	
		# We have some data to get
		cat("We have some data to get\n")
		Data<-ParseAMSurfToContourList(AMSurfFile,...)
		return(Data) 
	}
	return(-1)
}


ParseAMSurfToContourList<-function(filename,RegionNames="ALL",RegionChoice="Inner",Verbose=FALSE,FallbackRegionCol="grey"){
	.Deprecated('nat::read.hxsurf')
	nat::read.hxsurf(filename,RegionNames=if(RegionNames=="ALL") NULL else RegionNames,
		RegionChoice=RegionChoice,Verbose=Verbose,FallbackRegionCol=FallbackRegionCol)
}

WriteHxSurface=function(filename,Vertices,Indices=NULL,
		material=sub("^([^.]+)\\..*","\\1",basename(filename))){
		
	cat("# HyperSurface ASCII\nParameters {\n",file=filename)
	fc=file(filename,open="at") # ie append, text mode
	#cat("\tMaterials{\n\t\tInterior {\n\t\t\tid 0\n\t\t}\n",file=fc)
	cat("\t{color 0.83562 0.78 0.06,\nName \"",sep="",material,"\"}",file=fc)
	cat("\t\t",sep="",material," {\n\t\t\tid 1\n\t\t}\n",file=fc)
	cat("\t}\n",file=fc)
	cat("\tBoundaryIds {\n\t\tId0 {\n\t\t\tId 0,Info \"undefined\",Color 0.6 0.6 0.6\n\t\t}\n\t\tname \"BoundaryConditions\"\n\t}\n}\n",file=fc)


	nVertices=nrow(Vertices)
	cat("Vertices",nVertices,"\n",file=fc)
	write.table(Vertices,col.names=F,row.names=F,file=fc)
	cat("Patches 1\n{\tInnerRegion Interior \n OuterRegion", material,"\n BoundaryID 0\n BranchingPoints 0 \n\n Triangles",round(nVertices/3),"\n",file=fc)
	#cat("Markers { float[3] Coordinates } = @1\n@1\n",file=filename,append=T)
	Indices=matrix(1:nVertices,nrow=3)
	
	cat(apply(t(Indices),1,paste,collapse=" "),sep="\n",file=fc)
	cat("}\n",file=fc)
}

Write3DDensityToAmiraRectilinear<-function(filename,d){
		# Produces a Rectilinear format file -
		# that is one with an arbitrary x,y,z grid
		cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
		fc=file(filename,open="at") # ie append, text mode
		lattice=apply(d$eval.points,2,length)
		cat("define Lattice",lattice,"\n",file=fc)
		cat("define Coordinates",sum(lattice),"\n\n",file=fc)
		cat("Parameters {CoordType \"rectilinear\"}\n\n",file=fc)
		cat("Lattice { float ScalarField } = @1\n",file=fc)
		cat("Coordinates { float xyz } = @2\n\n",file=fc)
		cat("@1\n",file=fc)
		
		cat(as.vector(d$estimate),file=fc)
		cat("@2\n",file=fc)
		cat(d$eval.points[,1],"\n",file=fc)
		cat(d$eval.points[,2],"\n",file=fc)
		cat(d$eval.points[,3],"\n",file=fc)
}

Write3DDensityToAmiraLattice<-function(filename,dens,ftype=c("binary","text","hxzip"),
	dtype=c("float","byte", "short", "ushort", "int", "double"),WriteNrrdHeader=FALSE,endian=c('big','little')){
	# Produces a lattice format file -
	# that is one with a regular x,y,z grid
	# Can also write a detached Nrrd header that points to the AmiraMesh
	# data to allow it to be opened by a nrrd reader
	ftype=match.arg(ftype)
	endian=match.arg(endian)
	if(ftype=='text') cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
	else if(endian=='little') cat("# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n",file=filename)
	else cat("# AmiraMesh 3D BINARY 2.0\n\n",file=filename)
	
	fc=file(filename,open="at") # ie append, text mode
	cat("# Created by Write3DDensityToAmiraLattice - ",format(Sys.time(),usetz=T),"\n\n",file=fc)	

	if(!is.list(dens)) d=dens else d=dens$estimate
	# Find data type and size for amira
	dtype=match.arg(dtype)	
	dtypesize<-c(4,1,2,2,4,8)[which(dtype==c("float","byte", "short","ushort", "int", "double"))]
	# Set the data mode which will be used in the as.vector call at the
	# moment that the binary data is written out.
	if(dtype%in%c("byte","short","ushort","int")) dmode="integer"
	if(dtype%in%c("float","double")) dmode="numeric"
	
	
	#lattice=apply(d$eval.points,2,length)
	lattice=dim(d)
	cat("define Lattice",lattice,"\n",file=fc)

	cat("Parameters { CoordType \"uniform\",\n",file=fc)
	# note Amira's definition for the bounding box:
	# the range of the voxel centres.
	# So eval.points should correspond to the CENTRE of the
	# voxels at which the density is evaluated
	cat("\t# BoundingBox is xmin xmax ymin ymax zmin zmax\n",file=fc)
	BoundingBox=NULL
	if(!is.null(attr(dens,"BoundingBox"))){
		BoundingBox=attr(dens,"BoundingBox")
	} else if(is.list(d) && !is.null(d$eval.points)){
		BoundingBox=as.vector(apply(d$eval.points,2,range))
	}
	if(!is.null(BoundingBox)) cat("\t BoundingBox",BoundingBox,"\n",file=fc)
	cat("}\n\n",file=fc)
	
	if(ftype=="hxzip"){
		raw_data=writeBin(as.vector(d,mode=dmode),raw(),size=dtypesize,endian=endian)
		zlibdata=nat:::write.zlib(raw_data)
		cat("Lattice { ",dtype," ScalarField } = @1(HxZip,",length(zlibdata),")\n\n",sep="",file=fc)
	} else cat("Lattice {",dtype,"ScalarField } = @1\n\n",file=fc)

	cat("@1\n",file=fc)
	
	#cat(str(as.vector(d)))
	close(fc)

	# Write a Nrrd header to accompany the amira file if desired
	# see http://teem.sourceforge.net/nrrd/
	if(WriteNrrdHeader) {
		if(ftype=="hxzip") stop("Nrrd cannot cope with Amira's HxZip encoding (which is subtly different from gzip)")
		nrrdFilename=paste(filename,sep=".","nhdr")
		cat("NRRD0004\n",file=nrrdFilename)
		fc=file(nrrdFilename,open="at") # ie append, text mode
		nrrdType=ifelse(dtype=="byte","uint8",dtype)
		
		cat("encoding:", ifelse(ftype=="text","text","raw"),"\n",file=fc)
		cat("type: ",nrrdType,"\n",sep="",file=fc)
		cat("endian: ",endian,"\n",sep="",file=fc)
		# Important - this sets the offset in the amiramesh file from which
		# to start reading data
		cat("byte skip:",file.info(filename)$size,"\n",file=fc)
		cat("dimension: ",length(lattice),"\n",sep="",file=fc)
		cat("sizes:",lattice,"\n",file=fc)
		voxdims=voxdim.gjdens(dens)
		if(!is.null(voxdims)) cat("spacings:",voxdims,"\n",file=fc)
		if(!is.null(BoundingBox)){
			cat("axis mins:",matrix(BoundingBox,nrow=2)[1,],"\n",file=fc)
			cat("axis maxs:",matrix(BoundingBox,nrow=2)[2,],"\n",file=fc)
		}
		cat("data file: ",basename(filename),"\n",sep="",file=fc)
		cat("\n",file=fc)
		close(fc)
	}

	if(ftype=='text'){
		write(as.vector(d,mode=dmode),ncol=1,file=filename,append=TRUE)
	} else {
		fc=file(filename,open="ab") # ie append, bin mode
		if(ftype=="hxzip")
			writeBin(zlibdata,fc,size=1,endian=endian)
		else
			writeBin(as.vector(d,mode=dmode),fc,size=dtypesize,endian=endian)
		close(fc)
	}
}

Read3DDensityFromAmiraLattice<-function(filename,Verbose=FALSE){

	fc=file(filename,'rb')
	headerLines<-readLines(fc,1)
	if(length(grep("amiramesh",headerLines,ignore.case=T))!=1)
		stop("This is not an amiramesh file")
	while ( ( nextLine<-readLines(fc,1)) !="@1") {headerLines<-c(headerLines,nextLine)}
	
	endian='big'
	binary=FALSE
	# Figure out if the file is in binary format or not
	if(length(grep("LITTLE.ENDIAN",headerLines[1],ignore.case=TRUE))>0){
#		any(grep("^\\s*#\\s+amiramesh(\\s+3d)?\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
		binary = TRUE
		endian = 'little'
	} else if(any(grep("^\\s*#\\s+amiramesh(\\s+3d)?\\s+binary", headerLines[1],ignore.case=TRUE,perl=TRUE))){
		binary = TRUE
	}

	# Find the position of the header lines defining the lattice
	# this clearly assumes that the relevant data is in position
	# @1 in the file - which may or may not be the case.
	LatticeDefLine=grep("define\\s+lattice",headerLines,ignore.case=TRUE,perl=TRUE)
	LatticeTypeDefLine=grep("^Lattice.*}\\s*[=]{0,1}\\s*@1",headerLines,ignore.case=TRUE,perl=FALSE)
	#cat("LatticeTypeDefLine = ",LatticeTypeDefLine)
	LatticeBoundsLine=grep("^[^#]*BoundingBox",headerLines,ignore.case=TRUE,perl=TRUE)
	if(!any(LatticeDefLine)){
		warning(paste("No lattice definition line in file",filename))
		close(fc); return(-1)
	}
	
	if(!any(LatticeTypeDefLine)) {
		warning(paste("No lattice type definition line in file",filename))
		close(fc); return(-1)
	}

	strtrim<-function(str){
			str<-sub('\\s+$', '', str, perl = TRUE) ## Perl-style white space
			sub('^\\s+', '', str, perl = TRUE) ## Perl-style white space
	}
	
	# fetch the lattice dimensions
	latticeDims=sub(".*lattice([0-9 \\t]+).*","\\1",headerLines[LatticeDefLine],ignore.case=T)
	latticeDims=as.numeric(unlist(strsplit(strtrim(latticeDims),"\\s+")))
	dataLength=prod(latticeDims)
	
	# and the bounds
	latticeBounds=as.numeric(unlist(strsplit(strtrim(headerLines[LatticeBoundsLine]),"(\\s|,)"))[-1])		
	
	
	# and finally the data type
	ldl=headerLines[LatticeTypeDefLine]
	# dataTypeName=sub(".*\\{\\s*(\\w+)\\s+.?*\\s*\\}.*","\\1",ldl,ignore.case=T,perl=T)
	dataTypeName=sub(".*\\{\\s*(\\w+)\\s+.*?\\s*\\}.*","\\1",ldl,ignore.case=T,perl=T)
	# check if the data is encoded in some way
	#"Lattice { byte Labels } @1(HxByteRLE,391697)"
	postAt=strsplit(ldl,"@1")[[1]][2]
	dataEncoding=""
	if(isTRUE(grep("\\(Hx.*?\\)",postAt)==1)){
		postAt=sub(".*?\\(([^)]+)\\).*?","\\1",postAt)
		dataEncoding=toupper(strsplit(postAt,",")[[1]][1])
		compressedLength=as.integer(strsplit(postAt,",")[[1]][2])
	}
	
#     Amira docs: The primitive data types must be
# 		one of byte, short, int, float, double, or complex.  Vectors of
# 		primitive data types are allowed, aggregate structs are not, however.
	
# 		primType Returns the primitive data type of the field, i.e., the way
# 		how the values are represented internally.  A number with the following
# 		meaning is returned: 0 = bytes, 1 = 16-bit signed integers, 2 = 32-bit
# 		signed integers, 3 = 32-bit floating point values, 4 = 64-bit floating
# 		point values, 7 = 16-bit unsigned integers.
	
	# now read the data
	# note that  bytes are assumed to be unsigned
	# shorts could be either but will assume signed - don't know how
	# amiramesh specifies either way
	dataTypes=data.frame(name=I(c("byte", "ushort", "short", "int", "float", "double", "complex")),
			size=c(1,2,2,4,4,8,NA),what=I(c(rep("integer",4),rep("numeric",2),NA)),
			signed=rep(c(FALSE,TRUE),c(2,5)) )
	i=which(dataTypes$name==dataTypeName)
	if(!any(i==1:7)){
		close(fc)
		stop("Unrecognised data type")
	}
	
	if(Verbose) cat("dataLength =",dataLength,"dataType =",dataTypes$what[i],"size=",dataTypes$size[i],"\n")
	if(binary){
		if(dataEncoding=="HXBYTERLE"){
			d=readBin(fc,what=raw(0),n=compressedLength,size=1)
			d=nat:::decode.rle(d,dataLength)
			d=as.integer(d)
		} else if(dataEncoding == "HXZIP"){
		  uncompressed=nat:::read.zlib(fc, compressedLength=compressedLength)
		  d=readBin(uncompressed, n=dataLength, what=dataTypes$what[i],
		            size=dataTypes$size[i], signed=dataTypes$signed[i],
                endian=endian)
		} else if(dataEncoding==""){
			d=readBin(fc,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
				signed=dataTypes$signed[i],endian=endian)
		} else {
			stop("Unimplemented data encoding",dataEncoding,"in file",filename,"\n")
		}
		close(fc)
	} else {
		# this clearly assumes that the relevant data is in position
		# @1 in the file - which may or may not be the case.
		# cat("dataTypes$what[i]=",dataTypes$what[i],"\n")
		if(dataTypes$what[i]=='integer') whatVal=integer(0) else whatVal=double(0)
		d=scan(fc,what=whatVal,nmax=dataLength)
		close(fc)
	}
	dim(d)<-latticeDims
	if(length(latticeBounds)>0){
		attr(d,"BoundingBox")<-latticeBounds
		attr(d,"x")<-seq(latticeBounds[1],latticeBounds[2],len=latticeDims[1])
		attr(d,"y")<-seq(latticeBounds[3],latticeBounds[4],len=latticeDims[2])
		attr(d,"z")<-seq(latticeBounds[5],latticeBounds[6],len=latticeDims[3])
	} else {
		# No Bounding Box available
# 			attr(d,"BoundingBox")<-NULL
	}
	return(d)
}

ReadAmiraLandmarks<-function(filename,Verbose=FALSE,CoordinatesOnly=TRUE){
  .Deprecated('nat::read.amiralandmarks')
  nat::read.amiralandmarks(filename, CoordinatesOnly=CoordinatesOnly, 
    Verbose=Verbose)
}

WriteAmiraLandmarks<-function(filename,d){
  .Deprecated('nat::write.amiralandmarks')
  nat::write.amiralandmarks(file=filename, x=d)
}

WriteAmiraColormap<-function(filename,rgba,A=1,minmax=c(0,255)){
	if(is.vector(rgba)) rgba=t(col2rgb(rgba)/255)
	
	if(ncol(rgba)==3) {
		rgba=cbind(rgba,A)
	}
	if(ncol(rgba)!=4) stop("Colormap must have 4 columns")
	if(nrow(rgba)!=256) warning("Colormap should have 256 levels to be editable with colormap editor")
	
	cmaprange=range(rgba)
	if(cmaprange[1]<0 || cmaprange[2]>1) stop ("Colormap values must be between 0 and 1")

	cat("# AmiraMesh ASCII 1.0\n\ndefine Lattice",nrow(rgba),
	"\n\nParameters {\n\tContentType \"Colormap\",\n\tMinMax",minmax,"\n}\n",file=filename)
	cat("Lattice { float[4] Data } = @1\n",file=filename,append=T)

	cat("@1\n",file=filename,append=T)
	write.table(rgba,col.names=F,row.names=F,file=filename,append=TRUE)
	cat("\n",file=filename,append=T)
}

WriteGenericAmiramesh<-function(filename,d,ContentType){
	dataDef=attr(d,"dataDef")
	if(is.null(dataDef)) stop("Cannot write without data definition")
	
	cat("# AmiraMesh ASCII 1.0\n\n",file=filename)
	# definitions of lengths
	dimdefs=dataDef$Dims[unique(names(dataDef$Dims))]
	cat(paste("define",names(dimdefs),dimdefs,collapse="\n"),file=filename,append=TRUE)

	# Parameters - for now just handle content type
	if(missing(ContentType)) ContentType=attr(d,'Parameters')$ContentType
	if(!is.null(ContentType)){
		cat("\nParameters {\n\tContentType \"",, "\"\n}\n\n",file=filename,append=T)
	}

	with(dataDef,
		cat(paste(names(Dims)," { ",Type," ",DataName," } @",sep="",seq(nrow(dataDef)),collapse="\n"),
			file=filename,append=TRUE))
	
	cat("\n\n",file=filename,append=TRUE)
	if(!is.list(d)) d=list(d)
	for(i in seq(length(d))){
		if(length(d[[i]])==0) next 
		cat("@",i,"\n",sep="",file=filename,append=TRUE)
		write.table(d[[i]],row.names=FALSE,col.names=FALSE,sep="\t",file=filename,append=TRUE)
		cat("\n",file=filename,append=TRUE)
	}		
}


#' Read neuron in Amira's native lineset format
#' @param amfile Path to the amiramesh file
#' @param defaultDiameter If diameter information, missing use this default
#' @return A neuron object
#' @author jefferis
#' @export
#' @seealso \code{\link{read.neuron},\link{ReadNeuronFromAM3D}
ReadNeuronFromAM<-function(amfile,defaultDiameter=NA_real_){
  .Deprecated('nat::read.neuron.hxlineset')
  nat:::read.neuron.hxlineset(amfile, defaultDiameter=defaultDiameter)
}

