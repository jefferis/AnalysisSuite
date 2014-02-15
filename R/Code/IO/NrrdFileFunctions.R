ReadNrrdHeader<-function(filename,Verbose=TRUE,CloseConnection=TRUE){
  .Deprecated("nat::read.nrrd.header")
  nat::read.nrrd.header(filename, Verbose=Verbose)
}

NrrdDataFiles<-function(nhdr,ReturnAbsPath=TRUE){
	if(!is.list(nhdr)){
		# we need to read in the nrrd header
		if(length(nhdr)>1) return(sapply(nhdr,NrrdDataFiles))
		if(!is.nrrd(nhdr)) stop("This is not a nrrd file")
		h=ReadNrrdHeader(nhdr)
	} else h=nhdr
	if(is.null(h$datafile)){
		# straight nrrd without detached header
		dfs=attr(h,'path')
	} else if(length(h$datafile)>1){
		# list of files
		dfs=h$datafile[-1]
		firstlineparts=unlist(strsplit(h$datafile[1],"\\s+",perl=T))
		if(length(firstlineparts)>2) stop("invalid first line of LIST datafile specifier")
		if(length(firstlineparts)==2) attr(dfs,'subdim')=as.integer(firstlineparts[2])
	} else if(grepl("%",h$datafile)){
		# Format specifier TODO
		firstlineparts=unlist(strsplit(h$datafile[1],"\\s+",perl=T))
		if(length(firstlineparts)>5 || length(firstlineparts)<4)
			stop("invalid sprintf style datafile specifier")
		rangespec=as.integer(firstlineparts[2:4])
		dfs=sprintf(firstlineparts[1],
			seq(from=rangespec[1],to=rangespec[2],by=rangespec[3]))
		if(length(firstlineparts)==5) attr(dfs,'subdim')=as.integer(firstlineparts[5])
	} else dfs=h$datafile
	
	if(ReturnAbsPath){
		# check if paths begin with /
		relpaths=substring(dfs,1,1)!="/"
		if(any(relpaths)){
			nhdrpath=attr(h,"path")
			if(is.null(nhdrpath) && ReturnAbsPath)
				stop("Unable to identify nrrd header file location to return absolute path to data files")
			dfs[relpaths]=file.path(dirname(nhdrpath),dfs[relpaths])
		}
	}
	dfs
}

.standardNrrdType<-function(type){
  nat:::.standardNrrdType(type)
}

#' Read nrrd file into 3d array in memory
#' 
#' @details ReadByteAsRaw=unsigned (the default) only reads unsigned byte data
#'   as a raw array. This saves quite a bit of space and still allows data to be
#'   used for logical indexing.
#' @param filename Path to a 3D nrrd
#' @param Verbose Show data length (default FALSE)
#' @param ReadData When FALSE just return attributes (e.g. voxel size)
#' @param AttachFullHeader Include the full nrrd header as an attribute of the
#'   returned object (default FALSE)
#' @param ReadByteAsRaw Read 8 bit data as an R raw object rather than integer
#' @param origin Add a user specified origin (x,y,z) to the returned object
#' @return a 3D data array with attributes compatible with gjdens objects
#' @export
Read3DDensityFromNrrd<-function(filename,Verbose=FALSE,ReadData=TRUE,AttachFullHeader=!ReadData,
	ReadByteAsRaw=c("unsigned","all","none"),origin){
	ReadByteAsRaw=match.arg(ReadByteAsRaw)
	fc=file(filename,'rb')
	h=ReadNrrdHeader(fc,CloseConnection=FALSE)
	# store the path because ReadNrrdHeader couldn't do it 
	# TODO more elegant way of dealing with paths when connection sent to ReadNrrdHeader
	attr(h,'path')=filename

	# now read the data
	dataTypes=data.frame(name=I(c("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64",
	"float", "double", "block")),
			size=c(1,1,2,2,4,4,8,8,4,8,NA),what=I(c(rep("integer",8),rep("numeric",2),"raw")),
			signed=c(rep(c(T,F),4),rep(T,3)))
	if(ReadByteAsRaw=="all") dataTypes$what[1:2]='raw'
	else if(ReadByteAsRaw=="unsigned") dataTypes$what[2]='raw'
	
	i=which(dataTypes$name==.standardNrrdType(h$type))
	if(length(i)!=1){
		close(fc)
		stop("Unrecognised data type")
	}
	if(!is.null(h$datafile)){
		# detached nrrd
		if(!inherits(filename,"connection"))
			attr(h,'path')=filename
		datafiles=NrrdDataFiles(h)
		# TODO - handle more than one datafile!
		if(length(datafiles)!=1) stop("Can currently only handle exactly one datafile")
		close(fc)
		fc=file(datafiles,open='rb')
		filename=datafiles
	}
	dataLength=prod(h$sizes)
	endian=ifelse(is.null(h$endian),.Platform$endian,h$endian)
	if(Verbose) cat("dataLength =",dataLength,"dataType =",dataTypes$what[i],"size=",dataTypes$size[i],"\n")
	enc=tolower(h$encoding)
	if(ReadData){
		if(enc=="raw"){
			d=readBin(fc,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
				signed=dataTypes$signed[i],endian=endian)
			close(fc)
		} else if(enc%in%c("gz","gzip")) {
			# unfortunately gzcon seems to reset the connection 
			# rather than starting to read from the current location
			headerlength=seek(fc)
			close(fc)
			tf=tempfile()
			system(paste('tail -c +',sep="",headerlength+1," ",filename," > ",tf))
			gzf=gzfile(tf,'rb')
			d=readBin(gzf,what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
				signed=dataTypes$signed[i],endian=endian)
			close(gzf)
			unlink(tf)
		} else if(enc%in%c("ascii","txt","text")){
			if(dataTypes$what[i]=='integer') whatVal=integer(0) else whatVal=double(0)
			d=scan(fc,what=whatVal,nmax=dataLength,quiet=TRUE)
			close(fc)
		} else {
			stop("nrrd encoding ",enc," is not implemented")
		}
		dim(d)<-h$sizes
	} else {
		# don't read the data, we just wanted the (full) header information
		if(dataTypes$what[i]=='integer') d=integer(0) else d=double(0)
		attr(d,'datablock')=list(what=dataTypes$what[i],n=dataLength,size=dataTypes$size[i],
			signed=dataTypes$signed[i],endian=endian)
		attr(d,'datablock')$datastartpos=seek(fc)
		close(fc)
	}
	if(AttachFullHeader) attr(d,"header")=h
	voxdims<-NrrdVoxDims(h,ReturnAbsoluteDims = FALSE)
	if(any(is.na(voxdims))){
		# missing pixel size info, so just return
		return(d)
	}
	latticeBoundingBox=rbind(c(0,0,0),(h$sizes-1)*voxdims)
	if(!missing(origin)){
		latticeBoundingBox=t(origin+t(latticeBoundingBox))
	} else if('space origin'%in%names(h)){
		# FIXME should space origin interpretation depend on node vs cell?
	  # (and we are assuming that we will be working as node in R)
		latticeBoundingBox=t(h[['space origin']]+t(latticeBoundingBox))
	}
	attr(d,"BoundingBox")<-as.vector(latticeBoundingBox)
	attr(d,"x")<-seq(latticeBoundingBox[1],latticeBoundingBox[2],len=h$sizes[1])
	attr(d,"y")<-seq(latticeBoundingBox[3],latticeBoundingBox[4],len=h$sizes[2])
	attr(d,"z")<-seq(latticeBoundingBox[5],latticeBoundingBox[6],len=h$sizes[3])
	return(d)
}

WriteNrrdHeaderForAmirameshFile<-function(amfile,outfile=paste(amfile,sep=".","nhdr")){
	h=ReadAmiramesh.Header(amfile)
	hd=h$dataDef
	if(nrow(hd)==0) return(NULL)
	if(nrow(hd)>1) warning("Can only use first data block of Amira File")
	hd=hd[1,]
	if(hd$HxType!="raw") stop("Unable to make a nrrd header for compressed Amiramesh files")
	
	nrrdEncodings=structure(c("raw",NA,NA),names=c("byte","HxByteRLE","HxZip"))
	nrrdDataTypes=structure(c("uint8","uint16","int16","int32","float","double",NA),
		names=c("byte", "ushort", "short", "int", "float", "double", "complex"))
	nrrdDataType=nrrdDataTypes[hd$SimpleType]
	if(is.na(nrrdDataType)) stop("Unable to write nrrd header for data type: ",hd$SimpleType)
		
	cat("NRRD0004\n",file=outfile)
	cat("encoding: raw\ntype: ",nrrdDataType,"\n",sep="",append=TRUE,file=outfile)
	dims=unlist(hd$Dims)
	cat("dimension: ",length(dims),"\nsizes: ",paste(dims,collapse=" "),"\n",sep="",append=TRUE,file=outfile)
	invisible(outfile)
}

Write3DDensityToNrrd<-function(filename,dens,enc=c("raw","text","gzip"),
	dtype=c("float","byte", "short", "ushort", "int", "double"),endian=c('big','little')){
	# Produces a lattice format file -
	# that is one with a regular x,y,z grid
	# Can also write a detached Nrrd header that points to the AmiraMesh
	# data to allow it to be opened by a nrrd reader
	enc=match.arg(enc)
	endian=match.arg(endian)
	dtype=match.arg(dtype)

	nrrdDataTypes=structure(c("uint8","uint16","int16","int32","float","double"),
		names=c("byte", "ushort", "short", "int", "float", "double"))
	
	nrrdDataType=nrrdDataTypes[dtype]
	if(is.na(nrrdDataType)) stop("Unable to write nrrd file for data type: ",dtype)
	cat("NRRD0004\n",file=filename)
	cat("encoding: ",enc,"\ntype: ",nrrdDataType,"\n",sep="",append=TRUE,file=filename)
	cat("dimension: ",length(dim(dens)),"\nsizes: ",paste(dim(dens),collapse=" "),"\n",sep="",append=TRUE,file=filename)
	voxdims=voxdim.gjdens(dens)
	if(!is.null(voxdims)) cat("spacings:",voxdims,"\n",file=filename,append=TRUE)

	if(!is.list(dens)) d=dens else d=dens$estimate
	# Find data type and size for amira
	dtype=match.arg(dtype)	
	dtypesize<-c(4,1,2,2,4,8)[which(dtype==c("float","byte", "short","ushort", "int", "double"))]
	# Set the data mode which will be used in the as.vector call at the
	# moment that the binary data is written out.
	if(dtype%in%c("byte","short","ushort","int")) dmode="integer"
	if(dtype%in%c("float","double")) dmode="numeric"
	# record byte ordering if necessary
	if(enc!='text' && dtypesize>1) cat("endian: ",endian,"\n",sep="",file=filename,append=TRUE)
	# Single blank line terminates header
	cat("\n",file=filename,append=TRUE)
	
	if(enc=='text'){
		write(as.vector(d,mode=dmode),ncol=1,file=filename,append=TRUE)
	} else {
		if(enc=="gzip") fc=gzfile(filename,"ab")
		else fc=file(filename,open="ab") # ie append, bin mode
		writeBin(as.vector(d,mode=dmode),fc,size=dtypesize,endian=endian)
		close(fc)
	}
}

Read3DDensityFromHanchuanRaw<-function(filename){
	fc=file(filename,'rb')
	on.exit(close(fc))
	hanchuanMagic="raw_image_stack_by_hpeng"
	magic=readChar(fc,nchars=nchar(hanchuanMagic))
	if(magic != hanchuanMagic){
		stop("This does not appear to be a Hanchuan RAW file.")
	}
	endian=readChar(fc, nchars=1)
	if(endian=="L"){
		endian="little"
	}else if(endian=="B"){
		endian="big"
	}else{
		stop("Unknown endianess.")
	}
	dataTypeSize=readBin(fc,what=integer(),n=1,size=2,endian=endian)
	if(dataTypeSize==1){
		what="integer"
	}else if(dataTypeSize==2){
		what="integer"
	}else if(dataTypeSize==4){
		what="numeric"
	}else{
		stop("Unknown datatype.")
	}
	dims=readBin(fc,what=integer(),n=4,size=4,endian=endian)
	dens=readBin(fc,what=what,n=prod(dims),size=dataTypeSize,endian=endian)
	# Keep only dimensions with more than 1 voxel.
	dim(dens)<-dims[dims>1]
	return(dens)
}

Write3DDensityToHanchuanRaw<-function(filename,dens,dtype=c("float","byte","ushort"),
	endian=c('little',"big"),WriteNrrdHeader=FALSE){
	endian=match.arg(endian)
	dtype=match.arg(dtype)
	hanchuanDataTypes=structure(c(1,2,4),
		names=c("byte", "ushort", "float"))
	hanchuanDataType=hanchuanDataTypes[dtype]
	hanchuanMagic="raw_image_stack_by_hpeng"
	
	if(dtype%in%c("byte","ushort")) dmode="integer"
	if(dtype=="float") dmode="numeric"
	cat(hanchuanMagic,file=filename)
	cat(toupper(substring(endian,1,1)),file=filename,append=T)

	con=file(filename,open='ab')
	
	# Write data type
	writeBin(as.integer(hanchuanDataType),con,2,endian=endian)

	# Write dimensions
	dimstowrite=dim(dens)
	# make sure that there are 4 dimensions (0 padding as required)
	dimstowrite=as.integer(c(dimstowrite,rep(1,4-length(dimstowrite))))
	writeBin(dimstowrite,con,4,endian=endian)
	
	# Write data
	writeBin(as.vector(dens,mode=dmode),con,size=hanchuanDataType,endian=endian)
	
	headerLength=nchar(hanchuanMagic)+1+2+4*4
	close(con)
	# Write a Nrrd header to accompany the amira file if desired
	# see http://teem.sourceforge.net/nrrd/
	if(WriteNrrdHeader) {
		nrrdFilename=paste(filename,sep=".","nhdr")
		cat("NRRD0004\n",file=nrrdFilename)
		fc=file(nrrdFilename,open="at") # ie append, text mode
		nrrdType=ifelse(dtype=="byte","uint8",dtype)
		
		cat("encoding: raw","\n",file=fc)
		cat("type: ",nrrdType,"\n",sep="",file=fc)
		cat("endian: ",endian,"\n",sep="",file=fc)
		# Important - this sets the offset in the amiramesh file from which
		# to start reading data
		cat("byte skip: ",headerLength,"\n",sep="",file=fc)
		cat("dimension: ",length(dim(dens)),"\n",sep="",file=fc)
		cat("sizes:",dim(dens),"\n",file=fc)
		voxdims=voxdim.gjdens(dens)
		if(!is.null(voxdims)) cat("spacings:",voxdims,"\n",file=fc)
		BoundingBox=getBoundingBox(dens)
		if(!is.null(BoundingBox)){
			cat("axis mins:",matrix(BoundingBox,nrow=2)[1,],"\n",file=fc)
			cat("axis maxs:",matrix(BoundingBox,nrow=2)[2,],"\n",file=fc)
		}
		cat("data file: ",basename(filename),"\n",sep="",file=fc)
		cat("\n",file=fc)
		close(fc)
	}
}

ConvertNrrdToAmira<-function(infile,outfile=sub("\\.nrrd$",".am",infile),dtype,
	TypeConversion=c("scale","cast"),...){
	TypeConversion=match.arg(TypeConversion)
	d=Read3DDensityFromNrrd(infile,AttachFullHeader=T)
	h=attr(d,"header")
	
	nrrdDataTypes=structure(names=c("uint8","uint16","int16","int32","float","double"),
		c("byte", "ushort", "short", "int", "float", "double"))
	nrrdType=.standardNrrdType(h$type)
	oldtype=nrrdDataTypes[nrrdType]

	if(missing(dtype)) {
		dtype=oldtype
	} else if(TypeConversion=="scale") {
		# we may need to rescale, just do this for int types
		saveattrs=attributes(d)
		if(oldtype=="ushort" && dtype=="byte") {
			d=as.integer(d/257)
			attributes(d)<-saveattr
		} else if(oldtype=="float") {
			r=range(d)
			if(dtype=="byte"){
				d=as.integer((d-r[1])/(r[2]/255))
				attributes(d)<-saveattr
			} else if(dtype=="ushort"){
				d=as.integer((d-r[1])/(r[2]/65535))
				attributes(d)<-saveattrs
			} else stop("Don't yet know how to convert ",oldtype," to ",dtype)
		}
		else stop("Don't yet know how to convert ",oldtype," to ",dtype)
	}
	Write3DDensityToAmiraLattice(outfile,d,dtype=dtype,...)
}

#' Read a 1D nrrd histogram (as generated by unu hist)
#'
#' Details: There will be one more break than the number of samples in the 1D
#'   nrrd, whereas the midpoints will match the levels of the nrrd. 
#' @param filename Path to 1D nrrd histogram
#' @param ... Additional arguments passed to Read3DDensityFromNrrd
#' @return object of class histogram
#' @export
#' @seealso \code{\link{Read3DDensityFromNrrd}}
#' @examples
ReadHistogramFromNrrd<-function(filename,...){
	d=Read3DDensityFromNrrd(filename,AttachFullHeader=TRUE,...)
	h=attr(d, "header")
	# clear the header attributes
	attributes(d)<-NULL
	if(is.na(pmatch("histo",h$content)) || h$dimension!=1) {
		warning ("This does not appear to be a 1d nrrd histogram")
		return(d)
	}
	halfwidth=(h$axismaxs-h$axismins)/(h$sizes-1)/2
	mids=seq(h$axismins,h$axismaxs,len=h$sizes)
	breaks=seq(from=h$axismins-halfwidth,to=h$axismax+halfwidth,len=h$sizes+1)
	density=d/sum(as.numeric(d))

	# return it as an R histogram	
	structure(list(
		breaks = breaks, 
		counts = d, 
	intensities = density, 
	density = density, 
	mids = mids,
	xname = h$content, 
	    equidist = TRUE), 
	.Names = c("breaks", "counts", "intensities", 
	"density", "mids", "xname", "equidist"), class = "histogram")
}

FixSpaceOrigin<-function(f,origin, Verbose=TRUE, KeepOriginalModificationTime = FALSE)
{
	if(length(origin)==6) origin=origin[c(1,3,5)]
	if(length(origin)!=3) stop("Supply either an origin or a bounding box")
	
	KeepBackup=TRUE # for now, insist on this
	# function to change the space origin field in a nrrd file
	inh=ReadNrrdHeader(f)
	originalOrigin=c(0,0,0)
	if("space origin"%in%names(inh)) originalOrigin=inh[["space origin"]]
	if(Verbose) cat("Old origin was:",originalOrigin,"; new origin is:",origin,"\n")
	newOrigin=origin
	newOriginLine=paste("space origin: (",paste(newOrigin,collapse=","),")",sep="")
	
	oht=attr(inh,"headertext")
	if("space origin"%in%names(inh)){
		# replace existing space origin
		oht=sub("space origin: .*",newOriginLine,oht)
	} else {
		# just append
		oht=c(oht,newOriginLine)
	}
	# add a blank line
	oht=c(oht,"")
	tmpheader=tempfile()
	tmpfile=tempfile()
	writeLines(oht,tmpheader)
	oldfile=f
	if(KeepBackup){
		oldfile=paste(f,sep="",".bak")
		if(!file.rename(f,oldfile)) stop("Unable to rename",f,"to",oldfile)
	}
	system(paste("unu data",shQuote(oldfile),"| cat",tmpheader,"- >",shQuote(f)))
	if(KeepOriginalModificationTime  && KeepBackup){
		cmd=paste("touch -am -r",shQuote(oldfile),shQuote(f))
		system(cmd)
	}
	unlink(c(tmpfile,tmpheader))	
}
#' Add or replace lines in the header of a nrrd file
#'
#' Input is a named vector of new fields. Use unnamed fields for comments. 
#' Quote field names with spaces (conventionally with back ticks).
#' Note that this function will error out for invalid field names.
#' See http://teem.sourceforge.net/nrrd/format.html for nrrd field details
#' @param infile Path to input file
#' @param outfile Path to output file
#' @param newfields Named vector of fields to replace
#' @param Force Overwrite existing file (default FALSE)
#' @param Detached Write a detached header insted of a nrrd (default FALSE)
#' @param action addreplace (Default) addonly or replaceonly
#' @return TRUE or FALSE depending on success
#' @export
#' @seealso \code{\link{ReadNrrdHeader}}
#' @examples
#' \dontrun{
#' AddOrReplaceNrrdHeaderField(lhmaskfile,outfile=file.path(tmpdir,"LHMask.nrrd"),
#'   c("# My interesting comment",`space origin`="(2,2,2)"),Force=TRUE)
#' }
AddOrReplaceNrrdHeaderField<-function(infile,outfile,newfields,Force=FALSE,
	Detached=FALSE,action=c("addreplace","addonly","replaceonly")){
	# see if a given field exists and add or replace its value
	saveontop=ifelse(infile==outfile,TRUE,FALSE)
	if(!Force && file.exists(outfile)) stop("Use Force=TRUE to replace existing files")
  if(saveontop && Detached) stop("Unable to save a detached header on top of its nrrd file")
	action=match.arg(action)

	if(Detached){
		# make a detached header for the original file but don't write it
    tempheader=tempfile(tmpdir=dirname(outfile))
		oht=NrrdMakeDetachedHeaderForNrrd(infile,tempheader)
    inh=ReadNrrdHeader(tempheader)
    unlink(tempheader)
	} else {
		inh=ReadNrrdHeader(infile)
		oht=attr(inh,"headertext")
	}
	
	if(is.null(names(newfields))) names(newfields) <- rep("",length(newfields))
	for(i in seq(newfields)){
		field=names(newfields)[i]
		value=newfields[i]
		if(field==""){
			# this is a comment
			newFieldLine=value
		} else {
			if(!.validateNrrdFieldName(field))
				stop("Invalid nrrd field name: ",field)
			newFieldLine=paste(field,": ",value,sep="")
		}

		if(field%in%names(inh)) {
			# replace existing field
			if(action=="addonly") {
				warning("Unable to replace field in addonly mode")
				return(FALSE)
			}
			oht=sub(paste(field,": .*",sep=""),newFieldLine,oht)
		} else {
			if(action=="replaceonly") {
				warning("Unable to add field in replaceonly mode")
				return(FALSE)
			}
			# just append
			oht=c(oht,newFieldLine)
		}
	}
	
	# add a blank line
  if(Detached){
    writeLines(oht,outfile)
    return(TRUE)
  } 
  oht=c(oht,"")
	headerfile=tempfile()

  writeLines(oht,headerfile)
	if(saveontop){
		outfile=tempfile(pattern=basename(infile),tmpdir=dirname(infile))
	}
	rval=system(paste("unu data",shQuote(infile),"| cat",headerfile,"- >",shQuote(outfile)))
	unlink(headerfile)
	if(rval!=0){
		if(saveontop) unlink(outfile) # cleanup temporary nrrd
		stop("Error ",rval," saving file to: ",outfile)
	}
	# else success
	if(saveontop) file.rename(outfile,infile)
	return(TRUE)
}

.validateNrrdFieldName<-function(fieldname) {
	fieldname=.standardNrrdFieldName(fieldname)
	
	all(fieldname %in% c("space", "space dimension", "space units", "space origin", 
	"space directions", "measurement frame", "dimension", "type", 
	"blocksize", "encoding", "endian", "content", "min", "max", "oldmin", 
	"oldmax", "datafile", "lineskip", "byteskip", "number", "sampleunits", 
	"sizes", "spacings", "thicknesses", "axismins", "axismaxs", "centers", 
	"labels", "units", "kinds"))
}

.standardNrrdFieldName<-function(fieldname,Validate=FALSE)
{
	if(length(fieldname)>1) return(sapply(fieldname,.standardNrrdFieldName,Validate=Validate))
	if(!fieldname%in%c("space dimension","space units","space origin","space directions","measurement frame"))
		fieldname=gsub(" ","",fieldname,fixed=TRUE)
	if(Validate){
		# check that we have been given a valid field
		if(!.validateNrrdFieldName(fieldname))
			stop("Invalid nrrd field: ",fieldname)
	}
	fieldname
}

#' Make a detached header for a specified nrrd file
#'
#' If nhdr is not supplied defaults to <nrrd>.nhdr.
#' If nhdr=NA then new header is returned but not written.
#' @param nrrd Full path to a nrrd file
#' @param nhdr Full path nhdr file to be written
#' @return invisibly returned character vector with new header
#' @export
NrrdMakeDetachedHeaderForNrrd<-function(nrrd,nhdr=paste(nrrd,sep='.','nhdr')){
	h=ReadNrrdHeader(nrrd)
	# drop the directory if the nhdr will be next to the nrrd
	if(is.na(nhdr) || dirname(nrrd)==dirname(nhdr))
		nrrd=basename(nrrd)
	oht=attr(h,'headertext')
	# line skip should be length of old header + 1 for the blank line before data 
	nht=c(oht,paste("line skip:",length(oht)+1))
	nht=c(nht,paste("datafile:",nrrd))
	if(!is.na(nhdr)) writeLines(nht,nhdr)
	invisible(nht)
}

#' Return voxel dimensions (by default absolute voxel dimensions)
#' 
#' NB Can handle off diagonal terms in space directions matrix, 
#' BUT assumes that space direction vectors are orthogonal. 
#' @param f path to nrrd/nhdr file or a list containing a nrrd header
#' @param ReturnAbsoluteDims Defaults to returning absolute value of dims even if
#'        there are any negative space directions
#' @return voxel dimensions as numeric vector
#' @author jefferis
#' @seealso \link{\code{ReadNrrdHeader}}
#' @export
NrrdVoxDims<-function(f,ReturnAbsoluteDims=TRUE,Verbose=FALSE){
  .Deprecated('nat::nrrd.voxdims')
  nrrd.voxdims(f,ReturnAbsoluteDims=ReturnAbsoluteDims)
}
