ReadNrrdHeader<-function(filename,Verbose=TRUE,CloseConnection=TRUE){
	if(!inherits(filename,"connection")) con<-file(filename,open='rt')
	else con=filename
	if(CloseConnection) on.exit(close(con))
	# Look for empty line signifying end of header
	headerLines=readLines(con,1)
	NRRDMAGIC="NRRD000"
	if(substring(headerLines,1,nchar(NRRDMAGIC))!=NRRDMAGIC)
		stop("This does not appear to be a NRRD file: ",summary(con)$description)
	nrrdspec=list()
	nrrdkeyvals=vector('character')
	while( (l<-readLines(con,1))!=""){
		headerLines=c(headerLines,l)
		if(substring(l,1,1)=="#") next
		
		if(length(grep(": ",l))>0){
			# field
			hingepos=regexpr(": ",l,fixed=TRUE)
			fieldname=substring(l,1,hingepos-1)
			# make canonical name by removing spaces if required
			if(!fieldname%in%c("space dimension","space units","space origin","space directions","measurement frame"))
			fieldname=gsub(" ","",fieldname,fixed=TRUE)
			
			fieldval=substring(l,hingepos+2,nchar(l))
			
			if(substring(fieldval,1,1)=="("){
				# this is a vector, first remove all spaces
				fieldval=gsub(" ","",fieldval)
				# then remove first and last brackets
				fieldval=substring(fieldval,2,nchar(fieldval)-1)
				vectorstring=unlist(strsplit(fieldval,")(",fixed=TRUE))
				tc=textConnection(vectorstring)
				fieldval=scan(tc,sep=",",quiet=TRUE)
				if(length(vectorstring)>1)
				fieldval=matrix(fieldval,byrow=TRUE,nrow=length(vectorstring))
				close(tc)
			} else if(!fieldname%in%c("type")){
				if (length(grep("^[\\-+]{0,1}[0-9.]+",fieldval,perl=T))>0) what=0
				else what=""
				tc=textConnection(fieldval)
				fieldval=scan(tc,quiet=TRUE,what=what)
				close(tc)
			}
			nrrdspec[[fieldname]]=fieldval
			
		} else if(length(grep(":=",l))>0){
			# key val
			hingepos=regexpr(":=",l,fixed=TRUE)
			nrrdkeyvals[substring(l,1,hingepos-1)]=substring(l,hingepos+2,nchar(l))
		} else {
			warning("Skipping malformed line #",length(headerLines)," in NRRD header\n")
		}
	}
	attr(nrrdspec,'headertext')=headerLines
	attr(nrrdspec,'keyvals')=nrrdkeyvals
	nrrdspec
}

.standardNrrdType<-function(type){
	if(type%in%c("float","double","block")) return (type)
	if(type%in%c("signed char", "int8", "int8_t")) return("int8")
	if(type%in%c("uchar", "unsigned char", "uint8", "uint8_t")) return("uint8")
	if(type%in%c("short", "short int", "signed short", "signed short int", "int16", "int16_t")) return("int16")
	if(type%in%c("ushort", "unsigned short", "unsigned short int", "uint16", "uint16_t")) return("uint16")
	if(type%in%c("int", "signed int", "int32", "int32_t")) return("int32")
	if(type%in%c("uint", "unsigned int", "uint32", "uint32_t")) return("uint32")
	if(type%in%c("longlong", "long long", "long long int", "signed long long", "signed long long int", "int64", "int64_t"))
		return("int64")
	if(type%in%c("ulonglong", "unsigned long long", "unsigned long long int", "uint64", "uint64_t")) return("uint64")
	return(NULL)
}

Read3DDensityFromNrrd<-function(filename,Verbose=FALSE,AttachFullHeader=FALSE,
	ReadByteAsRaw=c("unsigned","all","none"),origin){
	ReadByteAsRaw=match.arg(ReadByteAsRaw)
	fc=file(filename,'rb')
	h=ReadNrrdHeader(fc,CloseConnection=FALSE)
		
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
	
	dataLength=prod(h$sizes)
	endian=ifelse(is.null(h$endian),.Platform$endian,h$endian)
	if(Verbose) cat("dataLength =",dataLength,"dataType =",dataTypes$what[i],"size=",dataTypes$size[i],"\n")
	enc=tolower(h$encoding)
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
	if(AttachFullHeader) attr(d,"header")=h
	
	if('space directions'%in%names(h)){
		voxdims=rowSums(sqrt(h[['space directions']]^2))
	} else if ('spacings'%in%names(h)){
		voxdims=h[["spacings"]]
	} else {
		# no pixel size info, so just return
		return(d)
	}
	latticeBoundingBox=rbind(c(0,0,0),(h$sizes-1)*voxdims)
	if(!missing(origin)){
		latticeBoundingBox=t(origin+t(latticeBoundingBox))
	} else if('space origin'%in%names(h)){
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

ReadHistogramFromNrrd<-function(filename,...){
	d=Read3DDensityFromNrrd(filename,AttachFullHeader=TRUE,...)
	h=attr(d, "header")
	# clear the header attributes
	attributes(d)<-NULL
	if(is.na(pmatch("histo",h$content)) || h$dimension!=1) {
		warning ("This does not appear to be a 1d nrrd histogram")
		return(d)
	}
	breaks=seq(from=h$axismins,to=h$axismax,len=h$sizes+1)
	density=d/sum(d)
	
	halfwidth=(h$axismaxs-h$axismins)/h$sizes/2

	# return it as an R histogram	
	structure(list(
		breaks = breaks, 
		counts = d, 
	intensities = density, 
	density = density, 
	mids = seq(h$axismins+halfwidth,h$axismaxs-halfwidth,len=h$sizes), 
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

AddOrReplaceNrrdHeaderField<-function(infile,outfile,fields,values,Force=FALSE,action=c("addreplace","addonly","replaceonly")){
	# see if a given field exists and add or replace its value
	if (infile==outfile) stop("AddOrReplaceNrrdHeaderField: Cannot currently save on top of existing file")
	if(!Force && file.exists(outfile)) stop("Use Force=TRUE to replace existing files")
	inh=ReadNrrdHeader(infile)
	action=match.arg(action)
	
	if(length(fields)!=length(values)) stop("Must supply same number of fields and values")
	
	for(i in seq(fields)){
		field=fields[i]
		value=values[i]
		newFieldLine=paste(field,": ",value,sep="")
		oht=attr(inh,"headertext")

		if(field%in%names(inh)) {
			# replace existing field
			if(action=="addonly") {
				warning("Unable to replace field in addonly mode")
				return(FALSE)
			}
			oht=sub(paste(field,": .*",sep=""),newFieldLine,oht)
		} else {
			if(action=="replaceonly") {
				warning("Unable to replace field in replaceonly mode")
				return(FALSE)
			}
			# just append
			oht=c(oht,newFieldLine)
		}		
	}
	
	# add a blank line
	oht=c(oht,"")
	tmpheader=tempfile()
	writeLines(oht,tmpheader)
	system(paste("unu data",shQuote(infile),"| cat",tmpheader,"- >",shQuote(outfile)))
	unlink(tmpheader)
}

.standardNrrdFieldName<-function(fieldname)
{
	if(length(fieldname)>1) return(sapply(fieldname,.standardNrrdFieldName))
	if(!fieldname%in%c("space dimension","space units","space origin","space directions","measurement frame"))
	fieldname=gsub(" ","",fieldname,fixed=TRUE)
	fieldname
}

is.nrrd<-function(f,ReturnVersion=FALSE,TrustSuffix=FALSE){
	# TrustSuffix => expect files to end in nrrd or nhdr
	if(TrustSuffix && ReturnVersion) 
		stop("Cannot use return nrrd version without reading file to check nrrd magic")

	if(TrustSuffix)
		return(grepl("\\.n(hdr|rrd)$",ff,ignore.case=TRUE))
	
	if(length(f)>1)
		return(sapply(f,is.nrrd,ReturnVersion=ReturnVersion))
	
	nrrd=as.raw(c(0x4e,0x52,0x52,0x44))
	magic=readBin(f,what=nrrd,n=8)
	if(any(magic[1:4]!=nrrd))
		return (FALSE)

	if(ReturnVersion)
		return(as.integer(magic[8]))

	TRUE
}