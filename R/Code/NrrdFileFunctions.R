ReadNrrdHeader<-function(filename,Verbose=TRUE,CloseConnection=TRUE){
	if(!inherits(filename,"connection")) con<-file(filename,open='rt')
	else con=filename
	if(CloseConnection) on.exit(close(con))
	# Look for empty line signifying end of header
	headerLines=readLines(con,1)
	NRRDMAGIC="NRRD000"
	if(substring(headerLines,1,nchar(NRRDMAGIC))!=NRRDMAGIC)
		stop("This does not appear to be a NRRD file: bad magic")
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

Read3DDensityFromNrrd<-function(filename,Verbose=FALSE,AttachFullHeader=FALSE,origin){

	fc=file(filename,'rb')
	h=ReadNrrdHeader(fc,CloseConnection=FALSE)
		
	# now read the data
	dataTypes=data.frame(name=I(c("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64",
	"float", "double", "block")),
			size=c(1,1,2,2,4,4,8,8,4,8,NA),what=I(c(rep("integer",8),rep("numeric",2),"raw")),
			signed=c(rep(c(T,F),4),rep(T,3)))
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
		voxdims=spacings
	} else {
		# no pixel size info, so just return
		return(d)
	}
	latticeBounds=rbind(c(0,0,0),h$sizes*voxdims)
	if(!missing(origin)){
		latticeBounds=t(origin+t(latticeBounds))
	} else if('space origin'%in%names(h)){
		latticeBounds=t(h[['space origin']]+t(latticeBounds))
	}
	attr(d,"BoundingBox")<-latticeBounds
	attr(d,"x")<-seq(latticeBounds[1],latticeBounds[2],len=h$sizes[1])
	attr(d,"y")<-seq(latticeBounds[3],latticeBounds[4],len=h$sizes[2])
	attr(d,"z")<-seq(latticeBounds[5],latticeBounds[6],len=h$sizes[3])

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

Nrrd2op<-function(infiles,outfile,unu2opfun=c("max","min","+", "-", "x", "/"),CreateDirs=TRUE,Verbose=TRUE,Force=FALSE){
	if(length(infiles)<1) return(NULL)
	# Do nothing if inputs are older than output unless Force=T
	if(!Force && !RunCmdForNewerInput(NULL,infiles,outfile)) return (NULL)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile))
	file.copy(infiles[1],outfile)
	if(length(infiles)==1) return(1)
	unu2opfun=match.arg(unu2opfun)
	for (f1 in infiles[-1]){
		cmd=paste("unu 2op",unu2opfun,
			shQuote(f1),shQuote(outfile),"-o",shQuote(outfile))
		system(cmd)
		if(Verbose) cat(".")
	}
}

NrrdMinMax<-function(filename,...){
	minmax=.callunu("minmax",shQuote(filename),...)
	minmax=sub("^.*")
}

NrrdResample<-function(infile,outfile,size,otherargs=NULL,...){
	if(is.integer(size)) size=paste("--size",paste(size,collapse=" "))
	else size=paste("--size",paste("x",size,sep="",collapse=" "))
	 
	.callunu("resample",paste(size, paste(otherargs,collapse=" "),
	 	"-i",shQuote(infile),"-o",shQuote(outfile)),...)
}

.callunu<-function(cmd,args,unu="unu",...){
	system(paste(unu,cmd,paste(args,collapse=" ")),intern=TRUE,...)
}
