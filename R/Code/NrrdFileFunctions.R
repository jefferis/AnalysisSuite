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

NrrdMinMax<-function(filename,...){
	minmax=.callunu("minmax",shQuote(filename),...)
	minmax=sub("^.*")
}

.callunu<-function(cmd,args,unu="unu",...){
	system(paste(unu,cmd,paste(args,collapse=" ")),intern=TRUE,...)
}
