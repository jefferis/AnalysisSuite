# Support for IO with ImageJ

# LUT functions checked agaist all current Fiji LUTs with:
# read followed by write-binary and write-text and re-read
# NB there seems to be no provision in imagej for rgba luts

WriteImageJLUT<-function(filename,rgb,minmax=c(0,255),ftype=c("text","binary")){
	
	ftype=match.arg(ftype)
	if(is.vector(rgb)) rgb=t(col2rgb(rgb))
	if(ncol(rgb)==4)
		rgb=rgb[,-4]
	if(ncol(rgb)!=3) stop("LUT must have 3 columns")
	if(nrow(rgb)!=256) warning("LUT should have 256 levels")
	
	cmaprange=range(rgb)
	if(cmaprange[1]<minmax[1] || cmaprange[2]>minmax[2]) 
		stop ("input values must be between ",minmax[1]," and ",minmax[2])
	# adjust to 0 minimum and 255 max
	rgb=scale(rgb,center=rep(minmax[1],3),scale=rep(diff(minmax)/255,3))

	if(ftype=='text')
		write.table(rgb,col.names=F,row.names=F,file=filename,sep="\t")
	else
		writeBin(as.integer(rgb),filename,size=1)
}

ReadImageJLUT<-function(filename,ftype=c('guess','binary','nihimage','text'),returnRaw=FALSE){
	
	ftype=match.arg(ftype)
	if(ftype=='guess'){
		fsize=file.info(filename)$size
		if(is.na(fsize)) stop("Unable to read file ",filename)
		magic=readBin(filename,what=integer(),size=4,endian='big')
		if(magic==1229147980) ftype='nihimage'
		else if(fsize==768 || fsize==970) ftype='binary'
		else ftype='text'
	}
	
	if(ftype=='text'){
		linesToSkip=0
		firstLineCheck=try(scan(filename,nmax=3,quiet=TRUE),silent=TRUE)
		if(inherits(firstLineCheck,"try-error")) linesToSkip=1
		
		lut=read.table(filename,skip=linesToSkip,nrow=256)
		if(nrow(lut)!=256) warning("Text LUT ",filename,' has wrong number of rows')
		if(nrow(lut)>256) lut=lut[1:256,]
		if(!(any(ncol(lut)==3:4))) stop("Text LUT ",filename,'has wrong number of columns')
		# 1st column contains indices, so discard
		if(ncol(lut)==4) lut=lut[,-1]		
	} else {
		con=file(filename,'rb')
		nColors=256
		if(ftype=="nihimage"){
			magic=readBin(con,what=integer(),size=4,endian='big')
			if(magic!=1229147980) stop("This is not an NIH Image LUT")
			
			version = readBin(con,what=integer(),size=2,endian='big')
			nColors = readBin(con,what=integer(),size=2,endian='big')
			start = readBin(con,what=integer(),size=2,endian='big')
			end = readBin(con,what=integer(),size=2,endian='big')
			fill1 = readBin(con,what=integer(),size=8,endian='big')
			fill2 = readBin(con,what=integer(),size=8,endian='big')
			filler = readBin(con,what=integer(),size=4,endian='big')
			
			if(nColors != 256) stop("Don't yet know how to interpolate NIH Image LUTs")
		}
		lut=matrix(nrow=nColors,ncol=3)
		lut[,1]=readBin(con,what=integer(),size=1,n=nColors,signed=FALSE)
		lut[,2]=readBin(con,what=integer(),size=1,n=nColors,signed=FALSE)
		lut[,3]=readBin(con,what=integer(),size=1,n=nColors,signed=FALSE)
		close(con)
	}
	if(returnRaw) lut
	else rgb(lut,maxColorValue=255)
}
