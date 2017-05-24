# Some Image Processing Functions making use of the unu command line tool
# unu = Utah Nrrd Utilities - ie they operate on nrrd files
# See http://teem.sourceforge.net/unrrdu/index.html
 
Nrrd2op<-function(infiles,outfile,fun=c("max","min","+", "-", "x", "/"),
	gzip=FALSE,CreateDirs=TRUE,Verbose=TRUE,Force=FALSE){
	if(length(infiles)<1) return(NULL)
	# Do nothing if inputs are older than output unless Force=T
	if(!Force && !RunCmdForNewerInput(NULL,infiles,outfile)) return (NULL)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile))
	fun=match.arg(fun)
	
	# 1 input: f1 -> output
	if(length(infiles)==1) file.copy(infiles[1],outfile)
	else if(length(infiles)>1) {
		# 2 or more inputs: f1 x f2 -> output
		cmd=paste("unu 2op",fun,
			shQuote(infiles[1]),shQuote(infiles[2]),"-o",shQuote(outfile))
		system(cmd)
		if(Verbose) cat(".")
	}
	# 3 or more inputs: f3 ... x output -> output
	if(length(infiles)>2) {
		for (f in infiles[-(1:2)]){
			cmd=paste("unu 2op",fun,
				shQuote(f),shQuote(outfile),"-o",shQuote(outfile))
			system(cmd)
			if(Verbose) cat(".")
		}		
	}
	if(length(infiles)>1 && gzip)
		system(paste("unu save --format nrrd --encoding gzip","-i",shQuote(outfile),"-o",shQuote(outfile)))
	return(length(infiles))
}

#' Return range of values in nrrd file (interface to unu minmax command)
#' @param filename nrrd file containing data
#' @param Blind8 Assume 8 bit data has range 0-255 (or -127-127), default F.
#' @param ... passed to \link{\code{system}} function
#' @return c(min,max) or c(NA,NA) like R's \link{\code{range}} function
#' @author jefferis
#' @export
NrrdMinMax<-function(filename, Blind8=FALSE, ...){
	args=shQuote(filename)
	if(Blind8) args=c(args,"-blind8 true")
	minmax=.callunu("minmax",args,intern=TRUE,...)
	if(length(minmax)==0) return(c(NA,NA))
	as.numeric(sub("(min|max): ","",minmax[1:2]))
}

#' Resample a nrrd with simple or complex filters
#'
#' wraps unu resample
#' 
#' When size is explicitly an integer (e.g. c(512L,512L,100L) then it is 
#' assumed to represent the size of the desired output in pixels.
#' 
#' unu resample defaults to cell centering in the asbence of information, but I 
#' prefer node. 
#' 
#' * downsampling a **cell** image (e.g. x2) will result in 
#'   * an origin shift **iff** an origin was specified in the input file
#'   * a change in the space directions that is not an exact multiple of the downsampling factor 
#' * downsampling a node image (e.g. x2) will result in:
#'   * no change in origin
#'   * a doubling in the space directions field
#' @param infile 
#' @param outfile 
#' @param size Numeric vector of scale factors (NA=> don't touch) or 
#'   integer vector of pixel dimensions. 
#' @param voxdims Target voxel dimensions
#' @param centering Whether to assume node or cell centering if not specified
#'   in input file.
#' @param otherargs Passed to unu resample
#' @param gzip 
#' @param suffix 
#' @param CreateDirs 
#' @param Verbose 
#' @param Force 
#' @param UseLock 
#' @param ... Additional params passed to .callunu
#' @return 
#' @author jefferis
#' @export
NrrdResample<-function(infile,outfile,size,voxdims=NULL,
		centering=c("cell","node"),otherargs=NULL,gzip=TRUE, suffix=NULL,
		CreateDirs=TRUE,Verbose=TRUE,Force=FALSE,UseLock=FALSE,...){

	centering=match.arg(centering)
	if(!file.exists(infile)) stop("infile: ",infile," does not exist")
	if(!is.null(voxdims)){
		# we have a target voxel size instead of a standard size specification
		# first fetch existing voxel size
		size=NrrdVoxDims(infile)/voxdims
	}
	if(is.integer(size)) size=paste("--size",paste(size,collapse=" "))
	else {
		size=paste("--size",paste("x",size,sep="",collapse=" "))
		# eg size=c(0.5,0.5,NA) => "--size x0.5 x0.5 ="
		size=sub("xNA","=",size)
	}
	 
	if (missing(outfile)) {
		# we only want to add a default suffix if we are putting the output
		# file into the same directory
		if(is.null(suffix))
			suffix=gsub(" ","",sub("--size","",size))
		outfile=sub("\\.nrrd$",paste(suffix,"\\.nrrd",sep=""),infile)
	} else {
		# here we have been passed an output directory rather than a file
		if(file.exists(outfile) && file.info(outfile)$isdir)
			outfile=file.path(outfile,
				sub("\\.nrrd$",paste(suffix,"\\.nrrd",sep=""),basename(infile)))
	}

	# return outfile to signal output exists (we just didn't make it now)
	if(!Force && !RunCmdForNewerInput(NULL,infile,outfile)) return (outfile)
	# NB createdirs only makes sense if we have directly specified outfile
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	lockfile=paste(outfile,sep=".","lock")
	# return FALSE to signal output doens't exist
	if(UseLock){
		if(makelock(lockfile))
			on.exit(unlink(lockfile))
		else
			return(FALSE)
	}
	otherargs=c(otherargs,"--center",centering)
	cmdargs=paste(size, paste(otherargs,collapse=" "),"-i",shQuote(infile))
	if(gzip)
		cmdargs=paste(cmdargs,"| unu save -f nrrd -e gzip")
	rval=.callunu("resample",paste(cmdargs,"-o",shQuote(outfile)),...)
	if(rval!=0) stop("error ",rval," in unu resample")
	return(outfile)
}

.callunu<-function(cmd,args,unu="unu",DryRun=FALSE,...){
	fullcmd=paste(unu,cmd,paste(args,collapse=" "))
	if(DryRun) {
		cat(fullcmd,"\n",sep="")
		return(0)
	}
	else system(fullcmd,...)
}

NrrdHisto<-function(infile,outfile=sub("\\.([^.]+)$",".histo.\\1",infile),
	maskfile,bins,min,max,blind8=TRUE,...){
	h=ReadNrrdHeader(infile)
	nrrdType=.standardNrrdType(h$type)
	unuhistooptions=''
	if(nrrdType%in%c("uint8","int8") && blind8 && missing(min) && missing(max)){
		# this is an 8 bit image
		# and we're going to use unu histo's blind8 option (default) to
		# assume uchar [0,255], signed char is [-128,127]
		bins=256
	} else {
	    if (missing(min) || missing(max)) {
			# calculate the minimum and maximum
			r=NrrdMinMax(infile)
			if(any(is.na(r))) stop("Unable to find min and max from: ",infile)
			if(missing(min)) min=r[1]
			if(missing(max)) max=r[2]
    	}  
		unuhistooptions=paste("-min",min,"-max",max)
	} 
	if(missing(bins)){
		# check if this is a float data type
		if(nrrdType%in%c("float","double"))
			bins=1000
		else {
			bins=as.integer((max-min)+1)
			if(bins>2^16) bins=1000
		}
	}
	unuhistooptions=paste(unuhistooptions,"-b",bins)
	if(!missing(maskfile)) unuhistooptions=paste(unuhistooptions,"-w",shQuote(maskfile))
	.callunu("histo",paste(unuhistooptions,"-i",shQuote(infile),"-o",shQuote(outfile)),...)
	return(outfile)
}

NrrdQuantize<-function(infile,outfile,min,max,bits=c("8","16","32"),
	gzip=FALSE,CreateDirs=TRUE,Verbose=TRUE,Force=FALSE,UseLock=FALSE){

	# Do nothing if inputs are older than output unless Force=T
	if(!Force && !RunCmdForNewerInput(NULL,infile,outfile)) return (FALSE)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	
	bits=match.arg(bits)
	
	lockfile=paste(outfile,sep=".","lock")
	if(UseLock && !makelock(lockfile)) return (FALSE)
	# unu quantize -b 8 -min 7397.386 -max 65535 -i AL-PNs_mask_MF_final_IS2_SAIA24-1_02.nrrd \
	# | unu save -f nrrd -e gz -o AL-PNs_mask_MF_final_IS2_SAIA24-1_02-quant.nrrd
	
	minmaxopts=paste(ifelse(missing(min),"",paste("-min",min)), 
		ifelse(missing(max),"",paste("-max",max)))
	cmd=paste("unu quantize -b",bits,minmaxopts,"-i",shQuote(infile))
	if(gzip) cmd=paste(cmd,"| unu save -f nrrd -e gzip -o",shQuote(outfile))
	else cmd=paste(cmd,"-o",shQuote(outfile))
	system(cmd)
	if(Verbose) cat(".")
	if(UseLock) unlink(lockfile)
	return(TRUE)
}

NrrdTestIntegrity<-function(infile,defaultReturnVal=TRUE){
	# Tests integrity of a compressed nrrd file using crc check
	if(!file.exists(infile)) return(NA)
	h=ReadNrrdHeader(infile)
	if(tolower(h$encoding)%in%c("gz","gzip")) testprog='gzip'
	else if(tolower(h$encoding)%in%c("bz2","bzip2")) testprog='bzip2'
	else {
		warning("Unable to test integrity of nrrd file with encoding",h$encoding)
		return(defaultReturnVal)
	}
	testval=system(paste("unu data",infile," | ",testprog,"-t"),ignore.stderr=TRUE)
	return(testval==0)
}

NrrdTestDataLength<-function(infile,defaultReturnVal=TRUE){
	# Tests integrity of a nrrd file by checking that the data block is as long
	# as it should be. For a gzip file, it checks that the last 4 bytes
	# encode the length of the uncompressed data. This says nothing about the contents
	# but does ensure that the file has not been truncated which is much the most 
	# common problem. 
	if(!file.exists(infile)) return(NA)
	# need to read this full header in order to get the size in bytes of
	# the uncompressed data
	fullh=Read3DDensityFromNrrd(infile,ReadData=F,AttachFullHeader=T)
	h=attr(fullh,'header')
	if(tolower(h$encoding)%in%c("gz","gzip")) enc='gzip'
	else if(tolower(h$encoding)=='raw') enc='raw'
	else {
		warning("Unable to test data length of nrrd file with encoding",h$encoding)
		return(defaultReturnVal)
	}
	# figure out how many bytes we are expecting
	dataLength=attr(fullh,"datablock")$n*attr(fullh,"datablock")$size
	datafile=nat:::nrrd.datafiles(h)
	if(length(datafile)>1){
		warning("Don't yet know how to handle detached headers with more than one data file")
		return(defaultReturnVal)
	}
	if(enc=='gzip'){
		# Check gzip file
		# last 4 bytes of gzip contain length of compressed stream
		con=file(datafile,open='rb')
		on.exit(close(con))
		seek(con,-4,origin='end')
		# how do we handle endianess?
		uncompressedBytes=readBin(con,what=integer(),size=4,n=1)
		if(uncompressedBytes>(2^32-1)) {
			warning("Don't know how to check integrity of files >= 2^32 bytes")
			return(defaultReturnVal)
		}
		return(dataLength==uncompressedBytes)
	} else {
		# TODO fix handling of detached nrrd files with raw encoding
		start=attr(fullh,'datablock')$datastartpos
		end=file.info(datafile)$size
		# TODO decide what we should say if we have MORE than enough data
		return(dataLength>=(end-start))
	}
}

NrrdCrc<-function(infile,UseGzip=FALSE,FastHeader=TRUE){
	# gets the CRC (hash) of a gzip encoded nrrd
	# Defaults to a quick method based on
	# knowledge of gzip file format from:
	# http://www.gzip.org/zlib/rfc-gzip.html
	# and assumption that there is only one member in gzip data
	# can also use gzip but this is much slower since have to copy
	# unu data to temporary file
	if(!file.exists(infile)) return(NA)

	ext=tolower(sub(".*\\.([^.]+)$","\\1",basename(infile)))
	if(ext=='nhdr') FastHeader=FALSE

	h=NULL
	if(FastHeader){
		# quick and dirty reading of header
		con<-file(infile,'rb')
		on.exit(close(con))
		magic=readBin(con,what=raw(),5)
		if(!isTRUE(all.equal(magic,as.raw(c(0x4e, 0x52, 0x52, 0x44, 0x30))))){
			warning("This is not a nrrd")
			return(NA)
		}
		while( length(l<-readLines(con,1))>0 && l!="" ){
			if(nchar(l)>12 && substr(l,1,10)=="encoding"){
				if(substr(encoding,11,12)!="gz"){
					warning("This is not a gzip encoded nrrd")
					return(NA)
				}
			}
		}
	} else {
		h=ReadNrrdHeader(infile)
		if(tolower(h$encoding)%in%c("gz","gzip")) {
			# testprog='gzip'
		} else {
			warning("This is not a gzip encoded nrrd")
			return(NA)
		}
	}

	if(UseGzip){
		tmp=tempfile()
		on.exit(unlink(tmp),add=TRUE)
		system(paste("unu data ",shQuote(infile)," > ",shQuote(tmp)))
		x=system(paste("gzip -lv",shQuote(tmp)),intern=TRUE)
		crc=try(strsplit(x[2],"[ ]+")[[1]][[2]])
		if(inherits(crc,'try-error')) crc=NA		
	} else {
		# TODO Fix handling of nhdr files
		if(!is.null(h) && !is.null(h$datafile)){
			if(length(h$datafile)>1)
				stop("Don't know how to handle more than 1 datafile")
			if(dirname(h$datafile)==".")
				h$datafile=file.path(dirname(infile),h$datafile)
			con=file(h$datafile,open='rb')
			on.exit(close(con),add=TRUE)
		} else if(!FastHeader){
			con=file(infile,open='rb')
			on.exit(close(con),add=TRUE)
		}
		seek(con,-8,origin='end')
		# TODO check endian issues (what happens if CRC was from opposite endian platform?)
		crc=readBin(con,integer(),size=4)
		crc=format(as.hexmode(crc),width=8)
	}	
	crc
}


#' Make projections along a nrrd axis
#'
#' @details gamma: Just as in xv, the gamma value here is actually the reciprocal
#'   of the exponent actually used to transform the values.
#' Note also that for \code{cropmin,cropmax} the special value M can be used
#'   to indicate the maximum index for that axis (i.e. n-1 when there are n 
#'   samples).
#' @param infile Path to input file
#' @param outfile Optional path to output file (constructed automatically when missing)
#' @param axis Number indicating 0-indexed axis or character "x", "y" or "z"
#' @param measure Character vector indicating summary function to apply to
#'   values in each column. Choose from \code{c("max", "min", "mean", "median",
#'   "mode", "variance", "skew", "intc", "slope", "error", "sd", "product", 
#'   "sum", "L1", "L2", "Linf")}
#' @return Logical indicating success
#' @export
#' @seealso \code{\link{NrrdSave},\link{NrrdMerge}}
#' @examples
#' \dontrun{
#' NrrdProject(infile,axis='z',cropmin='0 0 50',cropmax='0 0 M-50')
#' }
NrrdProject<-function(infile,outfile,axis,
	measure=c("max", "min", "mean", "median", "mode", "variance", "skew",
	"intc", "slope", "error", "sd", "product", "sum", "L1", "L2", "Linf"),
	scale="x0.3333 x0.333",kernel='cheap',gamma=NA,cropmin=NULL,cropmax=NULL,
	suffix=NULL,
	CreateDirs=TRUE,Verbose=TRUE,Force=FALSE,UseLock=FALSE){
	measure=match.arg(measure)
	if (missing(outfile)) {
		if(is.null(suffix)) suffix=paste("-",axis,measure,sep="")
		outfile=sub("\\.nrrd$",paste(suffix,".png",sep=""),infile)
	}
	if(!file.exists(infile)) stop("infile: ",infile," does not exist")
	if(is.character(axis)) {
		axis=match(tolower(axis),c("x",'y','z'))-1
		if(is.na(axis)) stop("Unable to recognise nrrd axis")
	}
	# return TRUE to signal output exists (we just didn't make it now)
	if(!Force && !RunCmdForNewerInput(NULL,infile,outfile)) return (TRUE)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	lockfile=paste(outfile,sep=".","lock")
	# return FALSE to signal output doens't exist
	if(UseLock && !makelock(lockfile)) return (FALSE)
	if(is.numeric(scale)) scale=paste(scale,collapse=" ")
	if(!is.null(cropmax)||!is.null(cropmin)){
		cmd=paste("unu crop",
			ifelse(is.null(cropmax),'',paste("-max",cropmax)),
			ifelse(is.null(cropmin),'',paste("-min",cropmin)),
			"-i",shQuote(infile),
			"| unu resample -s",scale," = -k ",kernel,
			"| unu project -a",axis,"-m ",measure," | unu quantize -b 8 ")
	} else cmd=paste("unu resample -s",scale," = -k ",kernel," -i",shQuote(infile),
		"| unu project -a",axis,"-m ",measure," | unu quantize -b 8 ")

	if(!is.na(gamma)){
		cmd=paste(cmd,"| unu gamma --gamma",gamma)
	}
	cmd=paste(cmd,"| unu save -f png -o",shQuote(outfile))
	rval = system(cmd)
	if(Verbose) cat(".")
	if(UseLock) unlink(lockfile)
	if(rval!=0) stop("unu error ",rval," in NrrdProject")
	return(TRUE)
}

NrrdMerge<-function(infiles,outdir=NULL,outfile=NULL,axis=0,
	CreateDirs=TRUE,Verbose=TRUE,Force=FALSE,UseLock=FALSE,...){
	if(is.null(outfile) && is.null(outdir)) {
		outfile=sub("\\.[^.]+$","-merge.png",infiles[1])
	} else if (is.null(outfile)){
		outfile=file.path(outdir,basename(infiles[1]))
		outfile=sub("\\.[^.]+$",".png",outfile)
	} else {
		# make sure this is a png
		outfile=sub("\\.[^.]+$",".png",outfile)
	}

	if(!all(file.exists(infiles))) stop("one or more infiles: ",infiles," do not exist")
	# return TRUE to signal output exists (we just didn't make it now)
	if(!Force && !RunCmdForNewerInput(NULL,infiles,outfile)) return (TRUE)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	lockfile=paste(outfile,sep=".","lock")
	# return FALSE to signal output doens't exist
	if(UseLock && !makelock(lockfile)) return (FALSE)

	cmd=paste("unu join -a",axis,"-incr -i",paste(shQuote(infiles),collapse=" "),"-o",shQuote(outfile))
	cat(cmd,"\n")
	rval = system(cmd)
	if(Verbose) cat(".")
	if(UseLock) unlink(lockfile)
	if(rval!=0) stop("unu error ",rval," in NrrdProject")
	return(TRUE)
}

NrrdFlip<-function(infile,outfile,axes,suffix=NULL,endian=c("big","little"),
	CreateDirs=TRUE,Verbose=TRUE,UseLock=FALSE, OverWrite=c("no","update","yes")){
	# TODO would be nice if we could 
	# a) have an absolute flip mode that checks the nrrd content field
	# b) similarly checks whether output image has been flipped accordingly
	
	if(is.logical(OverWrite)) OverWrite=ifelse(OverWrite,"yes","no")
	else OverWrite=match.arg(OverWrite)
	
	endian=match.arg(endian)
	if (missing(outfile)) {
		if(is.null(suffix)) suffix=paste("-flip",paste(axes,collapse=""),sep="")
		outfile=sub("\\.nrrd$",paste(suffix,".nrrd",sep=""),infile)
	}
	
	if(!file.exists(infile)) stop("infile: ",infile," does not exist")
	
	# return TRUE to signal output exists (whether or not we made it)
	if(file.exists(outfile)){
		if(OverWrite=="no"){
			if(Verbose) cat("Output",outfile,"already exists; use OverWrite=\"yes\" or \"update\" to overwrite or update\n")
			return(TRUE)
		} else if(OverWrite=="update"){
			# check modification times
			if(!RunCmdForNewerInput(NULL,infile,outfile)) return (TRUE)
		} else if(Verbose) cat("Overwriting",outfile,"because OverWrite=\"yes\"\n")
	}
	
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	lockfile=paste(outfile,sep=".","lock")
	# return FALSE to signal output doesn't (yet) exist
	if(UseLock && !makelock(lockfile)) return (FALSE)
	on.exit(unlink(lockfile))

	# First axis
	cmd=paste("unu flip -a",axes[1],"-i",shQuote(infile))
	# any additional axes
	for(axis in axes[-1]){
		cmd=paste(cmd,"| unu flip -a",axis)
	}
	# save
	cmd=paste(cmd," | unu save -f nrrd -e gz -en",endian,"-o",shQuote(outfile))
	
	rval = system(cmd)
	if(Verbose) cat(".")
	if(rval!=0) stop("unu error ",rval," in NrrdProject")
	return(TRUE)
}

#' Save a nrrd file in a different encoding / format
#'
#' In outfile is missing, will overwrite infile
#' @param infile, outfile Paths to input and output nrrds
#' @param format (one of 6 supported by unu save, default nrrd)
#' @param encoding (one of 5 supported by unu save, default gz when nrrd)
#' @param CreateDirs Whether to make missing direcories implied by output path
#' @param UseLock Whether to make a lockfile (useful for parallel processing)
#' @param DryRun Return command instead of running it
#' @return path to ouput file
#' @export
#' @seealso \code{\link{NrrdCrc}}
NrrdSave<-function(infile,outfile,format=c("nrrd","pnm","text","vtk","png","eps"),
  encoding=ifelse(format=='nrrd','gzip','raw'),
  CreateDirs=TRUE,UseLock=FALSE,DryRun=FALSE){
  format=match.arg(format)
  encodings=c("raw","ascii","hex","gzip","bzip2")
  encoding=encodings[pmatch(encoding,encodings)]
  if(is.na(encoding)) stop("Invalid encoding")
  
  if(missing(outfile)) {
    outfile=infile
    samefile=TRUE
  } else samefile=FALSE
  
  cmd=paste('unu save','--format',format,'-e',encoding,'-i',shQuote(infile),'-o',shQuote(outfile))
  if(DryRun) return(cmd)

  if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
  lockfile=paste(outfile,sep=".","lock")
  # return FALSE to signal output doesn't (yet) exist
  if(UseLock && !makelock(lockfile)) return (FALSE)
  on.exit(unlink(lockfile))

  rval=system(cmd)==0
  attr(rval,'outfile')=outfile
  rval
}
