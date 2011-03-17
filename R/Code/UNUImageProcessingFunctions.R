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

NrrdMinMax<-function(filename,...){
	minmax=.callunu("minmax",shQuote(filename),intern=TRUE,...)
	if(length(minmax)==0) return(c(NA,NA))
	as.numeric(sub("(min|max): ","",minmax))
}

NrrdResample<-function(infile,outfile,size,otherargs=NULL,...){
	if(is.integer(size)) size=paste("--size",paste(size,collapse=" "))
	else size=paste("--size",paste("x",size,sep="",collapse=" "))
	 
	.callunu("resample",paste(size, paste(otherargs,collapse=" "),
	 	"-i",shQuote(infile),"-o",shQuote(outfile)),...)
}

.callunu<-function(cmd,args,unu="unu",DryRun=TRUE,...){
	fullcmd=paste(unu,cmd,paste(args,collapse=" "))
	if(DryRun) print(fullcmd)
	else system(fullcmd,...)
}

NrrdHisto<-function(infile,outfile=sub("\\.([^.]+)$",".histo.\\1",infile),maskfile,bins,min,max,...){
	if (missing(min) || missing(max)) {
		# calculate the minimum and maximum
		r=NrrdMinMax(infile)
		if(any(is.na(r))) stop("Unable to find min and max from: ",infile)
		if(missing(min)) min=r[1]
		if(missing(max)) max=r[2]
	}
	if(missing(bins)){
		# check if this is a float data type
		if(ReadNrrdHeader(infile)$type%in%c("float","double"))
			bins=1000
		else {
			bins=as.integer((max-min)+1)
			if(bins>2^16) bins=1000
		}
	}
	options=paste("-b",bins,"-min",min,"-max",max)
	if(!missing(maskfile)) options=paste(options,"-w",shQuote(maskfile))
	.callunu("histo",paste(options,"-i",shQuote(infile),"-o",shQuote(outfile)),...)
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

NrrdProject<-function(infile,outfile,axis,
	measure=c("max", "min", "mean", "median", "mode", "variance", "skew",
	"intc", "slope", "error", "sd", "product", "sum", "L1", "L2", "Linf"),
	scale="x0.3333 x0.333",
	suffix=NULL,
	CreateDirs=TRUE,Verbose=TRUE,Force=FALSE,UseLock=FALSE){
	measure=match.arg(measure)
	if (missing(outfile)) {
		if(is.null(suffix)) suffix=paste("-",axis,measure,sep="")
		outfile=sub("\\.nrrd$",paste(suffix,".png",sep=""),infile)
	}
	if(!file.exists(infile)) stop("infile: ",infile," does not exist")
	# return TRUE to signal output exists (we just didn't make it now)
	if(!Force && !RunCmdForNewerInput(NULL,infile,outfile)) return (TRUE)
	if(CreateDirs && !file.exists(dirname(outfile))) dir.create(dirname(outfile),recursive = TRUE)
	lockfile=paste(outfile,sep=".","lock")
	# return FALSE to signal output doens't exist
	if(UseLock && !makelock(lockfile)) return (FALSE)
	if(is.numeric(scale)) scale=paste(scale,collapse=" ")
	cmd=paste("unu resample -s",scale," = -k cheap -i",shQuote(infile),
		"| unu project -a",axis,"-m ",measure," | unu quantize -b 8 | unu save -f png",
		"-o",shQuote(outfile))
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
