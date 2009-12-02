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
	minmax=.callunu("minmax",shQuote(filename),...)
	as.integer(sub("(min|max): ","",minmax))
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
