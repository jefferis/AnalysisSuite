# Some Image Processing Functions making use of the unu command line tool
# unu = Utah Nrrd Utilities - ie they operate on nrrd files
# See http://teem.sourceforge.net/unrrdu/index.html
 
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
