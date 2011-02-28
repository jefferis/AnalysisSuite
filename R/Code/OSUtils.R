# Functions to assist in coordinating parallel processing of data
# eg across many cpus/machines or other OS interactions

# inspired by makelock = munger.pl, which is supposed to be
# nfs safe (although I have my doubts)
makelock<-function(lockfile,lockmsg,CreateDirectories=TRUE){
	lockdir=dirname(lockfile)
	if(!file.exists(lockdir)){
		if(CreateDirectories) dir.create(lockdir,recursive=TRUE)
		else stop("Lock Directory for lockfile ",lockfile," does not exist")
	} 
	if(missing(lockmsg)) lockmsg=paste(system('hostname',intern=TRUE),Sys.getenv("R_SESSION_TMPDIR"))
	if (file.exists(lockfile)) return (FALSE)
	# note the use of paste makes the message writing atomic
	cat(paste(lockmsg,"\n",sep=""),file=lockfile,append=TRUE,sep="")
	firstline=readLines(lockfile,n=1)
	if(firstline!=lockmsg){
		# somebody else got there first
		return(FALSE)
	} else return(TRUE)
}

removelock<-function(lockfile){
	if(unlink(lockfile)!=0) {
		warning("Unable to remove ",lockfile)
		return (FALSE)
	}
	return (TRUE)
}

RunCmdForNewerInput<-function(cmd,infiles,outfile,Verbose=FALSE,...){
	if(!all(file.exists(infiles))){
		if(Verbose) cat("some input files missing: ",infiles[!file.exists(infiles)],"\n")
		return (FALSE)
	} else if(!file.exists(outfile)){
		# do nothing just fall through to end
		if(Verbose) cat("outfile: ",outfile,"missing\n")
	} else if(max(file.info(infiles)$mtime) < file.info(outfile)$mtime){
		# check times
		if(Verbose) cat("Skipping",outfile,"because input files are older; use OverWrite=\"yes\" to force\n")
		return(FALSE)	
	} else {
		if(Verbose){
			cat("Overwriting",outfile,"because 1 or more input files are newer\n")
			cat("Newest input mtime:",max(file.info(infiles)$mtime),
				"Output mtime:",file.info(outfile)$mtime,"\n")
		} 
	}
	if(!is.null(cmd)) system(cmd,...)
	return(TRUE)
}

ncpus<-function(default=1){
	os=R.version$os
	if(!is.na(pmatch("darwin",os)))
		return(as.integer(sub(".*:","",system("sysctl hw.ncpu",intern=T))))
	else if(!is.na(pmatch("linux",os)))
		return(as.integer(system("cat /proc/cpuinfo | grep processor | wc -l",intern=T)))
	else {
		warning("I don't know how to check cpu number on OS",os,"defaulting to",default)
		return(default)
	}
}

makehardlink=function(from,to,DryRun=FALSE,Force=FALSE){
	if(length(to)>1) stop("can only have one target")
	# fix any paths using e.g. ~ for home dir - necessary because we will 
	# quote the paths later on which will prevent the shell from expanding them
	to=path.expand(to);from=path.expand(from)
	# handle multiple froms
	if(length(from)>1){
		if(!file.info(".")$isdir) stop("target (to) must be a directory for multiple sources (from)")
		from=paste(shQuote(from),collapse=" ")
	} else from = shQuote(from)
	if(nchar(from)>20000) stop("Shell command length exceeded!")
	cmd=paste("ln",ifelse(Force,"-f",""),from,shQuote(to))
	if(DryRun) cat("I would run:",cmd,"\n")
	else system(cmd)
}

swapfilenames<-function(f1,f2){
	# quick function to swap filenames 
	if(length(f1)>1 || length(f2)>1) return(mapply(swapfilenames,f1,f2))
	
	if(!all(file.exists(f1),file.exists(f2))) stop("f1 and f2 must exist")
	
	tmpfile=tempfile(basename(f1),dirname(f1))
	rval=file.rename(from=f1,to=tmpfile) && file.rename(from=f2,to=f1) && file.rename(from=tmpfile,to=f2)
	return(rval)
}
