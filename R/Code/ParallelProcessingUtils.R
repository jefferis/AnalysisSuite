# Functions to assist in coordinating parallel processing of data
# eg across many cpus/machines

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
	if(!file.exists(outfile)){
		# do nothing just fall through to end
		if(Verbose) cat("outfile: ",outfile,"missing\n")
	} else if(!all(file.exists(infiles))){
		if(Verbose) cat("some input files missing: ",infiles[!file.exists(infiles)],"\n")
		return (FALSE)		
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
