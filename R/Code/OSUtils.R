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

#' Run a command if input files are newer than outputs
#'
#' cmd can be an R expression, which is evaluated if necessary,
#' a string to be passed to \code{\link{system}} or
#' NULL/NA in which cases the files are checked and TRUE or FALSE is returned
#' depending on whether action is required.
#'
#' When UseLock=TRUE, the lock file created is called outfiles[1].lock
#' 
#' @param cmd An \code{\link{expression}}, a string or NA/NULL
#' @param infiles Character vector of path to one or more input files 
#' @param outfiles Character vector of path to one or more output files 
#' @param Verbose Write information to consolse (Default FALSE)
#' @param UseLock Stop other processes working on this task (Default FALSE)
#' @return logical indicating if cmd was run or for an R expression, eval(cmd)
#' @export
#' @seealso \code{\link{makelock}}
#' @examples \dontrun{
#' RunCmdForNewerInput(expression(myfunc("somefile")))
#' }
RunCmdForNewerInput<-function(cmd,infiles,outfiles,Verbose=FALSE,UseLock=FALSE,...){
	# note that cmd can be an R expression as in 
	# RunCmdForNewerInput(expression(myfunc("somefile")))
	if(length(infiles)==0){
		if(Verbose) cat("no input files\n")
		return (FALSE)
	} else if(!all(fei<-file.exists(infiles))){
		if(Verbose) cat("some input files missing: ",infiles[!fei],"\n")
		return (FALSE)
	} else if(!all(feo<-file.exists(outfiles))){
		# do nothing just fall through to end
		if(Verbose) cat("outfiles: ",outfiles[!feo],"missing\n")
	} else if( (mit<-max(file.info(infiles)$mtime)) <=
						 (mot<-min(file.info(outfiles)$mtime)) ){
		# check times
		if(Verbose) cat("Skipping",outfiles,"because input files are older\n")
		return(FALSE)	
	} else {
		if(Verbose){
			cat("Overwriting",outfiles,"because 1 or more input files are newer\n")
			cat("Newest input mtime:",strftime(mit),
				"Oldest output mtime:",strftime(mot),"\n")
		} 
	}
	lockfile=paste(outfiles[1],sep=".","lock")
	# return FALSE to signal output doesn't exist
	if(UseLock){
		if(makelock(lockfile))
			on.exit(unlink(lockfile))
		else {
			if(Verbose) cat("Skipping ",outfiles," because someone else is working on ",
				ifelse(length(outfiles)==1,"it","them"),"\n",sep="")
			return(FALSE)
		}
	}
	if(is.expression(cmd)){
		return(eval(cmd))
	} else if(is.character(cmd)){
		system(cmd,...)
	}
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

abs2rel<-function(path,stempath,StopIfNoCommonPath=FALSE){
	# just return the part of path that is NOT in common with stempath
	# eg if path = "/Volumes/JData/JPeople/Sebastian/images"
	# and stempath = "/Volumes/JData/JPeople" (with or without slash)
	# returns "Sebastian/images"
	
	path=path.expand(path)
	stempath=fix.dir(path.expand(stempath))
	
	relpath=sub(stempath,"",path,fixed=TRUE)

	warnorstopfun=if(StopIfNoCommonPath) stop else warning
	if(relpath==path)
		warnorstopfun("stempath: ",stempath,"is not present in: ",path)

	relpath
}

is.dir<-function(dir){
	grepl("/$",dir)
}

fix.dir<-function(dir){
	ifelse(is.dir(dir),yes=dir,no=paste(dir,"/",sep=""))
}

rsync<-function(sourceDir, destinationDir,rsyncoptions="-va"){
	cmd<-paste("rsync",rsyncoptions,shQuote(sourceDir),shQuote(destinationDir))
	system(cmd)
}

#' Use unix touch utility to change file's timestamp
#' 
#' If neither a time or a reference file is provided then the current time is 
#' used. If the file does not already exist, it is created unless Create=FALSE.  
#' @param file Path to file to modify
#' @param time Absolute time in POSIXct format as 
#' @param reference Path to a reference file
#' @param timestoupdate "access" or "modification" (default both)
#' @param Create Logical indicating whether to create file (default TRUE)
#' @return TRUE or FALSE according to success
#' @author jefferis
#' @export
touch<-function(file,time,reference,timestoupdate=c("access","modification"),
    Create=TRUE){
	if(.Platform$OS.type!="unix") {
		warning("touch relies on the existence of a system touch command")
		return(FALSE)
	}
  if(!Create && !file.exists(file)) stop("Create=F and ",file," does not exist") 
	if(!missing(time) && !missing(reference))
		stop("Please supply either a time or a reference file but not both")
	args=paste("-",substr(timestoupdate,1,1),sep="")
	if(!missing(time)){
		if(!is.character(time)) time=strftime(time,"%Y%m%d%H%M.%S")
		args=c(args,"-t",time)
	} else if(missing(reference)) {
		# use current time
	} else {
		# use reference file to supply time
		if(!file.exists(reference)) stop("reference file: ",reference," missing")
		args=c(args,"-r",shQuote(reference))
	}
	
	cmd=paste("touch",paste(args,collapse=" "),shQuote(file))
	return(system(cmd)==0)
}
