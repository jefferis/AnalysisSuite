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
	cat(lockmsg,"\n",file=lockfile,append=TRUE)
	firstline=readLines(lockfile,n=1)
	if(firstline!=lockmsg){
		# somebody else got there first
		return(FALSE)
	} else return(TRUE)
}

removelock<-function(lockfile){
	if(unlink(lockfile)!=1) {
		warning("Unable to remove ",lockfile)
		return (FALSE)
	}
	return (TRUE)
}
