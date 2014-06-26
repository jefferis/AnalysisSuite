# Functions to assist in coordinating parallel processing of data
# eg across many cpus/machines or other OS interactions

# most of the code originally in this file has been moved to the new
# nat.utils package, on github/our local R repo
makehardlink<-function(...){
  .Deprecated('file.hardlink')
  file.hardlink(...)
}

swapfilenames<-function(...){
  .Deprecated('nat.utils::file.swap')
  nat.utils::file.swap(...)
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
