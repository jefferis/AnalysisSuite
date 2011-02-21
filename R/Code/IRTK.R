# Some functions to help use Daniel Rueckert's 
# Image Registration Toolkit (IRTK)
# See: http://www.doc.ic.ac.uk/~dr/software

irtk.dof2mat<-function(doffile,matfile,Invert=FALSE,...){
	# see http://www.doc.ic.ac.uk/~dr/software/usage.html#dof2mat
	# reads in a transformation matrix from doffile and converts
	# to 4 x4 homogeneous transformation matrix
	deleteoutput=FALSE
	if(missing(matfile)) {
		matfile=tempfile()
		deleteoutput=TRUE
	}
	
	rval=.callirtk('dof2mat',
		c(shQuote(doffile),"-matout",shQuote(matfile),ifelse(Invert,character(0),"-invert")),...)
	if(rval>0) return(NULL)
	
	mat=scan(matfile,skip=1)
	if(length(mat)!=16) stop("Unable to read a 4x4 matrix. Output was:",readLines(matfile))
	if(deleteoutput) unlink(matfile)
	matrix(mat,nrow=4,byrow=T)
}

.callirtk<-function(cmd,args,irtkdir=path.expand("~/dev/registration/irtk"),Verbose=FALSE,DryRun=FALSE,...){
	cmdline=paste(file.path(irtkdir,cmd),paste(args,collapse=" "))
	if(Verbose || DryRun) print(cmdline)
	if(DryRun) return(FALSE)
	system(cmdline,...)
}

WriteVTKLandmarks<-function(filename,d,title,datatype=c("float","double")){
	if(ncol(d)!=3) stop("Expect N rows x 3 cold of 3d points")
	nummarkers=nrow(d)
	datatype=match.arg(datatype)
	if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
	
	cat("# vtk DataFile Version 2.0",
		title,
		"ASCII",
		"DATASET POLYDATA",
		paste("POINTS",nummarkers,datatype),sep="\n",file=filename)

	write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
}
