# Functions for reading/writing VTK format data
# currently limited to VTK point formats
# for use with Daniel Rueckert's IRTK toolkit

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

ReadVTKLandmarks<-function(filename){
	if(!file.exists(filename)) stop("Cannot read: ",filename)
	con=file(filename,open='rb',encoding='ASCII')
	magic=readLines(con,n=1)
	if(regexpr("# vtk DataFile Version [23]",ignore=T)<0)
		stop("Bad header line in file: ",filename)
	
	title=readLines(con,1)
	encoding=readLines(con,1)
	if(regexpr("ASCII",ignore.case=TRUE)<0){
		
	}
	
	cat("# vtk DataFile Version 2.0",
		title,
		"ASCII",
		"DATASET POLYDATA",
		paste("POINTS",nummarkers,datatype),sep="\n",file=filename)

	write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
}
