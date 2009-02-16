# Simple support functions for Mark Longair's landmarks format from
# the VIB/landmarks package 

ReadLongairPoints<-function(f){
	x=readLines(f)
	numbers=sub(".*\\[(.*)\\].*","\\1",x)
	rnames=sub("(\\\"){0,1}(.?*):\\s+\\[.*","\\2",x)
	rnames=sub("\\\"$","",rnames)
	tc=textConnection(numbers)
	points=read.csv(tc,col.names=c("X","Y","Z"),header=FALSE,row.names=rnames,blank.lines.skip = TRUE)
	close(tc)
	points
}

require(XML)
ReadLongairTraceData<-function(f,Verbose=TRUE){
	doc<-xmlTreeParse(f)
	r<-xmlRoot(doc)
	if(xmlName(r)!="tracings") stop("This is not a Longair format tracing")
	
	l<-list()
	stopifnot(names(r)[1:2]==c("samplespacing","imagesize"))	
	attr(l,"samplespacing")<-r['samplespacing'][[1]]
	attr(l,"imagesize")<-r['imagesize'][[1]]
	
	tracings=r[-c(1:2)]
	if(length(tracings)==0) stop("No tracings in this file")
	if(Verbose) cat("There are",length(tracings),"tracings in this file\n")	
	
	for(i in 1:length(tracings)){
		l[[i]]=xmlSApply(tracings[[i]],function(x) as.numeric(xmlAttrs(x)[c("xd","yd","zd")]))
		l[[i]]=t(l[[i]])
		rownames(l[[i]])<-NULL
		colnames(l[[i]])<-c("X","Y","Z")
		attr(l[[i]],'pathAttributes')=xmlAttrs(tracings[[i]])
	}
	l
}

ReadNeuronsFromLongairTraces<-function(f,...){
	l=ReadLongairTraceData(f,...)
	dflist=list()
	for(i in seq(l)){
		d=l[[i]]
		df=data.frame(PointNo=1:nrow(d),Label=2)
		df=cbind(df,d)
		df$radius=1
		df$Parent=df[,1]-1
		pathAttributes=attr(l[[i]],"pathAttributes")
		
		if('startson'%in%names(pathAttributes)){
			# this path is actually joined to another
			# nb Mark's paths are 0 indexed (R is 1 indexed)
			parentStartIndex=as.numeric(pathAttributes['startsindex'])+1
			# set the parent of this new path to the point on parent path
			df$Parent[1]=parentStartIndex
			# adjust all other point ids to start by the number of rows
			# in the parent path data frame
			parentPathId=as.numeric(pathAttributes['startson'])+1			
			df$Parent[-1]=df$Parent[-1]+nrow(dflist[[parentPathId]])
			df$PointNo=df$PointNo+nrow(dflist[[parentPathId]])
			# now add these data to the dataframe for the parent path
			dflist[[parentPathId]]=rbind(dflist[[parentPathId]],df)
			dflist[[i]]=NA
		} else {
			df$Parent[1]=-1
			dflist[[i]]=df
		}
	}
	# dflist
	neuronList=list()
	for(df in dflist){
		if(!is.data.frame(df)) next
		neuronList[[length(neuronList)+1]]=ParseSWCTree(df,f)
	}
	neuronList
}