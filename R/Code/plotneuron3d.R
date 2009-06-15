plotneuron3d.simple<-function(ANeuron, WithLine=T,NeuronNames=FALSE,
	WithNodes=T,WithAllPoints=F,WithText=F,HighlightLongestSegment=FALSE,PlotSubTrees=T,ClearRGL=T,NeuronList=MyNeurons,col,...){
	# rewrite of plotneuron3d using updated rgl calls
	
	require(rgl)
	if(ClearRGL) rgl.clear()
    
	if(missing(col)) col='green'

	if ( is.character(ANeuron) || is.numeric(ANeuron) ){
		if(length(ANeuron)>1){
			
			if(is.function(col)) col=col(length(ANeuron))
			if(is.factor(col)) col=rainbow(nlevels(col))[as.integer(col)]
			if(is.logical(NeuronNames) && NeuronNames) NeuronNames=ANeuron
			return(invisible(mapply(plotneuron3d.simple,
			NeuronList[ANeuron],col=col,ClearRGL=FALSE,NeuronNames=NeuronNames,
			WithNodes=WithNodes,WithText=WithText,WithLine=WithLine,
			HighlightLongestSegment=HighlightLongestSegment,PlotSubTrees=PlotSubTrees,...)))
		} else ANeuron<-NeuronList[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
		warning("Cannot understand passed neuron")
		return(F)
    }
	# at this point we've only got one neuron, but we could still have col as a function
	if(is.function(col)) col=col(1)
	
	rglreturnlist=list()
	NodesOnly<-c(ANeuron$BranchPoints,
	    ANeuron$EndPoints[-which(ANeuron$EndPoints==ANeuron$StartPoint)],
	    ANeuron$StartPoint)
	NodeCols<-c(rep("red",length(ANeuron$BranchPoints)),
	    rep("green",length(ANeuron$EndPoints)-1),"purple" )
	
	if(WithNodes){
		Colour=col
		if(!WithLine) NodeCols=rep(Colour,length(NodeCols))
		rglreturnlist[["points"]]=points3d(ANeuron$d[NodesOnly,c("X","Y","Z")],color=NodeCols,size=3)
		if(WithText) # text labels for nodes
		rglreturnlist[["texts"]]=texts3d(ANeuron$d[NodesOnly,c("X","Y","Z")],text=NodesOnly,color=NodeCols,adj=c(0,0.5))
	}

	if(PlotSubTrees && !is.null(ANeuron$SubTrees)) ANeuron$SegList=unlist(ANeuron$SubTrees,recursive=FALSE)
	
	d=data.matrix(ANeuron$d[,c("X","Y","Z")])
	# xyzl=sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)})
	# NAs are used to break line segments
	xyzl<-do.call(rbind,sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)},simplify=FALSE))
	rglreturnlist[["lines"]]=lines3d(xyzl,col=col,...)
	if(HighlightLongestSegment){
		x=GetLongestSegment(ANeuron)
		rglreturnlist[["spheres"]]=spheres3d(x,col=col,...)
	}
	if(is.logical(NeuronNames) && NeuronNames) NeuronNames=ANeuron$NeuronName
	if(!is.logical(NeuronNames)){
		StartPoint=ifelse(is.null(ANeuron$StartPoint),1,ANeuron$StartPoint)
		texts3d(d[StartPoint,],texts=NeuronNames,col=col)
	}
	
	invisible(rglreturnlist)
}

IdentifyNeuronsMatchingRegion<-function(neurons,select3dfunc,EndPointsOnly=TRUE,Verbose=FALSE,...){
	## use a qeury function from rgl's select3d to outline a query region
	## in current rgl scene and see if any neurons are inside
	if(missing(select3dfunc)) {
		cat("Please choose a region using the left mouse button\n")
		select3dfunc=select3d(...)
	}
	l=list()
	for (n in neurons){
		if(EndPointsOnly)
			querypoints=n$d[n$EndPoints,c("X","Y","Z")]
		else querypoints=n$d[,c("X","Y","Z")]
		
		if(any(result<-select3dfunc(querypoints))){
			if(EndPointsOnly) positiveids=n$EndPoints[result]
			else positiveids=which(result)
			l[[n$NeuronName]]=positiveids
			if(Verbose){
				cat("neuron",n$NeuronName,"matches the selection region\n")
				cat("matching points =",positiveids,"\n");
			}
		}
	}
	if(length(l)==0) l=NULL
	invisible(l)
}

animate3d<-function(time=10,degrees=360,filestem=NULL){
	degreesPersecond=degrees/time
	start=proc.time()[3]
	while ((i <- degreesPersecond*(proc.time()[3]-start)) < degrees) {
		rot3d(y=i)
		if(!is.null(filestem))
			snapshot3d(paste(filestem,sep="",sprintf("%04.f",i),".png"))
	}	
}

set3d<-function(pos=c("front","left","back","right","ventral","dorsal"),zoom=.7,...){
	pos=match.arg(pos)
	m=diag(c(1,-1,-1,1)) # front

	if(pos=="left") {
		m=diag(c(0,-1,0,1))
		m[1,3]=1;m[3,1]=1
	}
	if(pos=="back") {
		m=diag(c(-1,-1,1,1))
	}
	if(pos=="right") {
		m=diag(c(0,-1,0,1))
		m[1,3]=m[3,1]=-1
	}
	view3d(userMatrix=m,zoom=zoom,...)
}
rot3d<-function(xangle=0,yangle=0,zangle=0,initialUserMatrix=get("r3dDefaults", envir=.GlobalEnv)$userMatrix){
	angle=max(xangle,yangle,zangle)
	if(angle!=0) {
		x=xangle/angle
		y=yangle/angle
		z=zangle/angle
	}	
	m=initialUserMatrix%*%rotationMatrix(x=x,y=y,z=z,angle=angle*2*pi/360)
	view3d(userMatrix=m)
	invisible(m)
}

a3d=animate3d
s3d=set3d
r3d=rot3d
