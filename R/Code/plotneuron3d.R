plotneuron3d.simple<-function(ANeuron, WithLine=T,
	WithNodes=T,WithAllPoints=F,WithText=F,HighlightLongestSegment=FALSE,PlotSubTrees=T,ClearRGL=T,NeuronList=MyNeurons,col,...){
	# rewrite of plotneuron3d using updated rgl calls
	
	require(rgl)
	if(ClearRGL) rgl.clear()
    
	if(missing(col)) col='green'

	if ( is.character(ANeuron) || is.numeric(ANeuron) ){
		if(length(ANeuron)>1){
			if(is.function(col)) col=col(length(ANeuron))
			if(is.factor(col)) col=rainbow(nlevels(col))[as.integer(col)]
			return(invisible(mapply(plotneuron3d.simple,
			NeuronList[ANeuron],col=col,ClearRGL=FALSE,
			WithNodes=WithNodes,WithText=WithText,WithLine=WithLine,
			HighlightLongestSegment=HighlightLongestSegment)))
		} else ANeuron<-NeuronList[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
		warning("Cannot understand passed neuron")
		return(F)
    }
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
	
	d=data.matrix(ANeuron$d[,c("X","Y","Z")])
	# xyzl=sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)})
	# NAs are used to break line segments
	xyzl<-do.call(rbind,sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)},simplify=FALSE))
	rglreturnlist[["lines"]]=lines3d(xyzl,col=col,...)
	if(HighlightLongestSegment){
		x=GetLongestSegment(ANeuron)
		rglreturnlist[["spheres"]]=spheres3d(x,col=col,...)
	}
	invisible(rglreturnlist)
}