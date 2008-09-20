plotneuron3d.simple<-function(ANeuron, WithLine=T,
	WithNodes=T,WithAllPoints=F,WithText=F,PlotSubTrees=T,ClearRGL=T,NeuronList=MyNeurons,...){
	# rewrite of plotneuron3d using updated rgl calls
	
	require(rgl)
	if(ClearRGL) rgl.clear()
    if (is.character(ANeuron)){
	ANeuron<-NeuronList[[GetNeuronNum(ANeuron)]]
    }
    if (is.numeric(ANeuron)){
	ANeuron<-NeuronList[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
	warning("Cannot understand passed neuron")
	return(F)
    }
	
	NodesOnly<-c(ANeuron$BranchPoints,
	    ANeuron$EndPoints[-which(ANeuron$EndPoints==ANeuron$StartPoint)],
	    ANeuron$StartPoint)
	NodeCols<-c(rep("red",length(ANeuron$BranchPoints)),
	    rep("green",length(ANeuron$EndPoints)-1),"purple" )
	
	if(WithNodes){
		if(!WithLine) NodeCols=rep(Colour,length(NodeCols))
		points3d(ANeuron$d[NodesOnly,c("X","Y","Z")],color=NodeCols,size=3)
		if(WithText) # text labels for nodes
		rgl.texts(ANeuron$d[NodesOnly,c("X","Y","Z")],NodesOnly,color=NodeCols,adj=c(0,0.5))		
	}
	
	d=data.matrix(ANeuron$d[,c("X","Y","Z")])
	# xyzl=sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)})
	# NAs are used to break line segments
	xyzl<-do.call(rbind,sapply(ANeuron$SegList,function(s) {rbind(d[s,],NA)},simplify=FALSE))
	lines3d(xyzl,...)
}