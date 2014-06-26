# Simplified version of original plotneuron3d function
# and related utility functions

# this means that objects of class neuron can be plotted by doing:
# plot3d(ANeuron)
# the plot3d generic is defined by rgl

plotneuron3d.simple<-function(ANeuron, WithLine=T,NeuronNames=FALSE,
	WithNodes=T,WithAllPoints=F,WithText=F,HighlightLongestSegment=FALSE,PlotSubTrees=T,
	ClearRGL=T,NeuronList=MyNeurons,col=NULL,...){
	
	.Deprecated('nat::plot3d.neuron')
	if(HighlightLongestSegment)
		stop("HighlightLongestSegment has been retired. See GetLongestSegment")
	plot3d(ANeuron)

	require(rgl)
	r3dDefaults$bg<-'grey'
	
	if ( is.character(ANeuron) || is.numeric(ANeuron) ){
		if(length(ANeuron)>1){
			if(missing(col)) col=rainbow
			if(is.function(col)) col=col(length(ANeuron))
			if(is.factor(col)) col=rainbow(nlevels(col))[as.integer(col)]
			if(is.logical(NeuronNames) && NeuronNames) NeuronNames=ANeuron
			return(invisible(mapply(plotneuron3d.simple,
			NeuronList[ANeuron],col=col,ClearRGL=FALSE,NeuronNames=NeuronNames,
			WithNodes=WithNodes,WithText=WithText,WithLine=WithLine,
			PlotSubTrees=PlotSubTrees,...)))
		} else ANeuron<-NeuronList[[ANeuron]]
	}
	
	if (!is.neuron(ANeuron)){
		warning("Cannot understand passed neuron")
		return(F)
	}
	invisible(nat::plot3d(ANeuron,col=col,add=!ClearRGL,NeuronNames=NeuronNames,
	WithNodes=WithNodes,WithText=WithText,WithLine=WithLine,
	PlotSubTrees=PlotSubTrees,...))
}

plot3dsurface<-function(material,d,VertexComponent="Vertices",col=rainbow,...){
	.Deprecated('nat::plot3d.hxsurf')
	invisible(plot3d(as.hxsurf(d),materials=material,col=col,...))
}

#' convert old style surface object into hxsurf
as.hxsurf<-function(x){
	if(!inherits(x,'hxsurf')){
		# check if this looks like a surface
		if(!is.list(x) || is.null(x$Vertices) || is.null(x$RegionList))
			stop("d does not look like an hxsurf object")
		class(x)<-c('hxsurf',class(x))
		x$Regions=x[x$RegionList]
		x[x$RegionList]=NULL
	}
	x
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

#' Set the 3d view using some standard Drosophila anatomical terms
#' 
#' @param pos The chosen view. Defaults to front (=anterior)
#' @param zoom Numeric zoom factor (default 0.7)
#' @seealso \code{\link[rgl]{view3d}}
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
	if(pos=="ventral") {
		m=diag(c(1,-1,-1,1))
	}
	if(pos=="dorsal") {
		m=diag(c(1,1,1,1))
	}
	view3d(userMatrix=m,zoom=zoom,...)
}

rot3d<-function(xangle=0,yangle=0,zangle=0,
	initialUserMatrix=get("r3dDefaults", envir=.GlobalEnv)$userMatrix,...){
	angle=max(xangle,yangle,zangle)
	if(angle!=0) {
		x=xangle/angle
		y=yangle/angle
		z=zangle/angle
	}	
	m=initialUserMatrix%*%rotationMatrix(x=x,y=y,z=z,angle=angle*2*pi/360)
	view3d(userMatrix=m,...)
	invisible(m)
}

a3d=animate3d
s3d=set3d
r3d=rot3d
