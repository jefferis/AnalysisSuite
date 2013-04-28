
# Wrapper function to call plotneuron2d
# if given a set of neurons (names or numbers)
plotneurons2d<-function(NeuronRef,MultiPlot=F,Ask=!MultiPlot,LineColours=NULL,...){
    # If there are several neurons to plot, it makes sense to pause
    # unless we have multiplot situation when we want to superimpose
    # neurons
    oldpar<-par(ask=Ask)
    # the ... should allow any additional arguments to be passed to plotneuron2d
    if(length(LineColours)!=length(NeuronRef)) LineColours=rep(LineColours,length(NeuronRef))
    if(length(NeuronRef)>1 && MultiPlot){
	plotneuron2d(NeuronRef[1],LineColour=LineColours[1],...)
 	for(i in 2:length(NeuronRef))
	    plotneuron2d(NeuronRef[i],Superimpose=T,LineColour=LineColours[i],...)
    } else for(i in 1:length(NeuronRef)){
	plotneuron2d(NeuronRef[i],LineColour=LineColours[i],...)
    }
    par(oldpar)
}

####################
#                  #
#   plotneuron2d   #
#                  #
####################

plotneuron2d<-function(ANeuron,WithLine=T,NodesOnly=T,EPsOnly=F,BPsOnly=F,WithText=F,NewWindow=F,
		ScalePlot=F,  # ScalePlot is irrelevant now that I have found par(asp=1)
		UseCurPalette=F,WithContour=T, WithContours=F, WithLHEPoint=T,PlotAxes="XY",axes=TRUE,asp=1,
		MainTitle=paste(ANeuron$NeuronName,ANeuron$CellType),RSOrders=F,xlim=NULL,ylim=NULL,AxisDirections=c(1,-1,1),
		Superimpose=F,LineColour=NULL,PointAlpha=1,tck=NA,lwd=par("lwd"),...){
		# NodesOnly=T will plot dots at endpoints and branchpoints
		# BPsOnly will override this and plot only branchpoints
		# PointAlpha specifies alpha value to use for all points
		
		
		
		#    WithLine<-NodesOnly<-WithLHEPoint<-ScalePlot<-WithContour<-T
		#    WithText<-NewWindow<-UseCurPalette<-WithContours<-F;WithLHEPoint<-T
		# PlotAxes<-c("X","Y")
		if (is.character(ANeuron)){
				ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
		}
		if (is.numeric(ANeuron)){
				ANeuron<-MyNeurons[[ANeuron]]
		}
		
		if (!is.list(ANeuron)){
				warning("Cannot understand passed neuron")
				return(F)
		}
		
		# ImageJ etc use the top left hand corner as the origin
		# whereas R uses cartesian co-ords
		# so default is to invert Y
		# Simplest just to multiply
		# data (contours and points by 1 or -1)
		# at the outset
		if(any(AxisDirections!=1)){
				ANeuron$d[,c("X","Y","Z")]=t(t(ANeuron$d[,c("X","Y","Z")])*AxisDirections)
		}
		
		if(PlotAxes=="XY") {PlotAxes<-c("X","Y");NumPlotAxes<-c(1,2)} else
		if(PlotAxes=="YZ") {PlotAxes<-c("Y","Z");NumPlotAxes<-c(2,3)} else
		if(PlotAxes=="XZ") {PlotAxes<-c("X","Z");NumPlotAxes<-c(1,3)} else 
		if(PlotAxes=="ZY") {PlotAxes<-c("Z","Y");NumPlotAxes<-c(3,2)}
		
		OldPalette<-palette()
		if(!UseCurPalette){
				palette(c("black",rainbow(6)))
		}
		
		if (NewWindow) macintosh(w=9.5,h=9.5)
		
		# This is so that the contour points don't fall off the plot
		if((WithContours || WithContour) ){  # Removed  & !is.null(ANeuron$c)
				myxlims<-range(c(ANeuron$d[PlotAxes[1]],ANeuron$c$d[PlotAxes[1]],
								ANeuron$LH$d[PlotAxes[1]],ANeuron$MB$d[PlotAxes[1]]),na.rm = TRUE)
				myylims<-range(c(ANeuron$d[PlotAxes[2]],ANeuron$c$d[PlotAxes[2]],
								ANeuron$LH$d[PlotAxes[2]],ANeuron$MB$d[PlotAxes[2]]),na.rm = TRUE)
		} else {
				myxlims<-range(ANeuron$d[PlotAxes[1]],na.rm = TRUE)
				myylims<-range(ANeuron$d[PlotAxes[2]],na.rm = TRUE)
		}
		
		# if xlim and ylim were set, just use them
		if (!is.null(xlim)){
				myxlims=xlim
		}
		if (!is.null(ylim)){
				myylims=ylim
		}    
		
		if(ScalePlot){
				deltax<-diff(myxlims)
				deltay<-diff(myylims)
				if (deltax>deltay){
						myylims[1]<-myylims[1]-(deltax-deltay)/2
						myylims[2]<-myylims[2]+(deltax-deltay)/2
				} else {
						myxlims[1]<-myxlims[1]-(deltay-deltax)/2
						myxlims[2]<-myxlims[2]+(deltay-deltax)/2
				}
		}
		if(EPsOnly){
				NodesOnly=setdiff(ANeuron$EndPoints,ANeuron$StartPoint)
				# NB Remove start Point  - not really sure in the end whether this is 
				# always treated as an end point or not
				mycols<-rep(rgb(0,1,0,PointAlpha),length(ANeuron$EndPoints))
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else if(BPsOnly){
				NodesOnly<-ANeuron$BranchPoints
				mycols<-rep("red",length(ANeuron$BranchPoints))
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else if(NodesOnly){
				NodesOnly<-c(ANeuron$BranchPoints,ANeuron$EndPoints,ANeuron$StartPoint)
				mycols<-c(rep("red",length(ANeuron$BranchPoints)),
						rep("green",length(ANeuron$EndPoints)),"purple" )
				PlottedPoints<-ANeuron$d[NodesOnly,c("PointNo",PlotAxes)]
		} else { # end if(NodesOnly)
				mycols<-rep("black",ANeuron$NumPoints)
				mycols[ANeuron$BranchPoints]<-"red"
				mycols[ANeuron$EndPoints]<-"green"
				mycols[ANeuron$StartPoint]<-"purple"
				PlottedPoints<-ANeuron$d[,c("PointNo",PlotAxes)]
		} # end if(BPsOnly)
		
		
		# Add the LHAnchorPoint that is the "major" branching point
		# if it exists
		if(!is.null(ANeuron$LHAnchorPoint)){
				PlottedPoints<-rbind(PlottedPoints,
						ANeuron$d[ANeuron$LHAnchorPoint,c("PointNo",PlotAxes)])
				mycols<-c(mycols,"blue")
		}
		
		# Add the AxonLHEP that is the point on the axon
		# closest to the lateral horn entry point defined
		# by Takaki
		if(!is.null(ANeuron$AxonLHEP)){
				PlottedPoints<-rbind(PlottedPoints,
						ANeuron$d[ANeuron$AxonLHEP,c("PointNo",PlotAxes)])
				mycols<-c(mycols,"purple")
		}
		
		#NOW DO THE PLOT
		#cat("myylims =",myylims,"\n")
		if(Superimpose) points(PlottedPoints[,PlotAxes],col=mycols,pch=20,asp=asp,...) 
		else plot(PlottedPoints[,PlotAxes],col=mycols,pch=20,xlim=myxlims,ylim=myylims,
	#		main=MainTitle,asp=1,axes=all(AxisDirections==1),tck=tck,...) 
			main=MainTitle,asp=asp,axes=axes && all(AxisDirections==1),tck=tck,...) 
#		if(!all(AxisDirections==1) && !Superimpose){
		if(axes && !all(AxisDirections==1) && !Superimpose){
				# need to provide special treatment for axes
				box()
				if(AxisDirections[NumPlotAxes][1]!=1){
						axis(1, at=axTicks(1),tck=tck,
								label=axTicks(1)*AxisDirections[NumPlotAxes][1])
				} else {
						axis(tck=tck,1)
				}
				if(AxisDirections[NumPlotAxes][2]!=1){
						axis(2, at=axTicks(2),tck=tck,
								label=axTicks(2)*AxisDirections[NumPlotAxes][2])
				} else {
						axis(tck=tck,2)
				}
		}
		
		if(WithText){
				text(PlottedPoints[,PlotAxes[1]],PlottedPoints[,PlotAxes[2]],
						PlottedPoints[,"PointNo"],col=mycols,pos=3)
		}
		
		if (WithLine){
				# All Black lines
				if(!is.null(LineColour) && length(LineColour==1)){
						MyCols=rep(LineColour,ANeuron$NumSegs)
				} else {
						
						MyCols<-rep(1,ANeuron$NumSegs)
						# Coloured by Orders from root
						if(!is.null(ANeuron$SegOrders)){
								MyCols<-ANeuron$SegOrders
						}
						# Coloured by Segment Types (Axon,LH,MB)
						if(!is.null(ANeuron$SegTypes)){
								MyCols<-ANeuron$SegTypes
						}
						# Coloured according to the reverse Strahler procedure
						if(RSOrders && !(is.null(ANeuron$RSOrders))){
								MyCols[ANeuron$LHSegNos]<-ANeuron$RSOrders+1
						}
						if(!is.null(LineColour)){
							MyCols=LineColour[MyCols]
						}
						
				}
				# NOW PLOT THE LINES!	    
				#for(i in 1:ANeuron$NumSegs){
				#    lines(ANeuron$d[ANeuron$SegList[[i]],PlotAxes]
				#	,col=MyCols[i])
				#}  
				############## 
				# Found a quicker way to do this using NAs
				# to interrupt plotting and start a new line
				# Unfortunately have to do this for each new line colour
				# since lines can only do one line colour at a time
				for(thisCol in unique(MyCols)){
						SegsToPlot=ANeuron$SegList[MyCols==thisCol]
						LinesToPlot=unlist(sapply(SegsToPlot,function(x) c(x,NA)))	
						lines(ANeuron$d[LinesToPlot,PlotAxes],col=thisCol,lwd=lwd)
				}
		} # end if (WithLine) 
		
		
		
		#See if any contours exist
		for(ContSet in list(ANeuron[c("c","LH","MB")])){
				if(is.null(ContSet$d)) next  # Seems to be optional, but better safe
				#Plot all contours if required
				if(WithContours){
						
						ContsToPlot=sort(unique(ContSet$d$ContourID))
						PointsToPlot=rep(0,nrow(ContSet$d)+length(ContsToPlot))
						PointsToPlot=NULL
						for(j in ContsToPlot){
								#ThisContourPoints<-which(ContSet$d$ContourID==j)
								# polygon will join the points (only works in 2D)
								#polygon(ContSet$d[ThisContourPoints,PlotAxes],
								#    type="l",lty="dotted",col="black")
								#
								PointsToPlot=c(PointsToPlot,which(ContSet$d$ContourID==j),NA)
						}
						#TRY A SPEEDIER VERSION
						polygon(ContSet$d[PointsToPlot,PlotAxes],lty="dotted",type='l')
						
				}#endif(WithContours){
				# Plot a boundary contour (2D convex hull)
				if(WithContour){
						polygon(ContSet$d[chull(ContSet$d[,PlotAxes]),PlotAxes])
				}
		}#end for
		
		if(WithLHEPoint & !is.null(ANeuron$MarkerPoints$LHEPoint)){
				points(ANeuron$MarkerPoints$LHEPoint[NumPlotAxes[1]],
						ANeuron$MarkerPoints$LHEPoint[NumPlotAxes[2]],col="purple",pch=22)
		}
		if(!is.null(ANeuron$c$GrandCent)){
				points(ANeuron$c$GrandCent[NumPlotAxes[1]],
						ANeuron$c$GrandCent[NumPlotAxes[2]],col="red",bg="red",pch=22)
		}
		
		palette(OldPalette)
		#nb makes an invisible return - doesn't print on command line
		#but can be stored in a variable.
		invisible(PlottedPoints)
}
