# InterpolateAlongNeuron.R

#RELEASE
#BEGINCOPYRIGHT
###############
# R Source Code to accompany the manuscript
#
# "Comprehensive Maps of Drosophila Higher Olfactory Centers: 
# Spatially Segregated Fruit and Pheromone Representation"
# Cell (2007), doi:10.1016/j.cell.2007.01.040
# by Gregory S.X.E. Jefferis*, Christopher J. Potter*
# Alexander M. Chan, Elizabeth C. Marin
# Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo
#
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# See flybrain.stanford.edu for further details
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################
#ENDMAINCOPYRIGHT

# Need to make a function that can take a neuron 
# and interpolate new points at a given step size along each line segment.
# 
# 2005-07-24: noticed that this version does not deal with branch points
# - it actually duplicates each branchpoint the second (and third) times 
# it encounters it.
# In general this is part of failure to pay any attention to the structure
# of the neuron.  Probably need to keep a table of the old and new branch
# points

# do this as follows
# For each head of a segment, ask if we have an entry in our table.
# If yes, then give the point in the seg list that number
# else increment counter, and add old/newID pair
# Then add as many internalpoints as required
# then add the tail (which should never have appeared before)
# and add an old/newID pair

# TO DO: interpolate width as well
#        calculate parents for new points

# source(file.path(CodeDir,"InterpolateAlongNeuron.R"))
InterpolateAlongNeuron<-function(ANeuron,stepSize=0.5){
		# Scheme
		d=matrix(unlist(ANeuron$d[,c("X","Y","Z")]),ncol=3)

		# calculate seglengths if we haven't 
		if(is.null(ANeuron$SegLengths)){
				warning(paste("Calculating SegLengths for",ANeuron$NeuronName))
				ANeuron$SegLengths=SegLengths(ANeuron)
		}
		
		
		oldID=NULL; newID=NULL
		
		newseglist=ANeuron$SegList
		
		totalPoints=sum(sapply(ANeuron$SegLengths,function(x) 2+floor((x-1e-9)/stepSize)))
		pointsSoFar=0
		pointArray=matrix(0,ncol=3,nrow=totalPoints)
		for(i in seq(len=length(ANeuron$SegList))){
				
				# length in microns of this segment
				l=ANeuron$SegLengths[i]
				
				if(l>stepSize){
						# new internal points, measured in length along segment
						internalPoints=seq(stepSize,l,by=stepSize)
						nInternalPoints=length(internalPoints)
						# if the last generated one is actually in exactly the same place 
						# as the endpoint then discard it
						if(internalPoints[nInternalPoints]==l) internalPoints=internalPoints[-length(internalPoints)]
						
						# find lengths between each original point on the segment
						diffs=diff(d[ANeuron$SegList[[i]],])
						indSegLens=sqrt(rowSums(diffs*diffs))
						cs=c(0,cumsum(indSegLens))
						
						#idxs=sapply(internalPoints,function(x) max(which(x>=cs)))
						#startPos=0;nextPos=internalPoints[1]
						idxs=rep(0,length(internalPoints))
						for(j in seq(len=length(cs))){
								idxs[idxs==0 & internalPoints<cs[j]]=j-1
						}
								
						newPoints=matrix(0,ncol=3,nrow=nInternalPoints+2)
						newPoints[1,]= d[ANeuron$SegList[[i]][1],]

						froms=d[ANeuron$SegList[[i]][idxs],]
						deltas=diffs[idxs,]
						fracs=(internalPoints-cs[idxs])/indSegLens[idxs]
						newPoints[-c(1,nrow(newPoints)),]=froms+(deltas*fracs)
						nNewPoints=nrow(newPoints)
						newPoints[nNewPoints,]=d[ANeuron$SegList[[i]][length(ANeuron$SegList[[i]])],]
				} else {
						nNewPoints=2
						newPoints=d[ANeuron$SegList[[i]][c(1,length(ANeuron$SegList[[i]]))],]
				}
				
				newseg=NULL
				# have we seen the headpoint of this seg before?
				if(any(ANeuron$SegList[[i]][1]==oldID)){
						# yes 
						nNewPoints=nNewPoints-1
						newPoints=newPoints[-1,] # prevent this head from being readded to point array
						newseg=c(newID[ANeuron$SegList[[i]][1]==oldID],
								seq(from=pointsSoFar+1,by=1,len=nNewPoints))
				} else {
						# no, make a note of it and add it to the array
						oldID=c(oldID,ANeuron$SegList[[i]][1])
						newID[length(oldID)]=pointsSoFar+1
						newseg=	seq(from=pointsSoFar+1,by=1,len=nNewPoints)
				}
				
				# add the tail to the table we are keeping track of
 				oldID=c(oldID,ANeuron$SegList[[i]][length(ANeuron$SegList[[i]])])
				newID=c(newID,pointsSoFar+nNewPoints)
				
				newseglist[[i]]=newseg
				pointArray[(pointsSoFar+1):(pointsSoFar+nNewPoints),]=newPoints
				#cat("i = ",i," pointsSoFar=",pointsSoFar,"nNewPoints=",nNewPoints,"\n")
				pointsSoFar=pointsSoFar+nNewPoints
		}
		pointArray=pointArray[1:pointsSoFar,]
		colnames(pointArray)=c("X","Y","Z")

		#return(oldID,newID,pointArray,newseglist)
		# OK now return a new neuron
		ANeuron$NumPoints=pointsSoFar
		ANeuron$StartPoint=newID[oldID==ANeuron$StartPoint]
		ANeuron$BranchPoints=newID[match(ANeuron$BranchPoints,oldID)]
		ANeuron$EndPoints=newID[match(ANeuron$EndPoints,oldID)]
		if(any(is.na(c(ANeuron$EndPoints,ANeuron$BranchPoints)))){
				stop("Problem matching up old & new end/branchpoints")
		}
							
		ANeuron$SegList=newseglist
		ANeuron$d=data.frame(PointNo=1:pointsSoFar,X=pointArray[,1],Y=pointArray[,2],Z=pointArray[,3])
		
		return(ANeuron)
}

NurbSmoothNeuron<-function(ANeuron, NurbKnotSpacing=5, NurbKnots=NA, ... ){
	
	# use a smoothing nurbs curve to interpolate individual neurons
	d=do.call(cbind,ANeuron$d[,c("X","Y","Z")])

	# calculate seglengths if we haven't 
	if(is.null(ANeuron$SegLengths)){
		warning(paste("Calculating SegLengths for",ANeuron$NeuronName))
		ANeuron$SegLengths=SegLengths(ANeuron)
	}
	segsWithMoreThan4Points=which(sapply(ANeuron$SegList,length)>3)
	if(is.na(NurbKnots)){
		# Linearly interpolate along neuron using new (coarser) spacing
		interpNeuron=InterpolateAlongNeuron(ANeuron,stepSize=NurbKnotSpacing)
		di=matrix(unlist(interpNeuron$d[,c("X","Y","Z")]),ncol=3)
	}
	for(i in segsWithMoreThan4Points){
		l=d[ANeuron$SegList[[i]],]
		if(!is.na(NurbKnots)){
			# Divide each segment into a fixed number of knots
			linterp=DivideLineIntoNEqualSubLines(l,NurbKnots)
		} else {
			# just skip this segment if we can't get at least 3 knots in
			if(ANeuron$SegLengths[i]<(NurbKnotSpacing*2)) next
			linterp=di[interpNeuron$SegList[[i]],]
			# drop any duplicated points (this can happen in Longair tracings)
			linterp=linterp[!duplicated(linterp),]
		}
		# note first and last points will be the same
		replacementMatrix=jgeom.nurbsInterpolate(linterp,seq(from=0,to=1,len=nrow(l)))
		d[ANeuron$SegList[[i]],]=replacementMatrix
	}
	ANeuron$d[,c("X","Y","Z")]=d
	ANeuron
}

MovingAverageSmoothNeuron<-function(ANeuron, filter=c(1/4,1/2,1/4)){
	# do a weighted moving average smoothing

	pointsToAverage=length(filter)
	# normalise filter if required
	if(sum(filter)!=1){
#		warning("Normalising filter in MovingAverageSmoothNeuron")
		filter=filter/sum(filter)
	}
	
	d=do.call(cbind,ANeuron$d[,c("X","Y","Z")])

	# calculate seglengths if we haven't 
	if(is.null(ANeuron$SegLengths)){
		warning(paste("Calculating SegLengths for",ANeuron$NeuronName))
		ANeuron$SegLengths=SegLengths(ANeuron)
	}
	segsWithMoreThanNPoints=which(sapply(ANeuron$SegList,length)>pointsToAverage)

	for(i in segsWithMoreThanNPoints){
		l=d[ANeuron$SegList[[i]],]
		nrowl=nrow(l)
		
		fxyz=apply(l,2,filter,filter=filter)
		pointsToChange=!is.na(fxyz[,1])
		
		# note some of the first and last points will be kept
		d[ANeuron$SegList[[i]][pointsToChange],]=fxyz[pointsToChange,]
	}
	ANeuron$d[,c("X","Y","Z")]=d
	# recalculate seglenths (since these have now changed)
	ANeuron$SegLengths=SegLengths(ANeuron)
	ANeuron
}

SplineSmoothNeuron<-function(ANeuron, ... ){
	# Will take a supplied neuron and do succesive smoothing splines for Y and Z coords
	# ... will be passed to smooth.spline
	d=matrix(unlist(ANeuron$d[,c("X","Y","Z")]),ncol=3)

	# calculate seglengths if we haven't 
	if(is.null(ANeuron$SegLengths)){
			warning(paste("Calculating SegLengths for",ANeuron$NeuronName))
			ANeuron$SegLengths=SegLengths(ANeuron)
	}
	segsWithMoreThan4Points=which(sapply(ANeuron$SegList,length)>3)
	for(i in segsWithMoreThan4Points){
		
		l=d[ANeuron$SegList[[i]],]
		l10=DivideLineIntoNEqualSubLines(l,10)
		# Weights for points - set first and last to very high to keep the same
		w=rep(1,nrow(l))
		w[1]=1e6
		w[length(w)]=w[1]
		# for(j in seq(from=(n-1),len=nrow(newPoints)-(n-1))){
		# 	pointsToAverage=oldPoints[seq(n)+j-1]
		# 	newPoints[j,]
		# }
		# Sort the 3 axes according to which has the largest total vector length
		# this helps prevent nasty spline behaviour
		# axisOrders=order(colSums(sqrt(diff(ll)^2)),decreasing=TRUE)
		# axisOrders=3:1
		chosen=list()
		chosen$combinedFitScore=Inf
		for(currentAxis in 1:3){
			otherAxes=(1:3)[-currentAxis]
			t=try(lxy<-smooth.spline(l[,c(currentAxis,otherAxes[1])],w=w,...), silent=TRUE)
			if(inherits(t,'try-error')) next
			t=try(lxz<-smooth.spline(l[,c(currentAxis,otherAxes[2])],w=w,...), silent=TRUE)
			if(inherits(t,'try-error')) next
			
			cat(i,":",lxy$crit,lxz$crit,"\n")
			# str(l[,1]);str(lxy$y);str(lxz$y)
			
			# Assess goodness of fit
			newMatrix=cbind(l10[,currentAxis],predict(lxy,l10[,currentAxis])$y,predict(lxz,l10[,currentAxis])$y)
			axisOrders<-c(currentAxis,otherAxes)
			combinedFitScore=sum(abs(newMatrix-l10[,axisOrders]))			
			# combinedFitScore=sum(abs(c(lxy$crit,lxz$crit)))
			cat("combinedFitScore =",combinedFitScore,"\n")
			if(combinedFitScore<chosen$combinedFitScore){
				chosen$combinedFitScore=combinedFitScore
				chosen$lxy=lxy
				chosen$lxz=lxz
				chosen$currentAxis=currentAxis
				chosen$otherAxes=otherAxes
			}
		}
		if(is.null(chosen$otherAxes)) stop(paste("Unable to find axis suitable for interpolation in segment",i))
		replacementMatrix<-with(chosen,cbind(l[,currentAxis],predict(lxy,l[,currentAxis])$y,predict(lxz,l[,currentAxis])$y))
		axisOrders<-with(chosen,c(currentAxis,otherAxes))
		d[ANeuron$SegList[[i]],axisOrders]=replacementMatrix
	}
	ANeuron$d[,c("X","Y","Z")]=d
	ANeuron
}

DivideLineIntoNEqualSubLines<-function(xyz,n){
	# Expects a matrix of p points with 3 columns and p rows
	# giving the 3d position of a series of line segments
	# NB a line consisting of 1 segment will have 2 points
	
	# this function will interpolate this into n segments of equal length (i.e. n+1 points)

	xyz=data.matrix(xyz)
	l=seglength(xyz)
	
	if(n<1){
		warning("Can't divide a line into less than 1 segment")
		return(null)
	} else if(n==1){
		# Just take the start and end points 
		return(xyz[c(1,nrow(xyz)),])
	}
	# ... else there is more than 1 segment to consider
	# new internal points, measured in length along segment
	
	stepSize=l/n
	
	internalPoints=seq(from=stepSize,by=stepSize,length.out=n-1)
	nInternalPoints=length(internalPoints)
	# if the last generated one is actually in exactly the same place 
	# as the endpoint then discard it

	# find lengths between each original point on the segment
	diffs=diff(xyz)
	indSegLens=sqrt(rowSums(diffs*diffs))
	cs=c(0,cumsum(indSegLens))
	
	# find the indices of the old points that will contribute to the new sublines
	idxs=rep(0,length(internalPoints))
	for(j in seq(len=length(cs))){
		idxs[idxs==0 & internalPoints<cs[j]]=j-1
	}	
	
	# Make tha array to store the results
	newPoints=matrix(0,ncol=3,nrow=n+1)
	
	# The first point is the same as the original first point 
	newPoints[1,]= xyz[1,]
	# and the last point is the same as the original last point  
	newPoints[nrow(newPoints),]=xyz[nrow(xyz),]

	# Now interpolate along the original sequence of sublines to find the positions of the new points
	froms=xyz[idxs,]
	deltas=diffs[idxs,]
	fracs=(internalPoints-cs[idxs])/indSegLens[idxs]
	newPoints[-c(1,nrow(newPoints)),]=froms+(deltas*fracs)
	return(newPoints)
}
