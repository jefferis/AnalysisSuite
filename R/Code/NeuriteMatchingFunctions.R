# Functions for matching 3d linesets / neurons
# source(file.path(DicksonRoot,"JaiScripts/NeuriteMatchingFunctions.R"))

require(RANN) # for finding nearest neighbour points
require(dtw) # for dynamic programming best matching of lines

MedianMinPointDistance<-function(n1,n2,processfun=GetLongestSegment,summaryfun=median,BothDirections=FALSE,BothDirectionsFun=mean,...){
	# this function takes 2 neurons (or 2 point lists)
	# and finds the nearest neighbours 
	if(is.null(processfun)){
		if(is.list(n1)) n1=n1$d[,c("X","Y","Z")]
		if(is.list(n2)) n2=n2$d[,c("X","Y","Z")]
	} else {
		n1=processfun(n1)
		n2=processfun(n2)		
	}
	d1=summaryfun(nn2(n1,n2,k=1,...)$nn.dists)
	if(!BothDirections) return(d1)
	d2=summaryfun(nn2(n2,n1,k=1,...)$nn.dists)
	return(BothDirectionsFun(c(d1,d2)))
}

VectorAndDistanceMatch<-function(n1,n2,processfun=GetLongestSegment,summaryfun=median,BothDirections=FALSE,BothDirectionsFun=mean,...){
	# TODO
	# this function takes 2 neurons (or 2 point lists)
	# finds the nearest neighbours for all points
	# and the vector cross products for matching points
	if(is.null(processfun)){
		if(is.list(n1)) n1=n1$d[,c("X","Y","Z")]
		if(is.list(n2)) n2=n2$d[,c("X","Y","Z")]
	} else {
		n1=processfun(n1)
		n2=processfun(n2)		
	}
	d1=summaryfun(nn2(n1,n2,k=1,...)$nn.dists)
	
	if(!BothDirections) return(d1)
	d2=summaryfun(nn2(n2,n1,k=1,...)$nn.dists)
	return(BothDirectionsFun(c(d1,d2)))
}

WeightedNNComparison<-function(a,b,normalise=FALSE,costfn=function(x) dnorm(x,sd=4),summaryfn=mean,...){
	# TODO
	# expects 2 sets of points
	# calculates a score which will be a weighted function of the distance
	# to emphasise close matches
	n<-nn2(a,b)
	mean(costfn(n$nn.dists,...))
}


BestContinuousAlignment<-function(a,b,missingCost,k=10,...){
	# use dtw package + a cost function consisting of distances for 
	# k=10 or 20 nearest neighbours and e.g. 100?s everywhere else

	# a=reference, b=query according to dtw terminology
	# use nn2 to find the nearest neighbours in b for all points of a
	abnn=nn2(b,a,k=k)
	if(missing(missingCost)) missingCost=max(abnn$nn.dists)
	dmat=matrix(missingCost,nrow=nrow(a),ncol=nrow(b))
	for (i in 1:nrow(a)){
		dmat[i,abnn$nn.idx[i,]]=abnn$nn.dists[i,]
	}
	# for dtw, x=query, y=reference
	# Element [i,j] of the local-distance matrix is understood as the distance between element x[i] and y[j]
	# The distance matrix has therefore n=length(x) rows and m=length(y) columns	
	dtw(dmat,...)
}

DtwToScore<-function(bca){
	# takes a dtw object (calculated with keep.internals=TRUE)
	# and computes a score based on the alignment indices and distances
	# the score should be additive based on the idea of 
	# the sum of the goodness of fits for each matched point 
}

CompareAllSegs<-function(n1,n2,...){
	ns1=length(n1$SegList)
	ns2=length(n2$SegList)
	d1=data.matrix(n1$d[,c("X","Y","Z")])
	d2=data.matrix(n2$d[,c("X","Y","Z")])
	smat=matrix(NA,nrow=ns1,ncol=ns2)
	l<<-list()
	for(s in 1:ns1){
		for(t in 1:ns2){
			bca<-try(BestContinuousAlignment(d1[n1$SegList[[s]],],d2[n2$SegList[[t]],],...))
			if(!inherits(bca, "try-error")){
				smat[s,t]=bca$normalizedDistance
				l[[paste(s,t)]]<<-bca				
			} else smat[s,t]=NA
		}
	}
	smat
}

NNBasedLinesetMatching<-function(n1,n2,...){
	# first find nearest neighbours in both directions

	# accept either neurons or just the point dataframes
	if(is.list(n1) & !is.data.frame(n1)) n1=data.matrix(n1$d[,c("X","Y","Z","Parent")])
	if(is.list(n2) & !is.data.frame(n2)) n2=data.matrix(n2$d[,c("X","Y","Z","Parent")])

	a=n1[,c("X","Y","Z")]
	b=n2[,c("X","Y","Z")]
	
	nnn1=nn2(a,b,k=1,...)
	#nnn2=nn2(b,a,k=1,...)
	
	idxArray=cbind(nnn1$nn.idx,seq(length(nnn1$nn.idx)))
	# Need to supply a set of pairs of points.
	# will use the parent of each chosen point.
	# if parent undefined, then ignore that point
	
	# Calculate the direction vectors
	dvs=findDirectionVectorsFromParents(n1,n2,idxArray,ReturnAllIndices=TRUE)
	if(length(attr(dvs,"badPoints"))>0) nnn1$nn.dists=nnn1$nn.dists[-attr(dvs,"badPoints")]
	# Calculate segment lengths
	l1.seglengths=normbyrow(dvs[,1:3])
	l2.seglengths=normbyrow(dvs[,4:6])
	# normalise the direction vectors
	dvs[,1:3]=dvs[,1:3]/l1.seglengths
	dvs[,4:6]=dvs[,4:6]/l2.seglengths
	# Calculate the point displacement term
	m1=l1.seglengths*(nnn1$nn.dists)^2

	# Calculate the line angle mismatch term	
	m2=l1.seglengths^3/6*(1-dotprod(dvs[,1:3],dvs[,4:6]))
	mismatch=sum(m1+m2)
	cat("sum m1:",sum(m1),"sum m2:",sum(m2),"\n")
	if(mismatch < -1e-6) {
		stop("Negative line mismatch score!")
	} else if(mismatch<0){
		mismatch=0
	}
	return(mismatch)	
}

WeightedNNBasedLinesetDistFun<-function(nndists,dotproducts,sd=3,...){
	summaryfun=function(x) 1-mean(sqrt(x),na.rm=T)
	sapply(sd,function(sd) summaryfun(dnorm(nndists,sd=sd)*dotproducts/dnorm(0,sd=sd)))
}

WeightedNNBasedLinesetDistFun.Sum<-function(nndists,dotproducts,sd=3,...){
	summaryfun=function(x) sum(sqrt(x),na.rm=T)
	sapply(sd,function(sd) summaryfun(dnorm(nndists,sd=sd)*dotproducts/dnorm(0,sd=sd)))
}

WeightedNNBasedLinesetDistFun.Score<-function(nndists,dotproducts,sd=3,threshold=2*sd,...){
	# another distance function, this time with a notion of a threshold
	# typically 2 sigma.  If you are further away than this, then you get a negative score
	# most points will have a small negative score, so the expected score is negative
	summaryfun=function(x) sum(x,na.rm=T)
	# (dnorm(seq(0,10,0.1),sd=3)-dnorm(6,sd=3))/(dnorm(0,sd=3)-dnorm(6,sd=3))
	# transformed threshold
	sdt=dnorm(threshold,sd=sd)
	sapply(sd,function(sd) summaryfun( (dnorm(nndists,sd=sd)-sdt) /(dnorm(0,sd=sd)-sdt) * dotproducts))
}

WeightedNNBasedLinesetDistFun.separate<-function(nndists,dotproducts,sd=3){
	summaryfun=function(x) 1-mean(x)
	c(summaryfun(dnorm(nndists,sd=sd)/dnorm(0,sd=sd)),summaryfun(dotproducts))
}

WeightedNNBasedLinesetMatching<-function(n1,n2,dvs1=NULL,dvs2=NULL,
	NNDistFun=WeightedNNBasedLinesetDistFun,Verbose=FALSE,
	BothDirections=FALSE,BothDirectionsFun=list,OnlyClosestPoints=FALSE,...){
	# my hybrid version
	# returns a score based on the similarity of nearest neighbour location
	# and the dot product of the direction vectors

	# accept either neurons or just the point dataframes
	if(is.list(n1) & !is.data.frame(n1)) 
		n1=data.matrix(n1$d[,c("X","Y","Z","Parent")])
	if(is.list(n2) & !is.data.frame(n2))
		n2=data.matrix(n2$d[,c("X","Y","Z","Parent")])
	
	NNDistFun=match.fun(NNDistFun)
	BothDirectionsFun=match.fun(BothDirectionsFun)
	if(BothDirections){
		f=WeightedNNBasedLinesetMatching(n1,n2,dvs1,dvs2,NNDistFun=NNDistFun,Verbose=Verbose,BothDirections=FALSE,...)
		b=WeightedNNBasedLinesetMatching(n2,n1,dvs1,dvs2,NNDistFun=NNDistFun,Verbose=Verbose,BothDirections=FALSE,...)
		if(length(f)==1 && length(b)==1) return (BothDirectionsFun(f,b))
		if(length(dim(f))==1 && length(f)==length(b)) return (cbind(f,b))
		return(BothDirectionsFun(f,b))
	}
	
	a=n1[,c("X","Y","Z")]
	b=n2[,c("X","Y","Z")]
	
	nnn1=nn2(a,b,k=1)

	
	idxArray=cbind(nnn1$nn.idx,seq(length(nnn1$nn.idx)))
		
	# Need to supply a set of pairs of points.
	# will use the parent of each chosen point.
	# if parent undefined, then ignore that point
	
	if(is.null(dvs1) || is.null(dvs2)){
		if(OnlyClosestPoints==TRUE)
			stop("OnlyClosestPoints is not yet implemented for neurons")
		# Calculate the direction vectors
		dvs=findDirectionVectorsFromParents(n1,n2,idxArray,ReturnAllIndices=TRUE,Verbose=Verbose)

		# Calculate segment lengths
		l1.seglengths=normbyrow(dvs[,1:3])
		l2.seglengths=normbyrow(dvs[,4:6])
		# normalise the direction vectors
		dvs[,1:3]=dvs[,1:3]/l1.seglengths
		dvs[,4:6]=dvs[,4:6]/l2.seglengths
		# Calculate absolute dot products
		# nb absolute, because we don't really care about directionality here
		dps=abs(dotprod(dvs[,1:3],dvs[,4:6]))
	} else {
		# OnlyClosestPoints prunes the list of query-target pairs so that no 
		# points in the target are duplicated (points in query are already unique)
		if(OnlyClosestPoints){
			# sort by increasing distance between pairs
			# remove duplicates in target
			targetdupes=duplicated(nnn1$nn.idx[order(nnn1$nn.dist)])
			idxArray=idxArray[!targetdupes,,drop=FALSE]
			nnn1$nn.dists=nnn1$nn.dists[!targetdupes]
		}
		dps=abs(dotprod(dvs1[idxArray[,1],],dvs2[idxArray[,2],]))
	}

	NNDistFun(nnn1$nn.dists,dps,...)
}

findDirectionVectorsFromParents<-function(d1,d2,idxArray,ReturnAllIndices=FALSE,Verbose=FALSE){
	# rather than dropping root points, just use the vector from them rather than to them
	if(Verbose) cat(".")
	p1=.CleanupParentArray(d1[,"Parent"])
	p2=.CleanupParentArray(d2[,"Parent"])
	parentPointsArray=cbind(p1[idxArray[,1]],p2[idxArray[,2]])
	# find any points with bad parents and instead use their offspring
	if(any(parentPointsArray[,1]<1 | parentPointsArray[,2]<1)){
		stop ("Some points do not have a parent: therefore impossible to calculate direction vector")
	}
	dvs=cbind(d1[idxArray[,1],c("X","Y","Z"),drop=FALSE]-d1[parentPointsArray[,1],c("X","Y","Z"),drop=FALSE],
		d2[idxArray[,2],c("X","Y","Z"),drop=FALSE]-d2[parentPointsArray[,2],c("X","Y","Z"),drop=FALSE])

	if(ReturnAllIndices){
		attr(dvs,"idxArray")=idxArray
		attr(dvs,"parentPointsArray")=parentPointsArray
	}
	dvs
}

.CleanupParentArray<-function(pa){
	# takes a list of parents for points and replaces any <1
	# with the first offspring of that point
	#if(length(dim(pa))>1) apply(pa,2,.CleanupParentArray)
	pointsNeedingWork<-which(pa<1)
	if(length(pointsNeedingWork)<1) return( pa )
	for(p in pointsNeedingWork){
		wp=which(pa==p)
		if(length(wp)>1){
			warning(cat("more than 1 point in .CleanupParentArray, choosing first from:",wp))
			pa[p]=wp[1]
		} else if(length(wp)<1){
			warning("no points to choose in .CleanupParentArray using original value")
		} else pa[p]=wp
	}
	pa
}

.RemoveNNRepeats<-function(nnlist){
	# takes a result from nn2
	# and removes all duplicates but the shortest
	df=data.frame(nnlist)
	df$query.idx=seq(nrow(df))
	# sort by neighbour number 
	df2=df[order(df[,1],df[,2]),]
	df3=subset(df2,!duplicated(nn.idx))
	df3[order(df3$query.idx),]
}
