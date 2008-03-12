# LineMatching.R

# Code for matching lines in 3d

# First set up a distance measure

# After Kamgar-Parsi and Kamgar-Parsi, 2004
# http://ieeexplore.ieee.org/iel5/34/28505/01273930.pdf?arnumber=1273930
LinesetEuclideanDistance<-function(l1, l2){
	# l1 and l2 are sets of l+1 points each defining l segments
	
	l1=data.matrix(l1)
	l2=data.matrix(l2)
	
	# Calculate the direction vectors
	l1.diffs=diff(l1)
	l2.diffs=diff(l2)
	
	# Calculate squared segment lengths
	l1.seglengths=normbyrow(l1.diffs)
	l2.seglengths=normbyrow(l2.diffs)
		
	# normalise the direction vectors
	l1.diffs=l1.diffs/l1.seglengths
	l2.diffs=l2.diffs/l2.seglengths
	
	# if(!identical(l1.seglengths,l2.seglengths)){
	# 	stop("Can only compute distance for lines with segments of equal length")
	# }
	
	# Calculate the point displacement term
	pointDiffs=l1[-nrow(l1),]-l2[-nrow(l2),]
	m1=l1.seglengths*normbyrow(pointDiffs)

	# Calculate the line angle mismatch term	
	m2=l1.seglengths^3/6*(1-dotprod(l1.diffs,l2.diffs))
	
	mismatch=sum(m1+m2)
	if(mismatch < -1e-6) {
		stop("Negative line mismatch score!")
	} else if(mismatch<0){
		mismatch=0
	}
	return(mismatch)
}

NeuronLongestSegmentDistance<-function(ANeuron,BNeuron,n=10){
	# calculates the LinesetEuclideanDistance between two segments of a pair of neurons
	# aseg = # of segment f
	a=DivideLineIntoNEqualSubLines(GetLongestSegment(ANeuron),n)
	b=DivideLineIntoNEqualSubLines(GetLongestSegment(BNeuron),n)
	
	LinesetEuclideanDistance(a,b)
}

GetLongestSegment<-function(ANeuron){
	seg=which.max(ANeuron$SegLengths)
	return(ANeuron$d[ANeuron$SegList[[seg]],c("X","Y","Z")])
}