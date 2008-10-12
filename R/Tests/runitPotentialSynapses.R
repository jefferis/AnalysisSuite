#@+leo-ver=4-thin
#@+node:jefferis.20060305220950:@thin R/runitPotentialSynapses.R
#@@language r
# runit tests for the Potential synapse code
# Created: GSXEJ 2006-03-05

# runTestFile(file.path(CodeDir,"runitPotentialSynapses.R"))
require(RUnit)


#@+others
#@+node:jefferis.20060305223034:Test DirectPotentialSynapses
test.DirectPotentialSynapses.parallelButNonOverlapping<-function(){
	
    s1=matrix( c( c(0,0,0), c(1,0,0) ), nrow=1)
    s2=matrix( c( c(2,0,0), c(3,0,0) ), nrow=1)
    
	checkEqualsNumeric(DirectPotentialSynapses(s1,s2),1)
}
#@-node:jefferis.20060305223034:Test DirectPotentialSynapses
#@+node:jefferis.20060305222740:Test dist3D_Segment_to_Segment
# These are the test for dist3D_Segment_to_Segment, the primitive function for
# calculating distances between two segments 
# presently passes all tests

test.dist3D_Segment_to_Segment.parallelButNonOverlapping<-function(){
	
    s1=c( c(0,0,0), c(1,0,0) )
    s2=c( c(2,0,0), c(3,0,0) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),1)
}


test.dist3D_Segment_to_Segment.parallelOverlapping<-function(){
	
    s1=c( c(0,0,0), c(1,0,0) )
    s2=c( c(0.5,0,0), c(1.5,0,0) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),0)
}


test.dist3D_Segment_to_Segment.orthogonalOverlapping<-function(){
	
    s1=c( c(0,0,0), c(1,0,0) )
    s2=c( c(0.5,0,0), c(0.5,1,0) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),0)
}

test.dist3D_Segment_to_Segment.diagonallyDisplaced<-function(){
	  # /
     # /
    # /   
        # \
         # \
          # \ 

    
    s1=c( c(0,0,0), c(1,1,0) )
    s2=c( c(1,0,0), c(2,-1,0) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),sqrt(2)/2)
}


test.dist3D_Segment_to_Segment.crossOver<-function(){
	# \ /
    #  X
    # / \ 
    
    s1=c( c(0,0,0), c(1,1,0) )
    s2=c( c(0,1,0), c(1,0,0) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),0)
}

test.dist3D_Segment_to_Segment.crossOverZdisplaced<-function(){
	# \ /
    #  X
    # / \ 
    
    s1=c( c(0,0,0), c(1,1,0) )
    s2=c( c(0,1,.3), c(1,0,.3) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),0.3)
}

test.dist3D_Segment_to_Segment.crossOverZFardisplaced<-function(){
	# \ /
    #  X
    # / \ 
    
    s1=c( c(0,0,0), c(1,1,0) )
    s2=c( c(0,1,10), c(1,0,10) )
    
	checkEqualsNumeric(dist3D_Segment_to_Segment(s1,s2),10)
}
#@nonl
#@-node:jefferis.20060305222740:Test dist3D_Segment_to_Segment
#@-others
#@nonl
#@-node:jefferis.20060305220950:@thin R/runitPotentialSynapses.R
#@-leo
