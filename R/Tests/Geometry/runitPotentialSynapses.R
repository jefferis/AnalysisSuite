#@+leo-ver=4-thin
#@+node:jefferis.20060305220950:@thin R/runitPotentialSynapses.R
#@@language r
# runit tests for the Potential synapse code
# Created: GSXEJ 2006-03-05

# runTestFile(file.path(TestDir,"Geometry","runitPotentialSynapses.R"))
require(RUnit)

test.PotentialSynapes.neuron<-function(){
  d1=structure(list(PointNo = 281:283, Label = c(2, 2, 2), X = c(114.158791, 
  115.422852, 116.922752), Y = c(94.342041, 95.47361, 97.593819
  ), Z = c(37.612789, 35.805931, 33.989944), W = c(1, 1, 1), Parent = c(280, 
  281, 282)), .Names = c("PointNo", "Label", "X", "Y", "Z", "W", 
  "Parent"), row.names = c("281", "282", "283"), class = "data.frame")
  
  d2=structure(list(PointNo = 297:299, Label = c(2, 2, 2), X = c(116.892075, 
  115.493431, 115.06131), Y = c(89.279106, 90.446068, 90.557793
  ), Z = c(36.778641, 36.041145, 36.880486), W = c(0.62, 0.62, 
  0.62), Parent = c(296, 297, 298)), .Names = c("PointNo", "Label", 
  "X", "Y", "Z", "W", "Parent"), row.names = c("297", "298", "299"
  ), class = "data.frame")

  n1=as.neuron(n=list(SegList=list(1:3),d=d1))
  n2=as.neuron(n=list(SegList=list(1:3),d=d2))
  # exactly one semgent pair in these fake neurons is within 5 microns
  checkEqualsNumeric(PotentialSynapses.neuron(n1,n2,s=5),2)
  checkEqualsNumeric(PotentialSynapses.neuron(n1,n2,s=4),1)
  checkEqualsNumeric(PotentialSynapses.neuron(n1,n2,s=3),0)
}

#@+others
#@+node:jefferis.20060305223034:Test DirectPotentialSynapses
test.DirectPotentialSynapses.parallelButNonOverlapping<-function(){
	
    s1=matrix( c( c(0,0,0), c(1,0,0) ), nrow=1)
    s2=matrix( c( c(2,0,0), c(3,0,0) ), nrow=1)
    
	checkEqualsNumeric(DirectPotentialSynapses(s1,s2),1)
  s3=matrix(
    c(0,0,0,1,0,0,
      10,10,10,11,0,0),ncol=6,byrow=T)
  s4=matrix(
    c(2,0,0,3,0,0,
      12,10,10,13,0,0),ncol=6,byrow=T)
  checkEqualsNumeric(DirectPotentialSynapses(s3,s4),2)
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
