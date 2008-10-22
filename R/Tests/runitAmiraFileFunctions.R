# runitAmiraFileFunctions.R

# Some test neurons to verify that all is well for the 
# AmiraFileFunctions rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(CodeDir,"runitAmiraFileFunctions.R"))
require(RUnit)

test.DoubleFromSingleEdgeList<-function(){
	
	el=rbind(c(1,2),c(2,3),c(2,4))
	
	# these should obviously be inverses of each other!
	checkEquals(el,SingleFromDoubleEdgeList(DoubleFromSingleEdgeList(el)))
}
	

test.ParseEdgeList<-function(){
	SimpleY <-
	structure(list(CurPoint = c(1, 2, 2, 2, 3, 4, 4, 5), Neighbour = c(2, 
	1, 3, 4, 2, 2, 5, 4)), .Names = c("CurPoint", "Neighbour"), row.names = c("1", 
	"2", "3", "4", "5", "6", "7", "8"), class = "data.frame")
	SimpleYResult<-
	list(c(1, 2), c(2, 3), c(2, 4, 5))

	checkEquals(ParseEdgeList(SimpleY),SimpleYResult)
	
	# QuadBranch
	# a branchpoint with 4 neighbours total - ie trifurcation
	QuadBranch <-
	structure(list(CurPoint = c(1, 2, 2, 2, 2, 3, 4, 4, 5, 6), Neighbour = c(2, 
	1, 3, 4, 6, 2, 2, 5, 4, 2)), .Names = c("CurPoint", "Neighbour"
	), row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
	"10"), class = "data.frame")

	QuadBranchResult <-
	list(c(1, 2), c(2, 3), c(2, 4, 5), c(2, 6))

	checkEquals(ParseEdgeList(QuadBranch),QuadBranchResult)
	
	NestedBranches <-
	structure(list(CurPoint = c(1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 5, 
	6, 7, 8), Neighbour = c(2, 1, 3, 4, 6, 2, 2, 5, 4, 7, 8, 2, 5, 
	5)), .Names = c("CurPoint", "Neighbour"), row.names = c("1", 
	"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
	"14"), class = "data.frame")

	# NestedBranchesResult=ParseEdgeList(NestedBranches)
	NestedBranchesResult<-
	list(c(1, 2), c(2, 3), c(2, 4, 5), c(5, 7), c(5, 8), c(2, 6))

	#test("NestedBranches")
	checkEquals(ParseEdgeList(NestedBranches),NestedBranchesResult)

	
	# a branch right after a previous branch
	el=rbind( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6) )
	RapidlyNestedBranchesResult<-
	list( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6) )

	checkEquals(ParseEdgeList(DoubleFromSingleEdgeList(el)),
		RapidlyNestedBranchesResult,msg="One Branch Point right after another")

	# a missing point number
	el=DoubleFromSingleEdgeList(rbind( c(1,2),c(2,3),c(2,5),c(5,6),c(5,7) ))
	MissingPointResult=list( c(1,2),c(2,3),c(2,5),c(5,6),c(5,7) )
	#checkException(ParseEdgeList(el), MissingPointResult)
	#DEACTIVATED("MissingPointResult")

	# unattached points
	el=DoubleFromSingleEdgeList(rbind( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6),c(7,8)))
	
	UnattachedPointResult=list( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6) )
	#checkException(ParseEdgeList(el), "unattached points")
	checkEquals(ParseEdgeList(el,Silent=T),UnattachedPointResult, "unattached points with Force=T")

	# starting from a different root
	el=rbind( c(1,2),c(2,3),c(2,4),c(4,5),c(4,6) )
	DifferentRootResult<-
	list( c(6,4),c(4,2),c(2,1),c(2,3),c(4,5) )

	checkEquals(ParseEdgeList(DoubleFromSingleEdgeList(el),RootPoint=6),
		DifferentRootResult,msg="Root other than point 1")

	
	# Nonsequential numbering
	el=rbind( c(1,2),c(2,6),c(2,4),c(4,5),c(4,3) )
	NonsequentialNumberingResult<-
	list( c(1,2),c(2,4),c(4,3),c(4,5),c(2,6) )

	checkEquals(ParseEdgeList(DoubleFromSingleEdgeList(el)),
		NonsequentialNumberingResult,msg="Nonsequential numbering")


}

test.ReadWrite3DDensityAmiraBinary<-function(){
	testData=array(rnorm(10^3),dim=rep(10,3))
	tmpfile=tempfile()
	
	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="binary",dtype='double')
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	unlink(tmpfile)
	checkEquals(testData,testData.new)
}

test.ReadWrite3DDensityAmiraText<-function(){
	testData=array(rnorm(10^3),dim=rep(10,3))
	tmpfile=tempfile()
	
	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="text")
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	unlink(tmpfile)
	checkEquals(testData,testData.new,tol=1e-6)
}

test.ReadWriteAmiraLandmarksSingle<-function(){
	testData=matrix(rnorm(15),ncol=3)
	tmpfile=tempfile()
	WriteAmiraLandmarks(tmpfile,testData)
	testData.new=ReadAmiraLandmarks(tmpfile)
	unlink(tmpfile)
	names(testData.new)<-NULL
	checkEquals(testData,testData.new,tol=1e-6)	
}

test.ReadWriteAmiraLandmarksPaired<-function(){
	testData=replicate(2,matrix(rnorm(15),ncol=3),simplify=FALSE)
	tmpfile=tempfile()

	WriteAmiraLandmarks(tmpfile,testData)
	testData.new=ReadAmiraLandmarks(tmpfile)
	unlink(tmpfile)
	names(testData.new)<-NULL
	checkEquals(testData,testData.new,tol=1e-6)	
}