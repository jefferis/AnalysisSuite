# runitImagesToDots.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitImagesToDots.R"))
require(RUnit)
require(R.matlab)

test.DotProperties<-function(){
	dots4=t(readMat(file.path(TestDir,"Geometry","SAKW13-1_dots4.mat"))[[1]])
	props4=readMat(file.path(TestDir,"Geometry","SAKW13-1_dots4_properties.mat"))
	# reorder matrices into R forms
	props4$alpha=as.numeric(props4$alpha)
	props4$vect=t(props4$vect)
	
	props4.new<-DotProperties(dots4)
	
	checkEqualsNumeric(props4.new$alpha,props4$alpha,tol=1e-6)
	# should probably just check that dot product is +/- 1
	checkEqualsNumeric(abs(props4.new$vect),abs(props4$vect),tol=1e-6)
	#checkEqualsNumeric(props4.new$vect,props4$vect,tol=1e-6)
}

test.ind2coord<-function(){
	# from matlab:
	# ind2coord([1024 768 100],[1 1034,10001,100001],[0.3,0.4,1])
	# 0.3000    3.0000  235.5000  201.9000
    # 0.4000    0.8000    4.0000   39.2000
    # 1.0000    1.0000    1.0000    1.0000
    
	i2c.matlab=t(structure(c(0.3, 0.4, 1, 3, 0.8, 1, 235.5, 4, 1, 201.9, 39.2, 
	1), .Dim = 3:4))
	
	i2c.r=ind2coord(c(1,1034,10001,100001),c(1024,768,100),c(0.3,0.4,1))
	
	checkEqualsNumeric(i2c.matlab,i2c.r,tol=1e-6)
	
	lh=Read3DDensityFromAmiraLattice(file.path(ObjDir,"LHMask.am"))
	x=ind2coord(lh)
	y=ind2coord(lh>0,voxdims=lh)
	checkEqualsNumeric(x,y,tol=1e-6)
	
	xinds=c(1223L, 4282L, 4907L, 5978L, 10327L, 12151L, 13042L, 19879L, 
	23842L, 27664L)
	x.correct=structure(c(133.5, 141.9, 127.9, 126.5, 127.9, 129.3, 129.3, 
	139.1, 137.7, 139.1, 87.3, 98.5, 92.9, 99.9, 111.1, 95.7, 90.1, 
	97.1, 73.3, 90.1, 13.3, 18.9, 20.3, 21.7, 27.3, 30.1, 31.5, 41.3, 
	48.3, 56.7), .Dim = c(10L, 3L), .Dimnames = list(NULL, c("X", 
	"Y", "Z")))
	
	checkEqualsNumeric(x[xinds,],x.correct)
}

test.sub2ind<-function(){
	# from matlab
	# sub2ind([512 768 112], 300, 200, 77)
	
	s2i.matlab = 29986604
	s2i = sub2ind(c(512,768,112),c(300,200,77))
	checkEqualsNumeric(s2i.matlab,s2i)
	
	s2i.correct = c(19762988,19974032)
	s2i = sub2ind(c(512,512,77),matrix(c(300,200,76,400,100,77),byrow=T,ncol=3))
	checkEqualsNumeric(s2i.correct,s2i)
}

