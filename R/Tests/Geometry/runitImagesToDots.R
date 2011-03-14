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
}