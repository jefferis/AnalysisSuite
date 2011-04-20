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

	# from matlab with axperm
	# ind2coord([1024 768 100],[1 1034,10001,100001],[0.3,0.4,1],[1 2 3])
	#          0    2.7000  235.2000  201.6000
	#          0    0.4000    3.6000   38.8000
	#          0         0         0         0
    
	i2c.matlab=t(structure(c(0, 0, 0, 2.7, 0.4, 0, 235.2, 3.6, 0, 201.6, 38.8, 0),
		.Dim = 3:4))
	
	i2c.r=ind2coord(c(1,1034,10001,100001),c(1024,768,100),c(0.3,0.4,1))
	
	checkEqualsNumeric(i2c.matlab,i2c.r,tol=1e-6)
	
	lh=Read3DDensityFromAmiraLattice(file.path(ObjDir,"LHMask.am"))
	x=ind2coord(lh)
	y=ind2coord(lh>0,voxdims=lh)
	checkEqualsNumeric(x,y,tol=1e-6)
	
	xinds=c(1223L, 4282L, 4907L, 5978L, 10327L, 12151L, 13042L, 19879L, 
	23842L, 27664L)
	x.correct=structure(c(132.1, 140.5, 126.5, 125.1, 126.5, 127.9, 127.9, 
	137.7, 136.3, 137.7, 85.9, 97.1, 91.5, 98.5, 109.7, 94.3, 88.7, 
	95.7, 71.9, 88.7, 11.9, 17.5, 18.9, 20.3, 25.9, 28.7, 30.1, 39.9, 
	46.9, 55.3), .Dim = c(10L, 3L), .Dimnames = list(NULL, c("X", 
	"Y", "Z")))
	
	checkEqualsNumeric(x[xinds,],x.correct)
}

test.coord2ind<-function(){
	# from matlab
	# c2i.matlab=1184259
	# c2i=coord2ind(c(1,3,4),imdims=c(512,768,112),voxdims=c(0.31,0.31,1.06))
	# checkEqualsNumeric(c2i.matlab,c2i)

	c2i.origin=1
	c2i=coord2ind(c(0,0,0),imdims=c(512,768,112),voxdims=c(0.31,0.31,1.06))
	checkEqualsNumeric(c2i.origin,c2i)

	c2i.origin=1
	c2i=coord2ind(c(40,50,60),imdims=c(512,768,112),voxdims=c(2,3,4),origin=c(40,50,60))
	checkEqualsNumeric(c2i.origin,c2i)

	c2i.origin=(80-60)/4*(768*512)+1
	c2i=coord2ind(c(40,50,80),imdims=c(512,768,112),voxdims=c(2,3,4),origin=c(40,50,60))
	checkEqualsNumeric(c2i.origin,c2i)

	c2i.origin=(80-60)/4*(768*512)+(60-40)/2+1
	c2i=coord2ind(c(60,50,80),imdims=c(512,768,112),voxdims=c(2,3,4),origin=c(40,50,60))
	checkEqualsNumeric(c2i.origin,c2i)
}

test.coord2ind.roundtrip<-function(){
	lh=Read3DDensityFromAmiraLattice(file.path(ObjDir,"LHMask.am"))
	coords=ind2coord(lh)
	inds=coord2ind(coords,lh)
	lhorigin=getBoundingBox(lh)[c(1,3,5)]
	coordsagain=ind2coord(inds,dim=dim(lh),voxdim=voxdim.gjdens(lh),origin=lhorigin)
	checkEqualsNumeric(coords,coordsagain,tol=1e-6)
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

