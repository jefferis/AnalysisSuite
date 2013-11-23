# runitAffine.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitAffine.R"))
require(RUnit)

checkM=function(m){
	m=matrix(m,ncol=3,byrow=TRUE)
	checkEquals(m,DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,]))
}

printM=function(m){
	if(!is.matrix(m)) m=matrix(m,ncol=3,byrow=TRUE)
	print(m)
	print(DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,]))
}


test.ComposeAffineNoShear<-function(){
	params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0,0,0,10,50,50),
		ncol=3,byrow=TRUE)
	affmat=matrix(c(1.09699,0.0494996,0.0494536,94.0825,
		-0.0574907,0.897406,-0.0549995,58.4546,
		-0.0575696,0.0470378,0.997261,8.36076,
		0,0,0,1),ncol=4,byrow=TRUE)
	checkEqualsNumeric(ComposeAffineFromIGSParams(params),affmat,tol=1e-6)
}

test.ComposeAffineWithShear<-function(){
	params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0.03,0.1,0.05,10,50,50),
		ncol=3,byrow=TRUE)
	affmat=matrix(c(1.09878,0.0599299,0.104303,90.8005,
		-0.0307421,0.891618,-0.0578741,58.6202,
		-0.0531753,0.146476,0.994382,3.48883,
		0,0,0,1),ncol=4,byrow=TRUE)
	checkEqualsNumeric(ComposeAffineFromIGSParams(params),affmat,tol=1e-6)
}


test.ReCompositionAffineNoShear<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0,0,0,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear1<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0,0,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear2<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0,0.05,0,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear3<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0,0,0.05,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear12<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0.05,0,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear12DiffShears<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0.13,0,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear13<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0,0.05,10,50,50)
	m=matrix(m,ncol=3,byrow=TRUE)
	mm=DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,])
	checkEquals(m,mm,
	msg=paste(m,'not equal',mm) )
}
test.ReCompositionAffineShear13DiffShears<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0,0.13,10,50,50)
checkM(m)}

test.ReCompositionAffineShear23<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0,0.05,0.05,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear123<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0.05,0.05,10,50,50)
	m=matrix(m,ncol=3,byrow=TRUE)
	mm=DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,])
	checkEquals(m,mm,paste(m,'not equal',mm) )
}

test.ReCompositionAffineShear123NoCentre<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0.05,0.05,0,0,0)
	checkM(m)
}
