# runitAffine.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitAffine.R"))
require(RUnit)

checkM=function(m){
	m=matrix(m,ncol=3,byrow=TRUE)
	checkEquals(m,DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,]))
}

checkM2=function(m){
  m=matrix(m,ncol=3,byrow=TRUE)
  checkEquals(m,DecomposeAffineToIGSParams(HomogenousAffineFromCMTKParams(m),cent=m[5,]))
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
test.DecomposeAffineWithShear<-function(){
  params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0.03,0.1,0.05,0,0,0),
                ncol=3,byrow=TRUE)
  affmat=matrix(c(1.09878,0.0599299,0.104303,90.8005,
                  -0.0307421,0.891618,-0.0578741,58.6202,
                  -0.0531753,0.146476,0.994382,3.48883,
                  0,0,0,1),ncol=4,byrow=TRUE)
  checkEqualsNumeric(DecomposeAffineToIGSParams(affmat),params,tol=1e-6)
}
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

test.HomogenousAffineFromCMTKParams<-function(){
  params=matrix(c(100,50,10, 3,3,3, 1.1,0.9,1, 0,0,0, 0,0,0), ncol=3, byrow=TRUE)
  m=structure(c(1.09699, -0.0574907, -0.0575696, 0, 0.0494996, 0.897406, 
              0.0470378, 0, 0.0494536, -0.0549995, 0.997261, 0, 100, 50, 10, 
              1), .Dim = c(4L, 4L))
  m2=HomogenousAffineFromCMTKParams(params)
  checkEqualsNumeric(m,m2)
}

#' Compare result of ComposeAffineFromIGSParams vs CMTK dof2mat
test.ComposeAffineFromIGSParamsvsCMTK<-function(){ 
  params=matrix(c(100,50,50,3,3,3,1.1,0.9,1,0,0,0,0,0,0),
                ncol=3,byrow=TRUE)
  m=ComposeAffineFromIGSParams(params)
  # nb this calls dof2mat
  m2=HomogenousAffineFromCMTKParams(params)
  checkEqualsNumeric(m,m2,'cmtk and ComposeAffineFromIGSParams disagree',
                     tolerance=1e-5)
}

#' Check DecomposeAffineToIGSParams vs CMTK mat2dof
test.DecomposeAffineFromIGSParamsvsCMTK<-function(){
  m=matrix(c(1.1,0,0,50,
             0,1.2,0,60,
             0,0,1.1,20,
             0,0,0,1),ncol=4,byrow=TRUE)
  
  params=DecomposeAffineToIGSParams(m,centre=c(0,0,0))
  params2=cmtk.mat2dof(m)
  checkEqualsNumeric(params,params2,paste(params,'not equal',params2),
                     tolerance=1e-5)
}

#' round trip test of mat2dof/dof2mat
test.cmtk.mat2dof.dof2mat.simple<-function(){
  m=matrix(c(1.1,0,0,50,
             0,1.2,0,60,
             0,0,1.1,20,
             0,0,0,1),ncol=4,byrow=TRUE)
  tf<-tempfile(fileext='.list')
  dir.create(tf)
  on.exit(unlink(tf,recursive=TRUE))
  cmtk.mat2dof(m,f=tf)
  m2=cmtk.dof2mat(tf)
  checkEquals(m,m2,"Failed CMTK mat2dof/dof2mat roundtrip test for simple matrix (no shears)")
}

test.cmtk.mat2dof.dof2mat.wshears<-function(){
  m=structure(c(0.911236, 0.00295678, 0.00483363, 0, -0.045134, 1.105, 
                -0.104863, 0, -0.0475781, 0.0580714, 0.997261, 0, -88.3912, -56.1266, 
                -5.21284, 1), .Dim = c(4L, 4L))
  tf<-tempfile(fileext='.list')
  dir.create(tf)
  on.exit(unlink(tf,recursive=TRUE))
  cmtk.mat2dof(m,f=tf)
  m2=cmtk.dof2mat(tf)
  checkEqualsNumeric(m,m2,"Failed CMTK mat2dof/dof2mat roundtrip test for matrix with shears")
}

test.cmtk.mat2dof<-function(x){
  m=structure(c(0.993768, -0.0869434, -0.0697565, 0, 0.199117, 1.08527, 
                0.0504537, 0, 0.303757, 0.211115, 1.19715, 0, 100, 50, 50, 1),
              .Dim = c(4L, 4L))
  params=cmtk.mat2dof(m)
  params_base=matrix(c(100,50,50, 3,4,5, 1,1.1,1.2, 0.1,0.2,0.3, 0,0,0), ncol=3,
                     byrow=T)
  checkEqualsNumeric(params,params_base,tolerance=1e-4)
}

test.cmtk.dof2mat<-function(x){
  
  reg=file.path(TestDir,'Data','dofv2.4wshears.list')
  params=matrix(c(100,50,50, 3,4,5, 1,1.1,1.2, 0.1,0.2,0.3, 0,0,0), ncol=3,
                byrow=T)
  m_base=structure(c(0.993768, -0.0869434, -0.0697565, 0, 0.199117, 1.08527, 
                0.0504537, 0, 0.303757, 0.211115, 1.19715, 0, 100, 50, 50, 1),
              .Dim = c(4L, 4L))
  checkEqualsNumeric(cmtk.dof2mat(reg),m_base,tolerance=1e-4)
  checkEqualsNumeric(cmtk.dof2mat(params),m_base,tolerance=1e-4)
}
