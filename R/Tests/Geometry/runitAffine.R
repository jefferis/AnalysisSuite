# runitAffine.R
# functions to test affine functions

# runTestFile(file.path(TestDir,"Geometry","runitAffine.R"))
require(RUnit)

checkM=function(m){
	m=matrix(m,ncol=3,byrow=TRUE)
	checkEqualsNumeric(m,DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,]),tolerance=1e-5)
}

cmtk_version=cmtk.dof2mat(version=TRUE)
cmtk_numeric_version=numeric_version(sub("([0-9.]+).*",'\\1',cmtk_version))

test.cmtk.version<-function(){
  checkTrue(cmtk_numeric_version>=numeric_version("2.4.0"),
            msg='AnalysisSuite now depends on CMTK>=2.4.0')
}

checkRoundTripFromMat=function(mat,
  test.cmtk=cmtk_numeric_version>=numeric_version("2.4.0")){
  params=DecomposeAffineToIGSParams(mat)
  m2=ComposeAffineFromIGSParams(params)
  checkEqualsNumeric(mat,m2,tolerance=1e-5)
  if(test.cmtk){
    # repeat with cmtk tools
    params2=cmtk.mat2dof(mat)
    m3=cmtk.dof2mat(params2)
    # can't absolutely rely on getting same params back in certain degenerate 
    # cases e.g. axis flip.
    # checkEqualsNumeric(params,params2,tolerance=1e-5)
    checkEqualsNumeric(mat,m3,tolerance=1e-5)
  }
}

#' Round trip using cmtk command line tools
checkM.cmtk=function(params){
  params=matrix(params,ncol=3,byrow=TRUE)
  affmat=cmtk.dof2mat(params)
  params2=cmtk.mat2dof(affmat,cent=params[5,])
  checkEqualsNumeric(params,params2,tolerance=1e-4)
}

printM=function(m){
	if(!is.matrix(m)) m=matrix(m,ncol=3,byrow=TRUE)
	print(m)
	print(DecomposeAffineToIGSParams(ComposeAffineFromIGSParams(m),cent=m[5,]))
}

test.ComposeAffineNoShear<-function(){
	params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0,0,0,10,50,50),
		ncol=3,byrow=TRUE)
	affmat=matrix(c(1.09698704245255, -0.0574906547972094, -0.0575695518672382, 
	                0, 0.0494995771563172, 0.897405837085788, 0.0470378084704441, 
	                0, 0.0494535530049303, -0.0549995301736858, 0.997260947684137, 
	                0, 94.0824730674121, 58.454591202367, 8.36075771094335, 1),
                ncol=4)
	checkEqualsNumeric(ComposeAffineFromIGSParams(params),affmat,tol=1e-6)
}

test.ComposeAffineWithShear<-function(){
	params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0.03,0.1,0.05,10,50,50),
		ncol=3,byrow=TRUE)
	affmat=matrix(c(1.09698704245255, -0.0574906547972094, -0.0575695518672382, 
	                0, 0.0794174055868413, 0.895837910136773, 0.0454677297831557, 
	                0, 0.151929624282028, -0.0103700734989691, 0.994640563641534, 
	                0, 87.462778082031, 56.3015147160819, 8.57028084743792, 1),
                ncol=4)
	checkEqualsNumeric(ComposeAffineFromIGSParams(params),affmat,tol=1e-6)
}

test.DecomposeAffineWithShear<-function(){
  params=matrix(c(100,50,10,3,3,3,1.1,0.9,1,0.03,0.1,0.05,0,0,0),
                ncol=3,byrow=TRUE)
  affmat=matrix(c(1.09698704245255, -0.0574906547972094, -0.0575695518672382, 
                  0, 0.0794174055868413, 0.895837910136773, 0.0454677297831557, 
                  0, 0.151929624282028, -0.0103700734989691, 0.994640563641534, 
                  0, 100, 50, 10, 1),ncol=4)
  checkEqualsNumeric(DecomposeAffineToIGSParams(affmat),params,tol=1e-6)
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
	checkM(m)
}
test.ReCompositionAffineShear13DiffShears<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0,0.13,10,50,50)
  checkM(m)
}

test.ReCompositionAffineShear23<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0,0.05,0.05,10,50,50)
	checkM(m)
}

test.ReCompositionAffineShear123<-function(x){
	m=c(100,50,10,3,3,3,1.1,0.9,1,0.05,0.05,0.05,10,50,50)
  checkM(m)
}

test.ReCompositionAffineShear123NoCentre<-function(x){
	m=c(100,50,10,3,4,5,1.1,0.9,1,0.05,0.02,0.03,0,0,0)
	checkM(m)
}

test.ReCompositionAffineShear123NoCentreMirrorX<-function(x){
  params=c(100,50,10,3,4,5,-1.1,0.9,1,0.05,0.02,0.03,0,0,0)
  mat=ComposeAffineFromIGSParams(params)
  checkRoundTripFromMat(mat)
}

test.ReCompositionAffineShear123NoCentreMirrorXY<-function(x){
  params=c(100,50,10,3,4,5,-1.1,-0.9,1,0.05,0.02,0.03,0,0,0)
  mat=ComposeAffineFromIGSParams(params)
  checkRoundTripFromMat(mat)
}

test.ReCompositionAffineShear123NoCentreMirrorXYZ<-function(x){
  params=c(100,50,10,3,4,5,-1.1,-0.9,-1,0.05,0.02,0.03,0,0,0)
  mat=ComposeAffineFromIGSParams(params)
  checkRoundTripFromMat(mat)
}

test.ReCompositionAffineShear123NoCentreMirrorY<-function(x){
  params=c(100,50,10,3,4,5,1.1,-0.9,1,0.05,0.02,0.03,0,0,0)
  mat=ComposeAffineFromIGSParams(params)
  checkRoundTripFromMat(mat)
}

test.ReCompositionAffineShear123NoCentreMirrorZ<-function(x){
  params=c(100,50,10,3,4,5,1.1,0.9,-1,0.05,0.02,0.03,0,0,0)
  mat=ComposeAffineFromIGSParams(params)
  checkRoundTripFromMat(mat)
}

test.ReCompositionSimpleMirrorX<-function(x){
  m=matrix(0,4,4)
  diag(m)=c(-1,1,1,1)
  checkRoundTripFromMat(m)
}

test.ReCompositionSimpleMirrorY<-function(x){
  m=matrix(0,4,4)
  diag(m)=c(1,-1,1,1)
  checkRoundTripFromMat(m)
}

test.ReCompositionSimpleMirrorZ<-function(x){
  m=matrix(0,4,4)
  diag(m)=c(1,1,-1,1)
  checkRoundTripFromMat(m)
}
if(cmtk_numeric_version>=numeric_version("2.4.0")){
#' Compare result of ComposeAffineFromIGSParams vs CMTK dof2mat
test.ComposeAffineFromIGSParamsvsCMTK<-function(){ 
  params=matrix(c(100,50,20,3,4,5,1.1,0.9,1,0.05,0.1,0.02,10,20,30),
                ncol=3,byrow=TRUE)
  m=ComposeAffineFromIGSParams(params)
  # nb this calls dof2mat
  m2=cmtk.dof2mat(params)
  checkEqualsNumeric(m,m2,'cmtk and ComposeAffineFromIGSParams disagree',
                     tolerance=1e-5)
}

#' Check DecomposeAffineToIGSParams vs CMTK mat2dof
test.DecomposeAffineFromIGSParamsvsCMTK<-function(){
  m=matrix(c(1.1,0.1,0.2,50,
             0.1,1.2,.25,60,
             0.3,0.4,1.1,20,
             0,0,0,1),ncol=4,byrow=TRUE)
  
  params=DecomposeAffineToIGSParams(m,centre=c(0,0,0))
  params2=cmtk.mat2dof(m)
  checkEqualsNumeric(params,params2,tolerance=1e-5)
}


}