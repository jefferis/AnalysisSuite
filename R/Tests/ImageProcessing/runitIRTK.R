# Tests for Daniel Rueckert's image registration toolkit (IRTK)
# See: http://www.doc.ic.ac.uk/~dr/software

# runTestFile(file.path(TestDir,"ImageProcessing","runitIRTK.R"))

require(RUnit)

RunIdentityRegistration<-function(regtype){
	testData=matrix(rnorm(15),ncol=3)
	testLandmarks=list(testData,testData)
	
	dofout=tempfile()
	on.exit(unlink(dofout))
	irtk.preg(testLandmarks,dofout=dofout,xformtype=regtype)
	mat=irtk.dof2mat(dofout)
	
	identity=structure(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1),
		.Dim = c(4L, 4L))
	
	checkEqualsNumeric(mat,identity,
		msg=paste("Failure in",regtype,"identity registration"),tol=1e-6)
}

test.RunIdentityRegistrations<-function(){
	RunIdentityRegistration("rigid")
	RunIdentityRegistration("affine")
}

test.ExactAffineLandmarksRegistration<-function(){
	# take a point set
	points=structure(c(154.351776123047, 60.5672950744629, 234.804016113281, 
	244.884216308594, 97.8584060668945, 128.448928833008, 210.034393310547, 
	166.5546875, 116.681953430176, 214.466766357422, 97.4351043701172, 
	94.9338836669922, 115.311538696289, 171.031753540039, 178.53173828125, 
	198.142440795898, 200.835556030273, 110.578590393066, 50.0540008544922, 
	49.9504699707031, 24.9999694824219, 47, 14, 13, 12, 17.9999694824219, 
	17.9999694824219, 51, 70, 70), .Dim = c(10L, 3L))
	
	# and an affine transform
	xform=structure(c(0.942915828297041, 0.0179328718440451, -0.0203888393514411, 
	0, -0.0136353310033790, 0.910056271341136, 0.169845521757469, 
	0, 0.0437687209705831, -0.323941288649830, 1.73923692084124, 
	0, -0.216405156503470, 18.1009782432420, -30.7167612528349, 1
	), .Dim = c(4L, 4L))
	
	# and apply to get a new point set exactly related by an affine xform
	testLandmarks=list(points,TransformPoints(points,xform))
	
	# Now feed those points to irtk pareg
	# which will find an affine registration to map
	# points from target -> src
	dofout=tempfile()
	on.exit(unlink(dofout))
	irtk.preg(testLandmarks,dofout=dofout,xformtype='affine')
	xformcalc=irtk.dof2mat(dofout,Invert=TRUE)
	
	# check that we get back the same affine transform
	checkEqualsNumeric(xformcalc,xform,
		msg=paste("Failed to recover exact affine transform matrix"),tol=1e-3)
	
	# and that applying that transform using ptransformation
	# can map point set 1 => point set 2
	ps1x=irtk.transformation(testLandmarks[[1]], dofout)
	checkEqualsNumeric(ps1x,testLandmarks[[2]],
		msg="Failed to map pointset 1 => pointset 2 with calculated xform",tol=1e-3)	
}

test.dof2mat<-function(){
	correctmat=structure(c(-0.2813, -1.0275, -0.0183, 0, -0.9583, 0.3093, -0.0147, 
	0, 0.12, -0.1181, -1.269, 0, 223.2372, 188.0956, 144.6631, 1), .Dim = c(4L, 
	4L))
	mat=irtk.dof2mat(file.path(TestDir,"ImageProcessing","pareg.dof"))
	checkEqualsNumeric(correctmat,solve(mat),tol=1e-3)
}

# debug(test.RunIdentityRegistration)
# debug(test.dof2mat)