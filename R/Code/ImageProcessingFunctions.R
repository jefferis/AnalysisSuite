# General Image Processing Functions combining CMTK / UNU or my R Native Code
# distinct from ImageAnalyisFunctions which are more statistical in nature
# this is really about munging files

ResampleAndFlipMasks<-function(masks,outdir,HorizontalBridgingReg,flipAxis=c("X","Y","Z"),targetspec){
	flipAxis=match.arg(flipAxis)
	if(!file.exists(outdir)) dir.create(outdir)
	for (infile in masks){
		resampledfile=file.path(outdir,
			sub(".nrrd$","-resampled.nrrd",basename(infile)))
		flippedresampledfile=file.path(outdir,
			sub(".nrrd$","-resampled-flip.nrrd",basename(infile)))
		orfile=file.path(outdir,
			sub(".nrrd$","-OR.nrrd",basename(infile)))

		if(!exists("identityReg") || !file.exists(identityReg))
			identityReg=WriteIdentityRegistration()
		# resample - use reformatx to do this, to ensure that we get the same result
		ReformatImage(infile,target=targetspec, registrations=identityReg,
			filesToIgnoreModTimes=identityReg, OverWrite='update',
			output=resampledfile,reformatoptions="-v --pad-out 0 --nn",dryrun=FALSE)

		# make flipping registration
		if(!exists("horizontalFlipReg") || !file.exists(horizontalFlipReg))
			horizontalFlipReg=WriteHorizontalFlipRegistration(resampledfile,axis=flipAxis)
		# and flip all masks, applying the horiz bridging registration
		ReformatImage(resampledfile,target=resampledfile,
			registrations=c(HorizontalBridgingReg,horizontalFlipReg),
			filesToIgnoreModTimes=horizontalFlipReg, OverWrite='update', Verbose=T,
			output=flippedresampledfile,reformatoptions="-v --pad-out 0 --nn",dryrun=FALSE)
		# and OR ing those results
		Nrrd2op(c(resampledfile,flippedresampledfile),orfile,'max')
	}
	unlink(horizontalFlipReg)
}
