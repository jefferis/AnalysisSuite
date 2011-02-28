# runitImageJFileFunctions.R

# Some test neurons to verify that all is well for the 
# ImageJ LUT reading rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(TestDir,"IO","runitImageJFileFunctions.R"))
require(RUnit)

ImageJDirectory=NULL

if(file.exists("/Applications/Fiji.app")){
	ImageJDirectory="/Applications/Fiji.app"	
} else if(file.exists("/Applications/Fiji/Fiji.app")){
	ImageJDirectory="/Applications/Fiji/Fiji.app"
} else if(file.exists("/Applications/ImageJ")){
	ImageJDirectory="/Applications/ImageJ"
}

test.ReadWriteImageJLUTs<-function(){
	tmpfile=tempfile()
	if(is.null(ImageJDirectory)){
		warning("Unable to find ImageJ test LUTs")
		return(1)
	}
	luts=list.files(file.path(ImageJDirectory,'luts'),patt='\\.lut$',full=T)
	for(lut in luts){
		cat("Testing with lut",lut,"\n")
		l=ReadImageJLUT(lut)
		WriteImageJLUT(tmpfile,l,ftype='binary')
		l2=ReadImageJLUT(tmpfile)
		checkEquals(l,l2,tol=1e-6)
	}
}