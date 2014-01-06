# ReadIGSFiles.R
#
# Created by Greg Jefferis 2006-01-08
#
# Functions to read the generic IGS TypedStream data
# and specifically the registration files (warp or affine)

#RELEASE
#BEGINCOPYRIGHT
###############
# R Source Code to accompany the manuscript
#
# "Comprehensive Maps of Drosophila Higher Olfactory Centers: 
# Spatially Segregated Fruit and Pheromone Representation"
# Cell (2007), doi:10.1016/j.cell.2007.01.040
# by Gregory S.X.E. Jefferis*, Christopher J. Potter*
# Alexander M. Chan, Elizabeth C. Marin
# Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo
#
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# See flybrain.stanford.edu for further details
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################
#ENDMAINCOPYRIGHT

#' Read CMTK(IGS) format registration
#'
#' @details Note that ReturnRegistrationOnly defaults to FALSE in the new
#' read.cmtkreg function but we have kept the old default of TRUE for calls
#' coming in via ReadIGSRegistration
ReadIGSRegistration <- function (filename, ReturnRegistrationOnly=NULL){
	.Deprecated('read.cmtkreg','nat')
	if(is.null(ReturnRegistrationOnly)) ReturnRegistrationOnly=TRUE
	read.cmtkreg(filename,ReturnRegistrationOnly)
}

#' Read CMTK TypedStream file to a list in memory
#' 
#' @details This is the default format used by CMTK for registration, studylist 
#'   and images files.
#' @param con Path to (optionally gzipped) file or (open) connection.
#' @param CheckLabel Check, fix and warn for invalid or duplicate labels (by
#'   default)
ReadIGSTypedStream<-function(con, CheckLabel=TRUE){
	.Deprecated('read.cmtk','nat')
	read.cmtk(con,CheckLabel=CheckLabel)
}

#' Write a suitable list to a CMTK TypedStream file on disk
#' 
#' @details NB a version specified on the command line overrides one encoded as 
#'   an attribute in the input list.
#' @param version TYPEDSTREAM version number, defaults to \code{"1.1"} if not 
#'   specified in the version attribute of \code{l}.
#' @export
WriteIGSTypedStream<-function(l, filename, gzip=FALSE, version=NA_character_){
	.Deprecated("write.cmtk",'nat')
	nat::write.cmtk(l,filename,gzip=gzip,version=version)
}

IGSParamsToIGSRegistration<-function(x,reference="dummy",model="dummy"){
	# Note that this produces a list that "should" be identical to what
	# would be produced by reading in an IGS registration file
	# This could now be written out wiith WriteIGSTypedStream
	affine_xform=unlist(apply(x,1,list),rec=F)
	names(affine_xform)=c("xlate","rotate","scale","shear","center")
	warning("Use of model_study entry in CMTK registrations is deprecated.",
	        " See CMTKParamsToCMTKRegistration")
	l=list(registration=list(reference_study=reference,
						model_study=model,
						affine_xform=affine_xform))
	attr(l$registration$reference_study,"quoted")=TRUE
	attr(l$registration$model_study,"quoted")=TRUE
	l
}

#' Produce in-memory list representation from CMTK affine parameters
#' 
#' @details Note that this uses the modern CMTK notation of floating_study 
#'   rather than model_study as used by IGSParamsToIGSRegistration (which 
#'   results in an implicit inversion by CMTK tools).
#' @details Note that the reference and floating fields have no impact on the 
#'   transformation encoded in the resultant .list folder and can be overridden 
#'   on the command line of CMTK tools.
#' @param x 5x3 matrix of CMTK registration parameters
#' @param reference,floating Path to refererence and floating images.
#' @return list of registration parameters suitable for 
#'   \code{\link{WriteIGSRegistration}}
#' @export
CMTKParamsToCMTKRegistration<-function(x,reference="dummy",floating="dummy"){
  affine_xform=unlist(apply(x,1,list),rec=F)
  names(affine_xform)=c("xlate","rotate","scale","shear","center")
  
  l=list(registration=list(reference_study=reference,
                           floating_study=floating,
                           affine_xform=affine_xform))
  attr(l$registration$reference_study,"quoted")=TRUE
  attr(l$registration$floating_study,"quoted")=TRUE
  version=attr(x,'version')
  if(is.null(version)) version=numeric_version('2.4')
  attr(l,'version')=version
  l
}

#' Convert homogeneous affine registration to CMTK registration list
#' 
#' @details Note that this uses the modern CMTK notation of a "floating" image 
#'   rather than model. It will also mark the resultant registration as version 
#'   "2.4", i.e. using the new correct Compose/Decompose functions of CMTK 
#'   >=2.4.0.
#' @param x 4x4 homogeneous affine matrix
#' @param centre Optional centre of rotation
#' @param reference, floating Optional paths to reference and floating images
#' @return list of CMTK registration parameters
#' @seealso 
#' \code{\link{CMTKParamsToCMTKRegistration},\link{WriteIGSRegistrationFolder}}
#' @export
AffineToCMTKRegistration<-function(x,centre,reference,floating){
  if(!missing(centre)) d=DecomposeAffineToIGSParams(x,centre=centre)
  else d=DecomposeAffineToIGSParams(x)
  arglist<-list(x=d)
  if(!missing(reference)) arglist$reference=reference
  if(!missing(floating)) arglist$floating=floating
  do.call(CMTKParamsToCMTKRegistration, arglist)
}

AffineToIGSRegistration<-function(x,centre,reference,model){
  warning("Use of model_study entry in CMTK registrations is deprecated.",
          " See AffineToCMTKRegistration")
	if(!missing(centre)) d=DecomposeAffineToIGSParams(x,centre=centre)
	else d=DecomposeAffineToIGSParams(x)
	IGSParamsToIGSRegistration(d,reference=reference,model=model)
}

WriteIGSRegistrationFolder<-function(reglist,foldername){
	# Makes a registration folder that could be used as the input to the
	# registration command or warp
	# A transformation in the forward direction (i.e. sample->ref)
	# e.g. as calculated from a set of landmarks where set 1 is the sample
	# is considered an inverse transformation by the IGS software.
	# So in order to use such a transformation as an initial affine with
	# the registration command the switch --initial-inverse must be used
	# specifying the folder name created by this function.
  warning("WriteIGSRegistrationFolder is deprecated in favour of write.cmtkreg")
	dir.create(foldername,showWarnings=FALSE,recursive=TRUE)
	if(!is.list(reglist)) reglist=IGSParamsToIGSRegistration(reglist)
	WriteIGSTypedStream(reglist,file.path(foldername,"registration"))
	
	studysublist= list(studyname=reglist$registration$reference_study)
	attr(studysublist$studyname,"quoted")=TRUE
	studysublist2= studysublist
        if ('model_study' %in% names(reglist$registration)) {
            studysublist2$studyname=reglist$registration$model_study
        } else {
            studysublist2$studyname=reglist$registration$floating_study
        }
	studylist=list(studylist=list(num_sources=2),
		source= studysublist,source=studysublist2)
	
	WriteIGSTypedStream(studylist, file.path(foldername,"studylist"))		
}

IGSLandmarkList<-function(xyzs){
	# IGS Landmark lists are unpaired ie contain information for only 1 brain
	xyzs=data.matrix(xyzs)
	ns=rownames(xyzs)
	ll=list()
	for(i in 1:nrow(xyzs)){
		ll[[i]]=list(name=paste("\"",ns[i],"\"",sep=""),location=xyzs[i,])		
	}
	names(ll)=rep('landmark',length(ll))
	ll
}

WriteIGSLandmarks<-function(xyzs,filename,Force=FALSE){
	ll=IGSLandmarkList(xyzs)
	if(file.exists(filename) && file.info(filename)$isdir) filename=file.path(filename,"landmarks")
	if(file.exists(filename) && !Force) {
		stop(paste(filename,"already exists, use Force=TRUE to replace"))
	}
	WriteIGSTypedStream(ll,filename)
}

ReadIGSLandmarks<-function(...){
	l=ReadIGSTypedStream(...,CheckLabel=FALSE)
	x=t(sapply(l,function(x) x[["location"]]))
	rn=sapply(l,function(x) x[["name"]])
	# nb this is necessary to avoid the names having names 
	# of the form landmarks.1, landmarks.1 ...
	names(rn)<-NULL
	rownames(x)=rn
	x
}

#' Convert a CMTK registration to Amira format (suitable for ResultViewer.hx)
#'
#' Prints matrix row-wise by default which is how Amira expects the file 
#' to look, but strangely enough not how it displays in the console.
#' @cmtkreg Path to registration file or folder 
#' @Transpose Transpose matrix to print out row-wise (default TRUE)
#' @Invert Invert the affine matrix (to go from Sample->Template, default FALSE)
#' @Overwrite Overwrite output file if already exists (Default TRUE)
#' @amiraregfile name of output file (ResultViewer expects hxtransform)
#' @return name of output file
#' @export
#' @seealso \code{\link{ComposeAffineFromIGSParams}}
AmiraRegFromCMTK<-function(cmtkreg,Transpose=TRUE,Invert=FALSE,Overwrite=TRUE,
	amiraregfile='hxtransform'){
	aff=cmtk.dof2mat(cmtkreg)
	if(Invert) aff=solve(aff)
	if(Transpose) aff=t(aff)
	
	# check if we were given registration directory or the actual file
	if(file.info(cmtkreg)$isdir){
		dir = cmtkreg
	}
	else {
		dir=dirname(cmtkreg)
	}
	amirareg=file.path(dir,amiraregfile)
	if(Overwrite || !file.exists(amirareg))
		write(aff,file=amirareg,ncolumns=16)
	amiraregfile
}

#' Convert an Amira registration (4x4 affine) to CMTK format
#'
#' Reads in matrix row-wise by default which is how Steffen Prohaska's
#' Amira ResultViewer.hx script expects the file to look,
#' but strangely enough not how it displays in the console.
#' @amirareg Path to registration file or folder 
#' @cmtkregfolder name of output folder (usually CMTK .list folder)
#' @Transpose Transpose Amira matrix ie read in row-wise (default TRUE)
#' @Invert Invert the affine matrix (to go from Sample->Template, default FALSE)
#' @return name of output file
#' @export
#' @seealso \code{\link{ComposeAffineFromIGSParams}}
CMTKRegFromAmira<-function(amirareg,cmtkregfolder=NULL,Transpose=TRUE,Invert=FALSE){
	if(file.info(amirareg)$isdir) amirareg=file.path(amirareg,"hxtransform")
	if(is.null(cmtkregfolder)) cmtkregfolder=dirname(amirareg)
	
	# backup any existing registration
	cmtkregpath=file.path(cmtkregfolder,"registration")
	if(file.exists(cmtkregpath))
		file.rename(cmtkregpath,paste(cmtkregpath,'bak',sep="."))
	# now make new CMTK format registration
	cmd='mat2dof'
	if(Invert) cmd=paste(cmd,'--inverse')
	if(!Transpose) cmd=paste(cmd,'--transpose')
	cmd=paste(cmd,"--list",shQuote(cmtkregfolder),"<",shQuote(amirareg))
	res=system(cmd)
	if(res!=0)
		stop("Unable to make CMTK registration from amira reg:",amirareg)
	cmtkregpath
}

