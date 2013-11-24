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

# source(file.path(CodeDir,"ReadIGSFiles.R"))


trim =function(t) sub('[[:space:]]+$', '', sub('^[[:space:]]+', '', t))

ReadIGSRegistration <- function (filename,ReturnRegistrationOnly=TRUE){
	# ReturnRegistrationOnly = FALSE will not attempt to extract the registration element from the
	# registration file
	if(!file.exists(filename)) {
		warning("filename ",filename," does not exist. ReadIGSRegistration Returning NULL")
		return(NULL)
	}
	if(file.info(filename)$isdir){
		# this is a directory, so see if we can find the registration
		dirname=filename
		filename=dir(dirname,patt="^registration(\\.gz){0,1}",full.names=T)[1]
		if(is.na(filename)) 
			stop(paste("Unable to read registration file in",dirname))
	}
	r=ReadIGSTypedStream(filename)
	if(!is.null(r$registration) && ReturnRegistrationOnly) return(r$registration)
	else return(r)
}

ReadIGSTypedStream<-function(con, CheckLabel=TRUE){
#  Reads Torsten's IGS TypedStream format which is what he uses for:
#  registration
#  studylist
#  images
#  Note that there are special methods to handle the 
#  coefficients and active members of a spline warp 

	if(is.character(con)) {
		filename=con
		if(any(grep("\\.gz$",filename))){
			con=gzfile(filename,'rb')
		} else con=file(con,'rb')

		t=readLines(con,1)
		if( !any(grep("! TYPEDSTREAM",t[1],fixed=TRUE)) ) 
			stop(paste("This doesn't appear to be an IGS TypedStream:",filename))
	}
	
	l=list()
	
	checkLabel=function(label) 	{
		if( any(names(l)==label)  ){
			newlabel=make.unique(c(names(l),label))[length(l)+1]
			warning(paste("Duplicate item",label,"renamed",newlabel))
			label=newlabel
		}
		label
	}
	# Should this check to see if the connection still exists?
	# in case we want to bail out sooner
	while ( isTRUE(isOpen(con)) ){
		thisLine<-readLines(con,1)
		# no lines returned - ie end of file
		if(length(thisLine)==0) break

		# trim and split it up by white space
		thisLine=trim(thisLine)
		
		# skip if this is a blank line
		if(nchar(thisLine)==0) next

		items=strsplit(thisLine," ",fixed=TRUE)[[1]]
		
		if(length(items)==0) next
		# get the label and items
		label=items[1]; items=items[-1]
		#cat("\nlabel=",label)
		#cat("; items=",items)

		# return list if this is the end of a section
		if(label=="}") {
			#cat("end of section - leaving this recursion\n")
			return (l)
		}		
		if(items[1]=="{"){
			# parse new subsection
			#cat("new subsection -> recursion\n")
			# set the list element!
			if(CheckLabel)
				label=checkLabel(label)

			l[[length(l)+1]]=ReadIGSTypedStream(con,CheckLabel=CheckLabel)
			names(l)[length(l)]<-label
			next
		}
		
		if(label == "coefficients"){
			# process coefficients
			numItems=prod(l[["dims"]])*3
			#cat("numItemsToRead =",numItems)
			remainingItems=scan(con,n=numItems-3,quiet=TRUE)
			l[[label]]=matrix(c(as.numeric(items),remainingItems),ncol=3,byrow=T)
			
		} else if (label == "active"){
			# process active flags
			numItems=prod(l[["dims"]])*3
			#cat("numItemsToRead =",numItems)
			# nb floor since we have read one line of 30 already
			numLinesToRead=floor(numItems/30)
			#cat("numLinesToRead=",numLinesToRead)
			x=c(items,trim(readLines(con,numLinesToRead,ok=FALSE)))
			#x=paste(x,collapse="")
			bits=strsplit(x,"")
			bits=as.raw(unlist(bits))
			l[[label]]=bits
		} else {
			# ordinary item
			# Check first item
			firstItemFirstChar=substr(items[1],1,1)		
			if(any(firstItemFirstChar==c("-",as.character(0:9)) )){
				# convert to numeric if not a string
				items=as.numeric(items)
			} else if (firstItemFirstChar=="\""){
				# dequote quoted string
				# can do this by using a textConnection
				tc=textConnection(thisLine)
				items=scan(tc,what="",quiet=TRUE)[-1]
				close(tc)
				attr(items,"quoted")=TRUE
			}
			# check if the list already has one of these
			
			# set the list element!
			if(CheckLabel)
				label=checkLabel(label)

			l[[length(l)+1]]=items
			names(l)[length(l)]<-label
		}
	}
	# we should only get here once if we parse a valid hierarchy
	close(con)
	if(isTRUE(try(file.exists(filename)))){
		attr(l,"file.info")=file.info(filename)
	}		
	return(l)
}

WriteIGSTypedStream<-function(l,filename,gzip=FALSE){
	# Will take a list in the form returned by ReadIGSTypedStream and
	# write it out to a text file
#	con=if(gzip) file(filename,'w') else gzfile(filename,'w')
	con=file(filename,'w')
	cat("! TYPEDSTREAM 1.1\n\n",file=con)
	.WriteIGSTypedStream.list(l,con)
	# iterate over list 
	close(con)
}

.WriteIGSTypedStream.list<-function(x,con,tablevel=0){
	nn <- names(x)
	ll <- length(x)
	tabs=""
	for(i in seq(len=tablevel)) tabs=paste(tabs,sep="","\t")
	if (length(nn) != ll) 
		nn <- paste("Component", seq(ll))
	for (i in seq(length = ll)) {
		# cat("i=",i,"name=",nn[i],"mode=",mode(x[[i]]),"\n")
		if (is.list(x[[i]])) {
			cat(sep="",tabs,nn[i]," {\n",file=con,append=TRUE)
			.WriteIGSTypedStream.list(x[[i]],con,tablevel+1)
			cat(sep="",tabs,"}\n",file=con,append=TRUE)
		} else if(is.matrix(x[[i]])){
			cat(sep="",tabs,nn[i]," ")
			write(t(x[[i]]),file=con,append=TRUE,ncolumns=ncol(x[[i]]))
		} else {
			qsep=""
			if(isTRUE(attr(x[[i]],"quoted"))) qsep="\""
			thisline=paste(nn[i],sep=""," ",paste(qsep,x[[i]],qsep,sep="",collapse=" "),"\n")
			cat(sep="",tabs,thisline,file=con,append=TRUE)
		}
	}
}

IGSParamsToIGSRegistration<-function(x,reference="dummy",model="dummy"){
	# Note that this produces a list that "should" be identical to what
	# would be produced by reading in an IGS registration file
	# This could now be written out wiith WriteIGSTypedStream
	affine_xform=unlist(apply(x,1,list),rec=F)
	names(affine_xform)=c("xlate","rotate","scale","shear","center")
	
	
	l=list(registration=list(reference_study=reference,
						model_study=model,
						affine_xform=affine_xform))
	attr(l$registration$reference_study,"quoted")=TRUE
	attr(l$registration$model_study,"quoted")=TRUE
	l
}

AffineToIGSRegistration<-function(x,centre,reference,model){
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

#' Convert CMTK registration parameters into homogeneous affine matrix
#' 
#' @param params 5 x 3 matrix of CMTK registration parameters
#' Note that this matrix will have the 4th row 0 0 0 1
HomogenousAffineFromCMTKParams<-function(params){
  tf<-tempfile(fileext='.list')
  on.exit(unlink(tf,recursive=TRUE))
  WriteIGSRegistrationFolder(params,foldername=tf)
  cmtk.dof2mat(tf,Transpose=TRUE)
}

#' Read CMTK registration with dof2mat and convert to homogeneous affine matrix
#' 
#' @details Transpose is true by default since this results in the orientation
#'   of cmtk output files matching the orientation in R.
#' @param reg Path to input registration file
#' @param Transpose ouput matrix so that form on disk matches R's convention
cmtk.dof2mat<-function(reg,Transpose=TRUE){
  cmd="dof2mat"
  if(Transpose) cmd=paste(cmd,'--transpose')
  cmd=paste(cmd,shQuote(reg))
  rval=system(cmd,intern=TRUE)
  numbers=as.numeric(unlist(strsplit(rval,"\t")))
  mat=matrix(numbers,ncol=4,byrow=TRUE)
}

#' Use CMTK mat2dof to convert homogeneous affine matrix into CMTK registration
#' 
#' @details If no output file is supplied, 5x3 params matrix will be returned 
#'   directly. Otherwise a logical will be returned indicating success or 
#'   failure at writing to disk.
#' @details Transpose is true by default since this results in an R matrix with 
#'   the transpose in the fourth column being correctly interpreted by cmtk.
#' @param m Homogenous affine matrix (4x4) last row 0 0 0 1 etc
#' @param f Output file (optional)
#' @param Transpose the input matrix so that it is read in as it appears on disk
#' @return 5x3 matrix of CMTK registration parameters or logical
cmtk.mat2dof<-function(m, f=NULL, Transpose=TRUE){
  if(!is.matrix(m) || nrow(m)!=4 || ncol(m)!=4) stop("Please give me a homogeneous affine matrix (4x4)")
  inf=tempfile()
  on.exit(unlink(inf),add=TRUE)
  
  write.table(m, file=inf, sep='\t', row.names=F, col.names=F)
  # always transpose because mat2dof appears to read the matrix with last column being 0 0 0 1
  cmd="mat2dof"
  if(Transpose) paste(cmd,'--transpose')
  if(is.null(f)){
    cmd=paste(cmd,sep="<",shQuote(inf))
    params=read.table(text=system(cmd,intern=T),sep='\t',comment.char="")[,2]
    if(length(params)!=15) stop("Trouble reading mat2dof response")
    return(matrix(params,ncol=3,byrow=TRUE))
  } else {
    cmd=paste(cmd,'--list',shQuote(f),"<",shQuote(inf))
    return(system(cmd)==0)
  }
}
