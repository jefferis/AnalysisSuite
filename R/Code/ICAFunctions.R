#@+leo-ver=4-thin
#@+node:jefferis.20051014173041.10:@thin R/ICAFunctions.R
#@@language r
# ICAFunctions.R
# Functions that support/adapt use  of AnalyzeFMRI package to do
# Independent Components Analysis
# I am not sure of the exact difference between this implementation
# and that in the fastICA package, but they don't produce identical
# results (and I _prefer_ AnalyzeFRMI)

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

# OK - after investigating the masking etc implemented in AnalyzeFMRI's 
# ICA I think that I have reproduced it fairly well in R code
# that I have put in fastICA.gj and that is my implementation of choice.

# source(file.path(CodeDir,"ICAFunctions.R" ))

#@+others
#@+node:jefferis.20051024092641:WriteAnalyzeFromDensity
WriteAnalyzeFromDensity<-function(d,file.name,pixdim=rep(1.4,3),vox.units='um',...){
	# make array from the list (if required - doesn't touch arrays)
	dd=densityArrayFromList(d)
	
	a=attributes(dd)[c("BoundingBox","x","y","z")]
	if(all(c("BoundingBox") %in% names(a))){ 
		pixdim=diff(matrix(a$BoundingBox,nrow=2))/(dim(dd)[1:3]-1)
	}
	
	# Trim suffix if reqd since f.write.analyze doesn't do this
	file.name<-sub("\\.(hdr|img)","",ignore.case=TRUE,file.name)
	f.write.analyze(dd,file.name,pixdim=pixdim, vox.units=vox.units,...)	
}
#@nonl
#@-node:jefferis.20051024092641:WriteAnalyzeFromDensity
#@+node:jefferis.20051024092641.1:fastICA.gj
fastICA.gj<-function(d,n.comp,alg='par',fun='logcosh',row.norm=T,meth='C',mask,...){
	# my customised fastICA function
	require(fastICA)
	dd<-densityArrayFromList(d)
	origdims=dim(dd)
	ndims<-length(origdims)
	lastdimnames<-dimnames(dd)[[ndims]]
	if(ndims<2) stop("Must have at least 2 dimensions")
	# make a 2d array
	dim(dd)<-c(prod(origdims[-ndims]),origdims[ndims])
	colnames(dd)<-lastdimnames
	
	if(!missing(mask)){
		# remove masked values before ICA
		ddd<-dd[mask>0.5,]
		rval.masked=fastICA(ddd,n.comp=n.comp,alg=alg,fun=fun,row.norm=row.norm,meth=meth,...)
		rval=rval.masked
		# then put them back as 0
		rval$S=matrix(0,nrow=nrow(dd),ncol=n.comp)	
		rval$S[mask>0.5,]<-rval.masked$S
	} else {
		rval=fastICA(dd,n.comp=n.comp,alg=alg,fun=fun,row.norm=row.norm,meth=meth,...)
	}	
	rval$A=t(rval$A)
	row.names(rval$A)=colnames(dd)
	dim(rval$S)<-c(origdims[-ndims],n.comp)
	rval
}
#@-node:jefferis.20051024092641.1:fastICA.gj
#@+node:jefferis.20051024092641.2:ica.gj
ica.gj<-function (d, n.comp, norm.col = TRUE, fun = "logcosh", 
	maxit = 1000, alg.type = "parallel", alpha = 1, tol = 1e-04, 
	mask.file.name = NULL, slices = NULL) 
{
	require(AnalyzeFMRI)
	# my customised version of the ICA function in the AnalyzeFMRI
	# package 
	
	dims=NULL
	if(!is.character(d)){
		useTempFile<-TRUE
		tmp<-tempfile(patt="")
        dd=densityArrayFromList(d)
        
		WriteAnalyzeFromDensity(dd,tmp)
		file.name<-paste(tmp,sep=".","img")
		
	} else {
		useTempFile<-FALSE
		file.name<-d
	}
	hdr <- f.read.analyze.header(file.name)
	dims=hdr$dim[2:5]

	ns <- prod(dims[1:3]) * n.comp
	na <- dims[4] * n.comp
	
	cat("dims=",dims)
	if (length(slices) == 0) 
		slices <- 2:(dims[3] - 1)
	else if (slices == "all") 
		slices <- seq(dims[3])
	else if (any(slices < 1 || slices > dims[3])) {
		return("some of selected slices out of allowable range")
	}

	mask.flag <- 1
	if (length(mask.file.name) == 0) {
		mask.flag <- 0
		mask.file.name <- ""
	}
	col.flag <- 1
	if (norm.col != TRUE) 
		col.flag <- 0
	fun.flag <- 1
	if (fun == "exp") 
		fun.flag <- 2
	def.flag <- 0
	if (alg.type == "deflation") 
		def.flag <- 1
	W <- matrix(rnorm(n.comp * n.comp), n.comp, n.comp)
	a <- .C("ica_fmri_JM", as.character(file.name), as.single(t(W)), 
		as.integer(n.comp), as.integer(1), as.integer(col.flag), 
		as.integer(fun.flag), as.integer(maxit), as.integer(def.flag), 
		as.single(alpha), as.single(tol), as.integer(mask.flag), 
		as.character(mask.file.name), as.integer(slices), as.integer(length(slices)), 
		S = single(ns), A = single(na), PACKAGE = "AnalyzeFMRI")
	S <- array(a$S, dim = c(dims[1:3], n.comp))
    attributes(S)<-c(attributes(S)["dim"],attributes(dd)[c("x","y","z","BoundingBox","class")])
	A <- matrix(a$A, dims[4], n.comp, byrow = TRUE)
    rownames(A)<-dimnames(dd)[[4]]
	# clean up if made temp files
	if(useTempFile) unlink(paste(tmp,sep='.',c('hdr','img')))
	return(list(A = A, S = S, file = file.name, mask = mask.file.name))
}
#@-node:jefferis.20051024092641.2:ica.gj
#@+node:jefferis.20051024092641.3:makeICAMask
makeICAMask<-function(d,thresh=0.1,statfun=max){
	# each pixel of d which > thresh*maximum is set to 1, 0 otherwise
	
	# make array from the list (if required - doesn't touch arrays)
	dd=densityArrayFromList(d)
	subdims=dim(dd)[-length(dim(dd))]
	lastdim=dim(dd)[length(dim(dd))]
	mask<-integer(prod(subdims))
	dim(dd)<-c(prod(subdims),lastdim)
	
	rsd=rowSums(dd)
	absThresh=thresh * statfun(rsd)
	mask[rsd>=absThresh]=1
	dim(mask)<-subdims
	a=attributes(dd)[c("BoundingBox","x","y","z")]
	attributes(mask)=c(attributes(mask),a)
	mask
}
#@nonl
#@-node:jefferis.20051024092641.3:makeICAMask
#@-others
#@nonl
#@-node:jefferis.20051014173041.10:@thin R/ICAFunctions.R
#@-leo
