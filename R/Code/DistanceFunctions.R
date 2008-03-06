# DistanceFunctions.R

# a set of support functions associated with distances and clustering

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

# source(file.path(CodeDir,"DistanceFunctions.R"))

plugindist.list=function (x, distfun, diag = FALSE, upper = FALSE, ...) 
{
	if(!is.list(x)) stop("plugindist.list expects a list of objects")
	distfun=match.fun(distfun)
	if (!is.function(distfun)) 
		stop("inappropriate distance function")
	N <- length(x)
	# distance is given as columns of lower diagonal 
	d=NULL
	for(c in 1:(N-1)){
		for(r in (c+1):N){
			if(c<=N)
				d=c(d,distfun(x[[c]],x[[r]],...))
		}
	}
	attr(d, "Size") <- N
	attr(d, "Labels") <- names(x)
	attr(d, "Diag") <- diag
	attr(d, "Upper") <- upper
#	attr(d, "method") <- METHODS[method]
	attr(d, "call") <- match.call()
	class(d) <- "dist"
	return(d)
}

# this allows subsetting of distance matrices
"[.dist"=function(d,x,y=NULL) {
	if(!is.null(y)) yy=y
	chararg=is.character(x)
	if(!is.null(y) || chararg) {
		if(chararg) y=x
		return(as.dist(as.matrix(d)[x,y]))
	}
	else as.vector(d)[x]
}

matchDistMatrices<-function(d1,d2){
	# save the names of the original dist matrices before we do anything else
	nx=c(substitute(d1),substitute(d2))
	d1=as.matrix(d1); d2<-as.matrix(d2)
	matchingNames=intersect(colnames(d1),colnames(d2))
	d1=d1[matchingNames,matchingNames]
	d2=d2[matchingNames,matchingNames]
	x=list(d1,d2)
	names(x)=nx
	class(x)=c("matchedDistMatrices",class(x))
	invisible(x)
}

print.matchedDistMatrices<-function(x,method='spearman',mantel=F,...){
	sharedGloms=colnames(x[[1]])
	cat(length(sharedGloms)," shared Glomeruli:",sharedGloms,"\n")
	if(isTRUE(mantel)){
		print(mantel(x[[1]],x[[2]],method=method,...))
	} else print(cor.test(as.dist(x[[1]]),as.dist(x[[2]]),method=method))	
	invisible(x)
}

RandomiseDistanceMatrix<-function(m,RandomiseNames=TRUE){
	nrows=nrow(m)
	if(nrows!=ncol(m)) stop("m must be a square matrix")
	s=sample(nrows)
	rval=m[s,s]
	# unfortunately original names are still attached
	if(RandomiseNames && !is.null(dimnames(m))){
		dimnames(rval)=dimnames(m)
	}
	rval
}
