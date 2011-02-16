#@+leo-ver=4-thin
#@+node:jefferis.20060305214152.1:@thin R/PotentialSynapses.R
#@@language r
# PotentialSynapses.R
# functions to calculate the number of potential synapses between 2
# neurons following the ideas of Stepanyants and Chklovskii

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

# source(file.path(CodeDir,"PotentialSynapses.R"))

# functions to calculate the number of potential synapses between 2
# neurons - expect data to arrive in the form 
# Start x y z , end x y z for each neuron  ie 6 col matrices

# Source ... For MakeStartEndList function
# (nb original version in this file had an error)
source(file.path(CodeDir,"ArbourDensity3D.R")) 
source(file.path(CodeDir,"dist3D_Segment_to_Segment.R")) 

#@+others
#@+node:jefferis.20060305214152.2:PotentialSynapses.neuron
PotentialSynapses.neuron<-function(a,b,bounds,method=c("direct","approx"),...){
	method=match.arg(method)
	if(!inherits(a,"list")){
		if(length(a)>1){
			cat("recursing a ... ")
			return(sapply(a,PotentialSynapses.neuron,b,method=method,...))
		} else a=GetNeuron(a)
	}
	
	if(!inherits(b,"list")){
		if(length(b)>1){
			cat("recursing b ... ")

			return(sapply(b,PotentialSynapses.neuron,a,method=method,...))
		} else b=GetNeuron(b)
	}
	a.sel=MakeStartEndList(a)
	b.sel=MakeStartEndList(b)
	if(!missing(bounds)){
		if(is.character(bounds)) bounds=getBounds(bounds)
		a.sel=restrictToBounds(a.sel,bounds)
		b.sel=restrictToBounds(b.sel,bounds)
	}
	if(method=="direct"){
		DirectPotentialSynapses(a.sel,b.sel,...)
	} else if (method=="approx") {
		PotentialSynapses(a.sel,b.sel,...)
	}
}
#@-node:jefferis.20060305214152.2:PotentialSynapses.neuron
#@+node:jefferis.20060305214152.3:PotentialSynapses
PotentialSynapses<-function(a,b,s=2,sigma=2){
	#Compare for matrices of input data rather than neurons
	
	# short circuit if there are no points to check in one list or other!
	if(nrow(a)*nrow(b)==0) return (0)

	# calculate midpoint position vectors
	ra=(a[,4:6]+a[,1:3])/2
	rb=(b[,4:6]+b[,1:3])/2
	
	# calculate segment  vectors
	na=a[,4:6]-a[,1:3]
	nb=b[,4:6]-b[,1:3]
	
	# calc lengths
	la=sqrt(rowSums(na*na))
	lb=sqrt(rowSums(nb*nb))	
	# deal with any zero length segments
	na=na[la>0,];	nb=nb[lb>0,]
	ra=ra[la>0,];	rb=rb[lb>0,]
	la=la[la>0];  lb=lb[lb>0]
	
	lab=outer(la,lb) # ie has na rows and nb cols
	
	# try replacing this apply by an explicit for loop
	# turns out to be about 30% quicker!
	#sintheta=sin(apply(nb,1,thetaBetween,na))
	# has na rows by nb cols
	sintheta=lab # so just use this to init
	h=nrow(na);w=ncol(nb)
	for(j in 1:h){
		sintheta[j,]=thetaBetween(nb,na[j,])
	}
	
	raminusrb=rowbyrow(ra,rb,"-")
	# find the squared magnitude of the difference vector ra-rb
	l2rab=rowSums(raminusrb*raminusrb)
	#l2rab.old=l2rab
	# reorder the vector
	dim(l2rab)=dim(lab)
	FourSigma2=4*sigma^2
	denom=(4*pi*sigma^2)^(3/2)
	expterm=exp(-l2rab/FourSigma2)/(denom)
	
	rval=2 * s * sum(lab * sintheta * expterm)
	
	#return(list(ra,rb,na,nb,la,lb,lab,sintheta,raminusrb,l2rab.old,l2rab,expterm,rval))
	return(rval)
}
#@+node:jefferis.20060305214152.4:DirectPotentialSynapses
DirectPotentialSynapses<-function(a,b,s=2,maxlineseglength=0.5,returnDistanceList=FALSE){
	# This takes a pair of matrices representing neuron segments
	# and _directly_calculates the number of potential synapses
	# ie regions of approach less than s
    
    #NB It is essential that a and b are formatted as matrices not vectors
    if (!is.matrix(a) || !is.matrix(b)) stop("a and b must be matrices")

	# 2) for each of the segs in a	
	# find the segs in b less than s+maxlineseglength away
    # since we only check the head of each segment with
    # the head of the ref segment we should look to see if they are
    # less than double the max seg length + max accecptable distance
    
    # o-----x      x-----o
    # eg in this arrangement we are measuring from o to o
    
	segswithin<-function(seg,segs,d){
		upperRange=seg[1:3]+d
		lowerRange=seg[1:3]-d
        
        #belowUpperRange <- segs[,1:3] <= upperRange
        #aboveLowerRange <- lowerRange >= segs[,1:3]
        
       # if(nrow(segs)==1) allWithinRange = all(belowUpperRange & aboveLowerRange)
       # else allWithinRange=apply(cbind(belowUpperRange , aboveLowerRange), 1, all)
        
        #return(segs[allWithinRange, ])
        
		 return(segs[segs[,1]<=upperRange[1] & segs[,2]<=upperRange[2]
			 & segs[,3]<=upperRange[3] & segs[,1]>=lowerRange[1]
			 & segs[,2]>=lowerRange[2] & segs[,3]>=lowerRange[3],] )
	}
	selectedBSegs=apply(a,1,segswithin,b,maxlineseglength*2+s*2)
	# This is redundant of course, but ...
	# selectedBSegs=apply(b,1,segswithin,a,maxlineseglength*2+s*2)
	if(length(selectedBSegs)==0)
        return (if(returnDistanceList) list() else 0)
	# This part could fairly easily be rewritten with a sparse matrix
	# rather than a list (which should be 
	distances=list()
	for(i in 1:nrow(a)){
        #cat("i =",i,"  ")
		if(is.matrix(selectedBSegs[[i]]) && nrow(selectedBSegs[[i]])>0){
			distances[[ rownames(a)[i] ]]=
				apply(selectedBSegs[[i]],1,dist3D_Segment_to_Segment,a[i,])
			
		} else if(length(selectedBSegs[[i]])>0){
			# handle case of one segment only	
			distances[[ rownames(a)[i] ]]=
				dist3D_Segment_to_Segment(selectedBSegs[[i]],a[i,])
			# apply handles setting names for vectors
			# but this needs to be done manually when there is only
			# one segment
			names(distances[[ rownames(a)[i] ]])=names(selectedBSegs[i])
		}
	}
    if(returnDistanceList) return(distances)
    else return(sum(unlist(distances)<=s))
}

SimpleDirectPotentialSynapses<-function(a,b,s=2,returnDistanceMatrix=FALSE){
    if (!is.matrix(a) || !is.matrix(b)) stop("a and b must be matrices")
	segNumsWithin<-function(seg,segs,d){
		upperRange=seg[1:3]+d
		lowerRange=seg[1:3]-d
        return(which(segs[,1]<=upperRange[1] & segs[,2]<=upperRange[2]
			 & segs[,3]<=upperRange[3] & segs[,1]>=lowerRange[1]
			 & segs[,2]>=lowerRange[2] & segs[,3]>=lowerRange[3]) )
	}
	segswithin<-function(seg,segs,d){
		upperRange=seg[1:3]+d
		lowerRange=seg[1:3]-d
        
        #belowUpperRange <- segs[,1:3] <= upperRange
        #aboveLowerRange <- lowerRange >= segs[,1:3]
        
       # if(nrow(segs)==1) allWithinRange = all(belowUpperRange & aboveLowerRange)
       # else allWithinRange=apply(cbind(belowUpperRange , aboveLowerRange), 1, all)
        
        #return(segs[allWithinRange, ])
        
		 return(segs[segs[,1]<=upperRange[1] & segs[,2]<=upperRange[2]
			 & segs[,3]<=upperRange[3] & segs[,1]>=lowerRange[1]
			 & segs[,2]>=lowerRange[2] & segs[,3]>=lowerRange[3],] )
	}

    selectedBSegs=apply(a,1,segswithin,b,maxlineseglength*2+s*2)
    
    rval=matrix(NA,nrow=nrow(a),ncol=nrow(b))
    rownames(rval)=rownames(a)
    colnames(rval)=rownames(b)

    for( rowOfA in seq(len=nrow(a)) ){
        rowsOfB=selectedBSegs[[rowOfA]]
        if(length(rowsOfB)>0){
            rval[rowOfA,rowsOfB]=apply(b[rowsOfB,],1,dist3D_Segment_to_Segment,a[rowOfA,])
        }
    }
    if(returnDistanceMatrix) return(rval)
    else return(sum(rval<=s))
}
#@-node:jefferis.20060305214152.4:DirectPotentialSynapses
#@-node:jefferis.20060305214152.3:PotentialSynapses
#@+node:jefferis.20060305214152.5:util functions
restrictToBounds<-function(a,bounds){
	a[  a[,1]>=bounds[1] & a[,1]<=bounds[2] &
	    a[,2]>=bounds[3] & a[,2]<=bounds[4] &
	    a[,3]>=bounds[5] & a[,3]<=bounds[6], ]
}

dotprod=function(a,b){
# expects 2 matrices with n cols each
	c=a*b
	if(length(dim(c))>1) 	rowSums(c)
	else sum(c)
}

rowbyrow<-function(X,Y,FUN="-",...){
	if(ncol(X)!=ncol(Y)) return(NA)
	FUN=match.fun(FUN)
	rX=nrow(X)
	rY=nrow(Y)
	#dX=dim(X)
	#dY=dim(Y)
	X=matrix(t(X),nrow=rX*rY,ncol=ncol(X),byrow=T)
	Y <- matrix(rep(Y, rep.int(rX, length(Y))),ncol=ncol(X))
	#Y=matrix(Y,nrow=rX*rY,ncol=ncol(Y),byrow=T)
	FUN(X,Y,...)
	#list(X=X,Y=Y)
}

normbyrow=function(a){
	# returns euclidean norm (by row if reqd)
	c=a*a
	if(length(dim(c))>1) 	sqrt(rowSums(c))
	else sqrt(sum(c))
}

sinthetafn<-function(a,b){
	if(is.vector(a) && is.matrix(b)){
		a=matrix(a,nrow=nrow(b),ncol=length(a),byrow=T)
	} else if(is.vector(b) && is.matrix(a)){
		b=matrix(b,nrow=nrow(a),ncol=length(b),byrow=T)
	}
	costheta= zapsmall(dotprod(a,b)/(normbyrow(a)*normbyrow(b)))
	sqrt(1-costheta*costheta)
}

# note use of zapsmall in case of rounding error
thetaBetween=function(a,b) {
	if(is.vector(a) && is.matrix(b)){
		a=matrix(a,nrow=nrow(b),ncol=length(a),byrow=T)
	} else if(is.vector(b) && is.matrix(a)){
		b=matrix(b,nrow=nrow(a),ncol=length(b),byrow=T)
	}
	acos(zapsmall(dotprod(a,b)/(normbyrow(a)*normbyrow(b))))
}
#@-node:jefferis.20060305214152.5:util functions
#@-others
#@-node:jefferis.20060305214152.1:@thin R/PotentialSynapses.R
#@-leo
