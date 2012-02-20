# source(file.path(CodeDir,"ArbourDensity3D.R"))

# Routines to handle production of 3D arbour density data 
# from neuronal traces - see also PotentialSynapses.R


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

ExtractAllSegMidPoints=function(gloms=levels(MNInfo$Glomerulus),neuronData=MyNeurons){
	z=do.call(rbind,lapply(neuronData[MNInfo$Glomerulus%in%gloms],function(x) {
		#data.frame( ExtractSegMidPointsWeights(x),Glomerulus=x$CellType,NeuronName=x$NeuronName )
	 ExtractSegMidPointsWeights(x) 
		}))
	as.data.frame(z,optional=TRUE)
}

ExtractSegMidPointsWeights<-function(ANeuron,mask=seq(len=length(ANeuron$SegList))){
		# For each individidual (seglet) segment in the neuron, find its midpoint 
		# and its length

		sel=MakeStartEndList(ANeuron,mask)
		midpoints=(sel[,4:6]+sel[,1:3])/2
		# calculate segment  vectors
		nsel=sel[,4:6]-sel[,1:3]
		# calc lengths
		lsel=sqrt(rowSums(nsel*nsel))
		rval=cbind(midpoints,lsel)
		colnames(rval)=c("X","Y","Z","L")
		rval
}

MakeStartEndList<-function(ANeuron,mask=seq(len=length(ANeuron$SegList))){
		# make an ordered list of the start/end points for each seg in form
		# Start x y z , end x y z   ie 6 col matrix

		StartIdxs=unlist(sapply(ANeuron$SegList[mask],function(x) x[-length(x)]))
		EndIdxs=unlist(sapply(ANeuron$SegList[mask],function(x) x[-1]))

		d=data.matrix(ANeuron$d[,c("X","Y","Z")])
		cbind(d[StartIdxs,],d[EndIdxs,])		
}

fracBelow<-function(x,below=0.490,comparator='<') {
		# little function to check what fraction of the segments lengths of a
		# neuron are below a specified value
		e=ExtractSegMidPointsWeights(x)
		FUN=match.fun(comparator)
		sum(FUN(e[,4],below))/nrow(e)
}

