# ShiftNeuron.R
# Simple function to displace a neuron object by specified or random XYZ
# distance

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


ShiftNeuron<-function(ANeuron,shiftfun=rnorm,ShiftIndividualPoints=FALSE,...){
	numPoints=nrow(ANeuron$d)
	if(ShiftIndividualPoints){
		shifts=matrix(shiftfun(3*numPoints,...),ncol=3)
	} else {
		shifts=matrix(shiftfun(3,...),ncol=3,nrow=numPoints,byrow=T)
	}
	ANeuron$d[,c("X","Y","Z")]=ANeuron$d[,c("X","Y","Z")]+shifts
	return(ANeuron)
}
