# QhullInterface.R
##################
# functions to prepare make a call to qhull
# via David Marchette's interface (only on unix)
# FindLHNeuronVolume can be used to find the volume enclose by the neurons's
# branching pattern in the LH

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

# source(file.path(CodeDir,"QhullInterface.R"))

source(file.path(CodeDir,"MarchetteQhullInterface.R"))

# FindLHNeuronVolumes
FindLHNeuronVolumes<-function(ANeuron){
    LHPoints<-unique(unlist(sapply(ANeuron$LHSegNos,function(x) {ANeuron$SegList[[x]]})))
    LHPoints.XYZ<-ANeuron$d[LHPoints,c("X","Y","Z")]
    LHEndPointsIdxs<-match(ANeuron$EndPoints,LHPoints,nomatch=0)
    LHEndPoints.XYZ<-LHPoints.XYZ[LHEndPointsIdxs,]
    
    # USE MY INTERFACE
    ##TempHullInFile<-WriteHullFileFromData(LHPoints.XYZ)
    # OR USE MARCHETTE'S
    #qhull(LHPoints.XYZ,total=T,silent=F)
    FullSAV<-qhull.va(LHPoints.XYZ)
    TipSAV<-qhull.va(LHEndPoints.XYZ)
    return(list(FullArea=FullSAV$SurfArea,FullVolume=FullSAV$Volume,
	    TipArea=TipSAV$SurfArea,TipVolume=TipSAV$Volume))
}


WriteHullFileFromData<-function(XYZData,hullfile=tempfile("hull")){     
    hullfile<-basename(hullfile)
    file.create(hullfile)
    cat(length(LHPoints.XYZ[,1]),"\n",file=hullfile)
    write.table(LHPoints.XYZ,file=hullfile,
    row.names=F,col.names=F,append=T)
    return(hullfile)
}

FindNeuronPointSetVolume<-function(Aneuron,PointSet){
    MyPoints.XYZ<-ANeuron$d[PointSet,c("X","Y","Z")]
    # USE MY INTERFACE
    ##TempHullInFile<-WriteHullFileFromData(LHPoints.XYZ)
    # OR USE MARCHETTE'S
    #qhull(LHPoints.XYZ,total=T,silent=F)
    qhull.va(MyPoints.XYZ)
}

FindAllContourVols<-function(recs=seq(MyNeurons)){
    CellTypes<-GetCellType(recs)
    df<-data.frame(CellType=CellTypes,matrix(0,length(recs),4))
    names(df)<-c("CellType","FullArea","FullVolume","TipArea","TipVolume")
    cat("Calling qhull to calculate convex hulls ")
    for(i in recs){
	SAV<-FindLHNeuronVolumes(MyNeurons[[i]])
	df$FullArea[i]<-SAV$FullArea
	df$FullVolume[i]<-SAV$FullVolume
	df$TipArea[i]<-SAV$TipArea
	df$TipVolume[i]<-SAV$TipVolume
	cat(".")
    }
    cat("\n")
    return(df)
}

CallFindAllLHVols<-function(){
    df2<-FindAllLHVols()
    oldwd<-getwd()
    setwd(file.path(HomeDir,"outfiles"))
    write.table(df2,file="NeuronLHVolumes.txt")
}
