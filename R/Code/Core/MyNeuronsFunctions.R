# MyNeuronsFunctions.R
# Utility functions associated with MyNeurons list
# Begun 2-March-2005

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

# source(file.path(CodeDir,"MyNeuronsFunctions.R"))


SortMyNeurons<-function(){
	# sort neurons in MyNeurons by brain name
	NeuronNames=getBrain(sapply(MyNeurons,function(x) x$NeuronName))
	MyNeurons<<-MyNeurons[order(NeuronNames)]
}

# This handy function removes the extra levels from all of the columns in
# a data frame which are factors.  This is useful because when dataframes
# are subsetted, factor levels are retained from the original dataframe
# even if the corresponding rows have been dropped.
# Originally from QDG Code Base
DropMissingFactorLevels<-function(df){
	.Deprecated('droplevels')
	droplevels(df)
}

MakeMNInfo<-function(){
	# Make MNInfo data frame describing contents of MyNeurons list
	#tMNInfo=MNInfo # commented out by GJ 2005-09-22 why reqd?
	NeuronNames=getBrain(sapply(MyNeurons,function(x) x$NeuronName))
	MyNeurons<<-MyNeurons[order(NeuronNames)]
	tMNInfo<-subset(Info,Brain%in%NeuronNames)
	tMNInfo<-tMNInfo[order(tMNInfo$Brain),]
	tMNInfo<-DropMissingFactorLevels(tMNInfo)

	# A bit of juggling so that we can get a seq for each glomerulus
	# Can use this to limit processing to first 5 etc.
	bb=by(tMNInfo$Glomerulus,tMNInfo$Glomerulus,function(x) seq(length(x)))
	tMNInfo$GlomSeq<-unsplit(bb,tMNInfo$Glomerulus)

	tMNInfo$NumNAs<-sapply(MyNeurons,function(x) sum(is.na(x$d$X)))
	tMNInfo$MBP1<-sapply(MyNeurons,function(x) ifelse(is.null(x$MBPoints[1]),NA,x$MBPoints[1]))
	tMNInfo$MBP2<-sapply(MyNeurons,function(x) ifelse(is.null(x$MBPoints[2]),NA,x$MBPoints[2]))
	tMNInfo$LHBP<-sapply(MyNeurons,function(x) ifelse(is.null(x$LHBranchPoint),NA,x$LHBranchPoint))
	tMNInfo$PNType<-PNType(tMNInfo$Glomerulus)
	tMNInfo$Seq<-seq(MyNeurons)
	tMNInfo$TraceFile<-sapply(MyNeurons,function(x) x$InputFileName)
	tMNInfo$nTrees<-sapply(MyNeurons,function(x) ifelse(is.null(x$nTrees),1,x$nTrees))
	tMNInfo$StartPoint<-sapply(MyNeurons,"[[","StartPoint")
	tMNInfo$CreatedAt<-sapply(MyNeurons,function(x) x$CreatedAt)
	class(tMNInfo$CreatedAt)<-"POSIXct" # gets dropped for some reason
	attr(tMNInfo,"CreatedAt")<-Sys.time()
	attr(tMNInfo,"CreatedOn")<-HostName
	MNInfo<<-tMNInfo
}

PNType<-function(x){

	if(is.numeric(x)) x=as.character(MNInfo$Glomerulus[x])
	if(is.factor(x)) x=as.character(x)
	if(is.list(x)) x=x$CellType
	else if(length(x)>1) return(sapply(x,PNType))
	else if(is.character(x)){
		if(x%in%MNInfo$Brain) x=MNInfo$Glomerulus[MNInfo$Brain==x]
	}
	if(length(grep("(NP|MZ|acj6)",x,ignore.case=TRUE))>0){
		return("LHN")
	}
	else if(length(grep("^v",x))>0){
		return("mPN")
	}
	else #if(x%in%unique(MNInfo$Glomerulus)){
		return("iPN")
	# commented out because there is a failure if MNInfo doesn't exist
	
	#} else {
	#	return("NA")
	#}
}

UpdateKeyPoints<-function(NeuronsToUpdate=NULL,NeuronTypes=c("iPN","mPN"),...){
	if(is.null(NeuronsToUpdate)){
		# First figure out which neurons do not have any key point data
		NeuronsToUpdate=which(sapply(MyNeurons,function(x) is.null(x$LHBranchPoint)))
	}
	NeuronTypes.MN=sapply(MyNeurons[NeuronsToUpdate],PNType)
	NeuronsToUpdate=NeuronsToUpdate[NeuronTypes.MN%in%NeuronTypes]
	if(length(NeuronsToUpdate)==0) return
	# then process them
	for( i in NeuronsToUpdate){
		if (PNType(MyNeurons[[i]])=="iPN") Components=c("Axon","MB","LH") else
			if (PNType(MyNeurons[[i]])=="mPN") Components=c("Axon","LH") else
			if (PNType(MyNeurons[[i]])=="LHN") Components="LH"
		cat(i," ")
		MyNeurons[[i]]<<-Short.GetKeyPoints(MyNeurons[[i]],Components=Components,...)
		MyNeurons[[i]]<<-PartitionNeuron(MyNeurons[[i]],Components=Components)
	}	
}

UpdateSegLengths<-function(NeuronsToUpdate=NULL,NeuronList=MyNeurons){
	if(is.null(NeuronsToUpdate)){
		# First figure out which neurons do not have any key point data
		NeuronsToUpdate=which(sapply(NeuronList,function(x) is.null(x$SegLengths)))
	} else if(NeuronsToUpdate=='all') NeuronsToUpdate=seq(NeuronList)


	for(i in NeuronsToUpdate){
		NeuronList[[i]]$SegLengths<-SegLengths(NeuronList[[i]])
	}
	NeuronList
}

MergeNeuron<-function(src,dest,dryrun=FALSE,verbose=T){
	# merge additional fields of src neuron into dest
	missingFields=setdiff(names(src),names(dest))
	if(length(missingFields)>0){
		if(verbose) cat(ifelse(dryrun,"Dry Run",""),"Merging fields",missingFields,"for",dest$NeuronName,"\n")
		if(!dryrun){
			oldClass=class(dest)
			dest=c(dest,src[missingFields])
			class(dest)=oldClass
		}
	} else if(verbose) cat("No fields to merge")
	invisible(dest)
}

SaveMyNeurons<-function(filestem="MyNeurons",fileext=".rda",
		list=c("MyNeurons", "MNInfo","InterpolatedNeurons"),Backup=TRUE){
		rdaFile=file.path(ObjDir,paste(filestem,fileext,sep=""))
		if(Backup && file.exists(rdaFile)){
				bakFile=file.path(ObjDir,paste(filestem,".bak",sep=""))
				file.rename(rdaFile,bakFile) # nb overwrites
		}
		save(list=list, file = rdaFile)
}

UpdateG5<-function(filestem="MyNeurons",fileext=".rda"){
		rdaFile=file.path(ObjDir,paste(filestem,fileext,sep=""))
		cmd=paste("rsync -Cavuz", rdaFile,"gjg5:projects/PN2/analysis/CombinedObjs")
		cat("Running:",cmd,"\n")
		system(cmd)
}
