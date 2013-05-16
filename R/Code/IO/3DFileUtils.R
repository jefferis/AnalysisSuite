# 3DFileUtils.R
# #############
# These routines write out a neuron into a POV file
# suitable for use with the POV ray tracer
# They also include a routine to write out a "Borst" file
# the format used by Axel Borst (eg in his TanBase) - I think
# it must orginally be due to Frederic Theunissen
# Finally there is also a routine to convert a SWC file directly to
# POV format
# ########
# depends on NeuronFunctions4.R
# GSXEJ 020301

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

# source(file.path(CodeDir,"3DFileUtils.R"))

# This works except for the fact that it writes "C" instead of C in the outfile
# OK - fixed that problem - I needed the quote=F option in write.table
# to prevent quoting.  I found the answer in the R Data Import/Export Manual.
WriteBorstFile<-function(ANeuron,AFileName=""){

    if (is.character(ANeuron)){
	ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
    }
    if (is.numeric(ANeuron)){
	ANeuron<-MyNeurons[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
	warning("Cannot understand passed neuron")
	return(F)
    }

    if(AFileName=="") AFileName<-paste(strsplit(ANeuron$NeuronName,'[.]')[[1]][1],sep="",".abf")
    
    if(!file.create(AFileName)) stop(paste("Couldn\'t create file",AFileName))
    
    #BorstArray<-ANeuron$d[,c("PointNo","X","Y","Z","W")]
    # Set each point to C, B or T for connect?, branch, terminal
    # The I is to avoid having the NodeType turned into a factor
    # Hmm removed the I() call and tried with the thing as a factor
    BorstArray<-cbind(PointNo=ANeuron$d$PointNo, NodeType=factor(rep('C',length(ANeuron$d[,1])),levels=c('B','C','T')),ANeuron$d[,c("X","Y","Z","W")])
    BorstArray$NodeType[ANeuron$BranchPoints]<-"B"
    BorstArray$NodeType[ANeuron$EndPoints]<-"T"
    BorstArray$NodeType<-factor(BorstArray$NodeType)
    
    # Should apparently count from 0 not 1
    BorstArray$PointNo<-BorstArray$PointNo-1
    
    write.table(BorstArray,AFileName,row.names = FALSE, col.names = FALSE,quote=F)
    
}

# TOFIX 6 Aug 2005
ReadBorstFile<-function(filename){
	# Read a file in axel borst's ".abf" file format
	# Probably do this by constructing an SWC array and then parsing that
	
	
# 	1  B  153.69  7.07097  -122.5  7.5
# 	2  C  155.67  6.6467  -122.5  7.5
# 	3  C  158.71  4.73751  -118.2  9.5
# 	4  B  161.185  4.38396  -117.4  9.5
# 	5  C  162.387  8.98015  -118  5
# 	6  B  163.236  9.54583  -118  5
# 	7  B  164.084  6.15172  -120.3  0.5
# 	8  T  163.095  5.58604  -120.3  0.5

	bdf=read.table(filename)
	colnames(bdf)<-c("PointNo","PointType","X","Y","Z","W")
	# now need to figure out parents - this will be based on the B, C or T fields
	# C and T nodes are easy - parent is just the preceding node
	#bdf$ParentNo[1]=-1
	#bdf$ParentNo[bdf$PointType!="B"]=bdf$PointNo[bdf$PointType!="B"]-1
	# Hmm but what about the nodes after T nodes - actually they must refer to
	# the branch before 

	
	SegList=parseBorstArray(bdf)
	ParsedNeuron<-list(NeuronName=NeuronNameFromFileName(filename),
		InputFileName=filename,
		CreatedAt=Sys.time(),
		NodeName=Sys.info()["nodename"],
		InputFileStat=file.info(filename)[1,],
		InputFileMD5=md5sum(path.expand(filename)),
		NumPoints=nrow(bdf),
		StartPoint=1, # NB always 1 ?
		BranchPoints=which(bdf$PointType=="B"),
		EndPoints=which(bdf$PointType=="T"),
		NumSegs=length(SegList),
		SegList=SegList,
		d=bdf	)
	return(ParsedNeuron)

	
	
}

parseBorstArray<-function(bdf){
	# iterate over array assuming that
	# all bifurcations are binary
	
	seglist=list()
	si=0
	openBranches=NULL
	processedBranches=NULL
	for(i in 1:nrow(bdf)){
		cat(i,".")
		if(bdf$PointType[i]=="B"){
			# if a branch, then close seg if one is open,
			if(si>0) {
				seglist[[si]]=c(seglist[[si]],i)
				lpb=length(processedBranches)
				processedBranches[lpb]=processedBranches[lpb]+1
			# pop if processedBranches=2
				if(processedBranches[lpb]>1){
					processedBranches=processedBranches[-lpb]
					openBranches=openBranches[-1]
				}
			} else {
				# anything to do if this was the root point?
			}
			# open a branch, 
			openBranches=c(openBranches,i)
			# set corresponding processedBranches to 0, 
			processedBranches=c(processedBranches,0)
			# make a new seg
			si=si+1
			seglist[[si]]=i
		}
		
		
		if(bdf$PointType[i]=="C"){
			# is this the root?
			if(si==0){
				si=si+1
				seglist[[si]]=i
				processedBranches=i
				openBranches=1
			} else {
				# regular continuation, so add to current seg
				seglist[[si]]=c(seglist[[si]],i)
			}			
		}
		
		if(bdf$PointType[i]=="T"){
		# if a t then close seg,
			seglist[[si]]=c(seglist[[si]],i)
		# inc processedBranches,
			lpb=length(processedBranches)
			processedBranches[lpb]=processedBranches[lpb]+1
		# pop if processedBranches=2
			if(processedBranches[lpb]>1){
				processedBranches=processedBranches[-lpb]
				openBranches=openBranches[-1]
			}
			
			
		# make a new seg if there are further points, adding the parental branchpoint
			if((i+1)<=nrow(bdf)){
				si=si+1
				lpb=length(processedBranches)
				processedBranches[lpb]=processedBranches[lpb]+1
				seglist[[si]]=openBranches[lpb]
				# pop if processedBranches=2
				if(processedBranches[lpb]>1){
					processedBranches=processedBranches[-lpb]
					openBranches=openBranches[-1]
				}

			}
		}
		
	}
	return(seglist)
}

	
BlendNeuron<-function(PreNeuron,PostNeuron,BlendFraction=0.5){
    # Function to return a new neuron which is a linear
    # interpolation between the pre and post transformation of
    # a single original tracing
    # if BlendFraction = 1 , it's all post if 0, all pre
    if(PreNeuron$NeuronName!=PostNeuron$NeuronName){
	warning("This function only works for different transformations of the same neuron")
	return(NULL)
    }
    if(BlendFraction<0 | BlendFraction>1){
	warning("BlendFraction must be between 0 (all pre) and 1 (all post)")
	return(NULL)
    }
    
    # Blend the Neuron
    RNeuron=PreNeuron
    if (!is.null(PreNeuron$d) & !is.null(PostNeuron$d)){
	RNeuron$d[,c("X","Y","Z","W")]=PreNeuron$d[,c("X","Y","Z","W")]+
	    BlendFraction*(PostNeuron$d[,c("X","Y","Z","W")]-PreNeuron$d[,c("X","Y","Z","W")])
    }
	
    # Blend the MB
    if (!is.null(PreNeuron$MB$d) & !is.null(PostNeuron$MB$d)){
	RNeuron$MB$d[,c("X","Y","Z")]=PreNeuron$MB$d[,c("X","Y","Z")]+
	    BlendFraction*(PostNeuron$MB$d[,c("X","Y","Z")]-PreNeuron$MB$d[,c("X","Y","Z")])
    }

    # Blend the LH
    if (!is.null(PreNeuron$LH$d) & !is.null(PostNeuron$LH$d)){
	RNeuron$LH$d[,c("X","Y","Z")]=PreNeuron$LH$d[,c("X","Y","Z")]+
	    BlendFraction*(PostNeuron$LH$d[,c("X","Y","Z")]-PreNeuron$LH$d[,c("X","Y","Z")])
    }

    return (RNeuron)
}

WriteBlendedPOVFiles<-function(PreNeurons,PostNeurons,AFileName="BlendedNeurons.pov",WithBlob=F,WithCone=TRUE,FieldOfView=100,Colour=3,InMem=F,...){
    if (!is.list(PreNeurons)| !is.list(PostNeurons)){
	warning("Cannot understand passed lists of neurons")
	return(F)
    }

    # OK Lets add a header
    cat(paste("//",AFileName),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=AFileName,append=F)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=AFileName,append=TRUE)
    MinAxes=c(0,0,0)
    MaxAxes=c(168,168,88)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)
    # Don't think I want to shift the origin
    NeuronCentre=c(0,0,0)

    # Set the camera
    CameraLoc=c(84,84,-80)
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(100,84,44),collapse=",")
    TCamera<-paste("camera { perspective orthographic location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
    cat(TCamera,"\n\n",file=AFileName,append=TRUE)
    
    # AND FINALLY THE LIGHT
	makeLight<-function(LightPos){
		TLightPos<-paste(LightPos,collapse=",")
		TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
		cat(TLight,"\n\n",file=AFileName,append=TRUE)
	}
	
	makeLight(c(84,84,-100))
	makeLight(c(168,84,-100))
    
    # This can be used as a way of moving the neurons for only part of the time
#	cat("#declare nrnClock=(1-(clock<0.5))*clock;\n",file=AFileName,append=TRUE)
	cat("#declare nrnClock=clock;\n",file=AFileName,append=TRUE)
        
    if (length(Colour)==1) Colour=rep(Colour,length(PreNeurons))
    for (i in 1:length(PreNeurons)){
	cat("//Neuron",PreNeurons[[i]]$NeuronName,"(",PreNeurons[[i]]$CellType,")\n",file=AFileName,append=TRUE)
	WriteBlendedPOVCore(PreNeurons[[i]],PostNeurons[[i]],AFileName,WithBlob,WithCone,Colour[i],InMem,NeuronCentre,...)
    }
    

}


WriteBlendedPOVCore<-function(PreNeuron,PostNeuron,AFileName,WithBlob,WithCone,Colour=NULL,InMem,NeuronCentre,Scale=c(1,-1,-1)){
    # This function represents the guts of WriteBlendedPovFile
    # I modularised it so that it could be used in WriteBlendedPOVFiles
    cat("union {\n",file=AFileName,append=TRUE)
    # FIGURE OUT IF WE NEED TO USE THE BLOB command
    if(WithBlob) {cat("blob { threshold .25\n",file=AFileName,append=TRUE)}

	if(!is.null(Scale)){
		PreNeuron=ScaleNeuronForPov(PreNeuron,Scale)
		PostNeuron=ScaleNeuronForPov(PostNeuron,Scale)
	}
	
	
    # OK GUTS START HERE
    # BASIC IDEA IS TO ITERATE THROUGH
    # Neuron$d starting a cylinder for every line which has non zero parent
    # and getting the position of the parent to finish the cylinder
    # Set the width to be the average of the TwoWidths
    TotalLines<-nrow(PreNeuron$d)
    LineBuffer<-""
    PrePointArray=PreNeuron$d[,c("X","Y","Z","W","Parent")]
    # NOW calculate the difference between post and pre
    DeltaPointArray=PostNeuron$d[,c("X","Y","Z")]
    DeltaPointArray=DeltaPointArray-PrePointArray[1:3]
    names(DeltaPointArray)=c("DX","DY","DZ")
    # and make a big df with everything
    PointArray=cbind(PrePointArray[1:4],DeltaPointArray)[,c("X","DX","Y","DY","Z","DZ","W")]
    
    # Find the points corresponding to the head and end of each segment
    ValidSegEnds=which(PrePointArray$Parent>0)
    ValidSegHeads=PrePointArray$Parent[ValidSegEnds]
    # Make a new array putting it all together
    # First preEnd, PreHead then PostEnd,PostHead
    NewArray=cbind(PointArray[ValidSegEnds,],PointArray[ValidSegHeads,])
    # Col Names to "X1"  "DX1" "Y1"  "DY1" "Z1"  "DZ1" "W1"  "X2"  "DX2" ...
    names(NewArray)=as.vector(t(sapply(c("X","DX","Y","DY","Z","DZ","W"),paste,c("1","2"),sep="")))
    if (is.null(Colour)){
	if(!is.null(PreNeuron$d$SegTypes)){
	    # set colour from seg types it exists
	    SegType<-PreNeuron$d$Label[ValidSegHeads]
	    # Set the segment type to Axon if the two ends are different
	    thisSegColour=SegType
	} else {
	    # else use a default of red
	    SegColours=rep(1,nrow(NewArray))
	}
    } else {
	# ie if a colour is specified, that over-rides the default colour
	SegColours=rep(Colour,nrow(NewArray))
    }
    
    FullInstructions=MakeBlendedPOVConeArray(NewArray,SegColours,WithBlob)
    cat(FullInstructions,file=AFileName,sep="\n",append=TRUE)
    cat(":")
        
    # NOW THE END OF THE NEURON
    cat("normal { bumps 0.4 scale 0.2 }","finish { phong 1 }}",ifelse(WithBlob,"}",""),sep="\n",file=AFileName,append=TRUE)
}
    

WriteBlendedPOVFile<-function(PreNeuron,PostNeuron,AFileName="",WithBlob=F,WithCone=TRUE,FieldOfView=100,Colour=3,InMem=F){
    
    if (!is.list(PreNeuron)| !is.list(PostNeuron)){
	warning("Cannot understand passed neurons")
	return(F)
    }

    if(AFileName=="") AFileName<-paste(strsplit(PreNeuron$NeuronName,'[.]')[[1]][1],sep="",".pov")
    
    if(!file.create(AFileName)) stop(paste("Couldn\'t create file",AFileName))
    
    # OK Lets add a header
    cat(paste("//",AFileName),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=AFileName,append=TRUE)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=AFileName,append=TRUE)

    # Now lets figure out where the centre of the neuron is:
    # OK these calculations work
    # We'll subtract NeuronCentre from the points written out later
    MinAxes<-sapply(rbind(PreNeuron$d[,c("X","Y","Z")],PostNeuron$d[,c("X","Y","Z")]),min)
    MaxAxes<-sapply(rbind(PreNeuron$d[,c("X","Y","Z")],PostNeuron$d[,c("X","Y","Z")]),max)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)

    # Now Position the camera
    AxisExtents<-abs(MaxAxes-MinAxes)
    CameraLoc<-NeuronCentre
    # becuase when I've set it all up, the neuron will be centred around 0,0,0
    CameraLoc[c("X","Y")]<-0
    # In order to see the whole neuron, I think the camera should be as far
    # back as the widest axis extent
    # Changed this 020306 so that we are now looking from Anterior 
    # i.e. inserted the -1
    # Also decided that 80% of the AxisExtents was enough to see all
    # of the neuron
    CameraLoc["Z"]<--1*(NeuronCentre["Z"]-(FieldOfView/100)*max(AxisExtents))
    
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(0,0,0),collapse=",")
    TCamera<-paste("camera { location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
    cat(TCamera,"\n\n",file=AFileName,append=TRUE)

    # AND FINALLY THE LIGHT
    LightPos<-CameraLoc
    # Decided I preferred the light to be at (0,0,SomeZVal)
    ##LightPos[c("X","Y")]<-LightPos["Z"]
    TLightPos<-paste(LightPos,collapse=",")
    TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
    cat(TLight,"\n\n",file=AFileName,append=TRUE)
    
    # DECIDED TO ADD SECOND LIGHT SOURCE TO BRIGHTEN THINGS UP
    cat("light_source { <0,0,0> color White } ",file=AFileName,append=TRUE)

    # One more thing:
    # Check if there are different segment labels for this neuron
    # and if so add colour instructions
    if (!is.null(PreNeuron$SegTypes)){
	for ( i in 1:PreNeuron$NumSegs){
	    segpoints<PreNeuron$SegList[[i]]
	    PreNeuron$d$Label[segpoints]<-PreNeuron$SegType[i]
	}
	
    } else {
	# If not, set the default colour
	PreNeuron$d$Label<-rep(Colour,length(PreNeuron$d$Label))
    }
    
    WriteBlendedPOVCore(PreNeuron,PostNeuron,AFileName,WithBlob,WithCone,Colour,InMem,NeuronCentre)
    
    # START OF THE NEURON ITSELF
}
    

# Function WritePOVFile
WriteMultiPOVFile<-function(Neurons,AFileName="Neurons.pov",WithBlob=F,WithCone=TRUE,FieldOfView=100,Colour=3,InMem=F,...){
	# Instead of being adapted from WritePOVFile this is actually adapted
	# from WriteBlendedPOVFiles and calls WritePOVCore to do the main work
    if (!is.list(Neurons)){
	warning("Cannot understand passed lists of neurons")
	return(F)
    }

    # OK Lets add a header
    cat(paste("//",AFileName),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=AFileName,append=F)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=AFileName,append=TRUE)
    MinAxes=c(0,0,0)
    MaxAxes=c(168,168,88)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)
    # Don't think I want to shift the origin
    NeuronCentre=c(0,0,0)

    # Set the camera
    CameraLoc=c(84,84,-160)
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(84,84,44),collapse=",")
    TCamera<-paste("camera { perspective orthographic location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
#	TCamera<-paste("camera { perspective orthographic location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")

    cat(TCamera,"\n\n",file=AFileName,append=TRUE)
    
    # AND FINALLY THE LIGHT
    LightPos<-CameraLoc
    # Decided I preferred the light to be at (0,0,SomeZVal)
    ##LightPos[c("X","Y")]<-LightPos["Z"]
    TLightPos<-paste(LightPos,collapse=",")
    TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
    cat(TLight,"\n\n",file=AFileName,append=TRUE)
    
    # This can be used as a way of moving the neurons for only part of the time
    cat("#declare nrnClock=clock;\n",file=AFileName,append=TRUE)
    
    # DECIDED TO ADD SECOND LIGHT SOURCE TO BRIGHTEN THINGS UP
    cat("light_source { <0,0,0> color White } \n",file=AFileName,append=TRUE)

    # Make a list of included files
    for(i in 1:length(Neurons)){
		cat("//Neuron",Neurons[[i]]$NeuronName,Neurons[[i]]$CellType,"Colour",Colour[i],")\n",file=AFileName,append=TRUE)
    }

    if (length(Colour)==1) Colour=rep(Colour,length(Neurons))
    for (i in 1:length(Neurons)){
		cat("//Neuron",Neurons[[i]]$NeuronName,"(",Neurons[[i]]$CellType,")\n",file=AFileName,append=TRUE)
		WritePOVCore(Neurons[[i]],AFileName,WithBlob,WithCone,Colour[i],InMem,NeuronCentre,...)
    }
    

}

WritePOVCore<-function(Neuron,AFileName,WithBlob,WithCone,Colour=NULL,InMem,NeuronCentre,AxisDirections=c(1,-1,1),WidthOverride){
	# This function is adapted from WriteBlendedPovCore
	# and is designed to be used with WritePOVFiles
	cat("union {\n",file=AFileName,append=TRUE)
	# FIGURE OUT IF WE NEED TO USE THE BLOB command
	if(WithBlob) {cat("blob { threshold .25\n",file=AFileName,append=TRUE)}
	
	# OK GUTS START HERE
	# BASIC IDEA IS TO MAKE A BIG DF WITH THE FOLLOWING 
	# cone { <56.6187612888,121.3017820678 ,48>, 1.13 <56.44802438116 ,121.3017820678 ,48 >,1.13  pigment{ Blue } }

	# Like RGL, POV uses a standard cartesian geometry system
# 	if(any(AxisDirections!=1)){
# 		Neuron$d[,c("X","Y","Z")]=t(t(Neuron$d[,c("X","Y","Z")])*AxisDirections)
# 	}
	Neuron$d$Y=168.78-Neuron$d$Y
	
	if(!missing(WidthOverride)) Neuron$d$W=WidthOverride
	TotalLines<-nrow(Neuron$d)
	LineBuffer<-""
	PointArray=Neuron$d[,c("X","Y","Z","W","Parent")]
	
	# Find the points corresponding to the head and end of each segment
	ValidSegEnds=which(PointArray$Parent>0)
	ValidSegHeads=PointArray$Parent[ValidSegEnds]
	# Make a new array putting it all together
	NewArray=cbind(PointArray[ValidSegEnds,c("X","Y","Z","W")],PointArray[ValidSegHeads,c("X","Y","Z","W")])
	# Col Names to "X1"  "Y1"  "Z1"  "DZ1" "W1"  "X2"  "DX2"  ...
	names(NewArray)=as.vector(t(sapply(c("X","Y","Z","W"),paste,c("1","2"),sep="")))
	if (is.null(Colour)){
		if(!is.null(Neuron$d$SegTypes)){
			# set colour from seg types it exists
			SegType<-Neuron$d$Label[ValidSegHeads]
			# Set the segment type to Axon if the two ends are different
			thisSegColour=SegType
		} else {
			# else use a default of red
			SegColours=rep(1,nrow(NewArray))
		}
	} else {
		# ie if a colour is specified, that over-rides the default colour
		SegColours=rep(Colour,nrow(NewArray))
	}
	
	FullInstructions=MakePOVConeArray(NewArray,SegColours,WithBlob)
	cat(FullInstructions,file=AFileName,sep="\n",append=TRUE)
	cat(":")
	
	# NOW THE END OF THE NEURON
	cat("normal { bumps 0.4 scale 0.2 }","finish { phong 1 }}",ifelse(WithBlob,"}",""),sep="\n",file=AFileName,append=TRUE)
}

WritePOVFile<-function(ANeuron,AFileName="",WithBlob=F,WithCone=TRUE,FieldOfView=100,Colour=3,InMem=F){
# FieldOfView determines how far back the camera is positioned.  100 => as far
# back as the max extent in any axis
# 020620 Added new option InMem which defaults to F
# idea is to try and do the cat in memory to speed things up a bit
    
    if (is.character(ANeuron)){
	ANeuron<-MyNeurons[[GetNeuronNum(ANeuron)]]
    }
    if (is.numeric(ANeuron)){
	ANeuron<-MyNeurons[[ANeuron]]
    }
    
    if (!is.list(ANeuron)){
	warning("Cannot understand passed neuron")
	return(F)
    }

    if(AFileName=="") AFileName<-paste(strsplit(ANeuron$NeuronName,'[.]')[[1]][1],sep="",".pov")
    
    if(!file.create(AFileName)) stop(paste("Couldn\'t create file",AFileName))
    
    # OK Lets add a header
    cat(paste("//",AFileName),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=AFileName,append=TRUE)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=AFileName,append=TRUE)

    # Now lets figure out where the centre of the neuron is:
    # OK these calculations work
    # We'll subtract NeuronCentre from the points written out later
    MinAxes<-sapply(ANeuron$d[,c("X","Y","Z")],min)
    MaxAxes<-sapply(ANeuron$d[,c("X","Y","Z")],max)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)

    # Now Position the camera
    AxisExtents<-abs(MaxAxes-MinAxes)
    CameraLoc<-NeuronCentre
    # becuase when I've set it all up, the neuron will be centred around 0,0,0
    CameraLoc[c("X","Y")]<-0
    # In order to see the whole neuron, I think the camera should be as far
    # back as the widest axis extent
    # Changed this 020306 so that we are now looking from Anterior 
    # i.e. inserted the -1
    # Also decided that 80% of the AxisExtents was enough to see all
    # of the neuron
    CameraLoc["Z"]<--1*(NeuronCentre["Z"]-(FieldOfView/100)*max(AxisExtents))
    
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(0,0,0),collapse=",")
    TCamera<-paste("camera { location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
    cat(TCamera,"\n\n",file=AFileName,append=TRUE)

    # AND FINALLY THE LIGHT
    LightPos<-CameraLoc
    # Decided I preferred the light to be at (0,0,SomeZVal)
    ##LightPos[c("X","Y")]<-LightPos["Z"]
    TLightPos<-paste(LightPos,collapse=",")
    TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
    cat(TLight,"\n\n",file=AFileName,append=TRUE)
    
    # DECIDED TO ADD SECOND LIGHT SOURCE TO BRIGHTEN THINGS UP
    cat("light_source { <0,0,0> color White } ",file=AFileName,append=TRUE)

    # One more thing:
    # Check if there are different segment labels for this neuron
    # and if so add colour instructions
    if (!is.null(ANeuron$SegTypes)){
	for ( i in 1:ANeuron$NumSegs){
	    segpoints<-ANeuron$SegList[[i]]
	    ANeuron$d$Label[segpoints]<-ANeuron$SegType[i]
	}
	
    } else {
	# If not, set the default colour
        ANeuron$d$Label<-rep(Colour,length(ANeuron$d$Label))
    }
    
    
    # START OF THE NEURON ITSELF
    cat("union {\n",file=AFileName,append=TRUE)
    # FIGURE OUT IF WE NEED TO USE THE BLOB command
    if(WithBlob) cat("blob { threshold .25\n",file=AFileName,append=TRUE)
        
    # OK GUTS START HERE
    # BASIC IDEA IS TO ITERATE THROUGH
    # Neuron$d starting a cylinder for every line which has non zero parent
    # and getting the position of the parent to finish the cylinder
    # Set the width to be the average of the TwoWidths
    TotalLines<-length(ANeuron$d[,1])
    LineBuffer<-""
    for (PresentLine in 1:TotalLines){
	# Check that there is a valid parent
	if( (ParentLine<-ANeuron$d$Parent[PresentLine]) > 0 ){
	    # Ok make a cylinder
	    Point1<-ANeuron$d[PresentLine,c("X","Y","Z")]-NeuronCentre
	    Point2<-ANeuron$d[ParentLine,c("X","Y","Z")]-NeuronCentre
	    Widths<-ANeuron$d$W[c(PresentLine,ParentLine)]
	    Width<-mean(Widths)
	    SegType<-ANeuron$d$Label[PresentLine]
	    # Set the segment type to Axon if the two ends are different
	    if (SegType!=ANeuron$d$Label[ParentLine]) SegType<-1	    
	    # Write out the cylinder
	    if (!WithCone){
		StrBuf<-MakePOVCylinder(Point1,Point2,Width,SegType,WithBlob)
		# Trying out a new in memory version
		if(InMem){
		    FileBuffer<-c(FileBuffer,StrBuf)
		} else {
		    cat(StrBuf,file=AFileName,sep="\n",append=TRUE)		    
		}
	    } else {
		# or write out the cone
		StrBuf<-MakePOVCone(Point1,Point2,Widths,SegType,WithBlob)
		if(InMem){
		    FileBuffer<-c(FileBuffer,StrBuf)
		} else {
		    cat(StrBuf,file=AFileName,sep="\n",append=TRUE)		    
		}
	    }
	}
	cat('.')
	
    } # End of for loop
    # Still need to figure out what is required by way of header / footer 
    # particularly setting up the camera.
    # also rotating etc.
    
    # If we were doing this in memory, now need to dump everything out
    if(InMem){
	cat(LineBuffer,sep="\n",AFileName,sep="\n",append=TRUE)
    }
    
    
    # NOW THE END OF THE NEURON
    cat("normal { bumps 0.4 scale 0.2 }","finish { phong 1 }}",ifelse(WithBlob,"}",""),sep="\n",file=AFileName,append=TRUE)
    
}

# Utility function to make a cylinder joining two points
MakePOVCylinder<-function(Point1,Point2,Width,CylCol=1,WithBlob){

    # Return an empty string if the points are the same,
    # since zero length cones/cylinders upset PovRay
    if(all(Point1==Point2)){
	return("")
    }

    # nb Point1 is an array so to separate the elements of the array,
    # I need to use collapse rather than sep
    TPoint1<-paste(Point1,collapse=",")
    TPoint2<-paste(Point2,collapse=",")    
    
    CylColText<-switch(toString(CylCol), '1'='Blue', '2'='Red','3'='Green')
    
    # Added an extra option on 020306 which adds ,1 to the line
    # before pigment{} if WithBlob is TRUE
    TCylinder<-paste("cylinder { <",TPoint1,"> <",TPoint2,"> ,",Width,ifelse(WithBlob," ,1 "," ")," pigment{ ",CylColText," } }",sep="")
   
    return(TCylinder)
}

# Utility function to make a cone joining two points
MakePOVCone<-function(Point1,Point2,Widths,ConeCol=1,WithBlob){

    # Return an empty string if the points are the same,
    # since zero length cones/cylinders upset PovRay
    if(all(Point1==Point2)){
	return("")
    }
    
    # nb Point1 is an array so to separate the elements of the array,
    # I need to use collapse rather than sep
    TPoint1<-paste(Point1,collapse=",")
    TPoint2<-paste(Point2,collapse=",")    
    
    ConeColText<-switch(toString(ConeCol), '1'='Blue', '2'='Red','3'='Green')
    
    # Added an extra option on 020306 which adds ,1 to the line
    # before pigment{} if WithBlob is TRUE
    TCone<-paste("cone { <",TPoint1,">, ",Widths[1]," <",TPoint2,">,",Widths[2],ifelse(WithBlob," ,1 "," ")," pigment{ ",ConeColText," } }",sep="")

    return(TCone)
}
# Utility function to make a cone joining two points
# This one generates a cone with a variable which can take a value from 
# 0 to 1 thereby transforming it from a cone at
# Point1Pre to Point1Post to a cone at Point1Post,Point2Post
MakeBlendedPOVCone<-function(p1pre,p2pre,p1post,p2post,Widths,ConeCol=1,WithBlob){

    # Return an empty string if the points are the same,
    # since zero length cones/cylinders upset PovRay
    if(all(p1pre==p2pre)){
	return("")
    }

    # Calculate the difference betwen pre and post points
    deltap1=p1post-p1pre
    deltap2=p2post-p2pre

    # nb Point1 is an array so to separate the elements of the array,
    # I need to use collapse rather than sep
    # Have also added a clock variable (defined in the header)
    # which will vary between 0 and 1 for animation
    TPoint1<-paste(paste(p1pre,"+ nrnClock *",deltap1),collapse=",")
    TPoint2<-paste(paste(p2pre,"+ nrnClock *",deltap2),collapse=",")    
    
    ConeColText<-switch(toString(ConeCol), '1'='Blue', '2'='Red','3'='Green','4'='Yellow','5'='Cyan','6'='Magenta','7'='White')
    
    # Added an extra option on 020306 which adds ,1 to the line
    # before pigment{} if WithBlob is TRUE
    TCone<-paste("cone { <",TPoint1,">, ",Widths[1]," <",TPoint2,">,",Widths[2],ifelse(WithBlob," ,1 "," ")," pigment{ ",ConeColText," } }",sep="")

    return(TCone)
}
MakeBlendedPOVConeArray<-function(df,ConeCols,WithBlob){
    # Based on an array where each cone in the eventual file is specified 
    # by a single line in the array in the order
    # P1xyz P2xyz Delta1xyz Delta2xyz Width
    # 
    
    # Need to remove any cones of 0 length since this upsets povray
    ZeroLengthCones=(df$X1==df$X2 & df$Y1==df$Y2  & df$Z1==df$Z2)
    df=df[!ZeroLengthCones,]

    # General form of cone:
    # cone { <END1>, RADIUS1, <END2>, RADIUS2 [open] }
    P1=paste("<",df$X1,",",df$Y1,",",df$Z1,"> + nrnClock * <",df$DX1,",",df$DY1,",",df$DZ1,"> ,",df$W1/2)
    P2=paste("<",df$X2,",",df$Y2,",",df$Z2,"> + nrnClock * <",df$DX2,",",df$DY2,",",df$DZ2,"> ,",df$W2/2)
    
	if(!is.character(ConeCols)){
		ConeColText<-sapply( ConeCols,function(x) switch(x, '1'='Blue', '2'='Red','3'='Green','4'='Yellow','5'='Cyan','6'='Magenta','White')  )
	} else {
		ConeColText=ConeCols
	}
	
	
    return(paste("cone {",P1,P2,ifelse(WithBlob," ,1 "," ")," pigment{ ",ConeColText," } }"))
}
MakePOVConeArray<-function(df,ConeCols,WithBlob){
    # Based on an array where each cone in the eventual file is specified 
    # by a single line in the array in the order
    # P1xyz P2xyz Width
    # 
    
    # Need to remove any cones of 0 length since this upsets povray
    ZeroLengthCones=(df$X1==df$X2 & df$Y1==df$Y2  & df$Z1==df$Z2)
    df=df[!ZeroLengthCones,]
    
    # General form of cone:
    # cone { <END1>, RADIUS1, <END2>, RADIUS2 [open] }
    
    P1=paste("<",df$X1,",",df$Y1,",",df$Z1,"> ,",df$W1/2)
    P2=paste("<",df$X2,",",df$Y2,",",df$Z2,"> ,",df$W2/2)

    if(is.character(ConeCols)) ConeColText=ConeCols
    else {
	ConeColText<-sapply( ConeCols,function(x) switch(x, '1'='Blue', '2'='Red','3'='Green','4'='Yellow','5'='Cyan','6'='Magenta','White')  )
    }
    return(paste("cone {",P1,P2,ifelse(WithBlob," ,1 "," ")," pigment{ ",ConeColText," } }"))
}

# Wrapper function to call WritePOVFile
# if given a set of neurons (names or numbers)
# NB renamed to make function clear
WriteToSeparatePOVFiles<-function(NeuronRef,...){    
    cat("Writing",length(NeuronRef),"pov files: ")
    for(i in 1:length(NeuronRef)){
	t<-WritePOVFile(NeuronRef[i],...)
	cat(".")	
    }
    cat(" Finished!\n")
}

# This function will take the name of .SWC file as input and write it 
# directly to a .pov file 
# It won't handle anything outside the currecnt directory at the moment
SWCToPov<-function(FileName,...){
    # Get the file name without the .SWC suffix
    ShortName<-unlist(strsplit(FileName,"[.]"))[1]
    # Read in the Data
    swcarray<-ReadSWCFile(FileName)
    # Parse the data and make a neuronobject out of it
    TmpNeuron<-TmpNeuron<-SWC2Neuron(swcarray,ShortName)
    
    WritePOVFile(TmpNeuron,...)
}

    
MakePOVPolygon<-function(df){
    # expects an n x 3 df or array of points
    NumPoints=nrow(df)
    PolyPoints=apply(df,1,paste,collapse=",")
    PolyPoints2=paste("<",PolyPoints,">",collapse=",\n")
    return(paste("polygon {",paste(NumPoints,","),PolyPoints2,"}",sep="\n"))
}


MakePOVPoints<-function(df){
    # expects an n x 3 df or array of points
    PolyPoints=apply(df,1,paste,collapse=",")
    PolyPoints2=paste("sphere { <",PolyPoints,"> 0.5}",collapse="\n")
    #return(paste("polygon {",paste(NumPoints,","),PolyPoints2,"}",sep="\n"))
}

ParseRegFile<-function(AFileName){
    # Parse a Registration file generated by Torsten's
    # warp registration
    
    # See if we've been given a gzip file
    if (any(grep(".gz$", AFileName))){
        myFile=gzfile(AFileName)
    } else {
        myFile=file(AFileName)
    }
    
    # Read in all the lines of the file
    
    t=readLines(myFile,-1)
    
    
    LastHeaderLine=grep("coefficients",t)
    if(LastHeaderLine<2 | LastHeaderLine>1000) return(NULL)
    Header=t[1:LastHeaderLine]
    ImageKeyLine=grep("model_study",Header)
    ImageKey=gsub(".*([A-Z]{2}[0-9]{1,2}[LRTB]{1}[1-3]{1}).*","\\1",Header[ImageKeyLine])
    
    DimLine=grep("dims",Header)
    dims=strsplit(Header[DimLine]," ")[2:4]
    domainLine=grep("domain",Header)
    domain=strsplit(Header[domainLine]," ")[2:4]
    originLine=grep("origin",Header)
    origin=strsplit(Header[originLine]," ")[2:4]
    
    splinewarpLine=grep("spline_warp",Header)
    # xformlines=Header[(splinewarpLine+2):(splinewarpLine+6)]
    
    if (any(grep(".gz$", AFileName))){
	 myFile=gzfile(AFileName)
     } else {
	 myFile=open(AFileName)
     }
   
    d=read.table(myFile,skip=splinewarpLine+1,col.names=c("Var","X","Y","Z"),nrows=5)
    d$Var=as.character(d$Var)
    df=as.data.frame(t(d[,-1]))
    names(df)=d$Var

    activeLine=grep("active",t)
    CloseBraces=grep("[}]+",t)
    # Line before the first close brace after the start of the active
    LastActiveLine=CloseBraces[CloseBraces>activeLine][1]-1
    ActivePointLines=t[activeLine:LastActiveLine]
    ActivePointLines[1]=gsub(".*active(.*)","\\1",ActivePointLines[1])
    ActivePointLines2=gsub("(000|111)","\\1 ",ActivePointLines)
    ActivePointLines2=gsub("^[ ]+(.*)","\\1",ActivePointLines2)
    ActivePoints=sapply(unlist(strsplit(ActivePointLines2," ")),as.integer)
    
   
    if (any(grep(".gz$", AFileName))){
	 myFile=gzfile(AFileName)
     } else {
	 myFile=open(AFileName)
     }

    x=scan(myFile,skip=LastHeaderLine-1,nlines=activeLine-LastHeaderLine,na.string="coefficients")[-1]
    x=matrix(x,ncol=3,byrow=TRUE)
    close(myFile)

    rval=cbind(x,ifelse(ActivePoints>0,1,0))
    colnames(rval)=c("X","Y","Z","Active")
    row.names(rval)=NULL
    attr(rval,"ImageKey")=ImageKey
    attr(rval,"Affine")=df
    attr(rval,"dims")=dims
    
    return(rval)
}


RotWriteReg<-function(InFile){
    p=ParseRegFile(InFile)
    p$Active[p$Active==0]=-7
    p$Active[p$Active>0]=-2
    OutFile=basename(dirname(InFile))
    OutFile=gsub("^(.*)_warp.*$","\\1_reg.rot",OutFile)
    if(!file.create(OutFile)) stop(paste("Couldn\'t create file",Outfile))
    # Write out some header information
    cat("#",basename(OutFile),"\n",file=OutFile,append=TRUE)
    cat("# created on",date(),"\n",file=OutFile,append=TRUE)    
    # slightly complicated but allows file to 
    # be written out with a header.
    #apply(t(p),2,cat,"\n",file=OutFile,append=TRUE)
    write.table(p,file=OutFile,col.names=F,row.names=F)
}

cubicMesh<-function(nx,ny,nz,w=168.45,h=168.45,d=87){
    # create a cubic mesh with the given
    # width, height, depth and 
    # spacings
    # nb nx equals number of point planes
    
    dx=w/(nx-1)
    dy=h/(ny-1)
    dz=d/(nz-1)
    X=seq(from=0,to=w,by=dx)
    Y=seq(from=0,to=h,by=dy)
    Z=seq(from=0,to=d,by=dz)

    x=rep(X,length(Y)*length(Z))
    y=rep(rep(Y,rep(length(X),length(Y))),length(Z))
    z=rep(Z,rep(length(X)*length(Y),length(Z)))

    All=cbind(x,y,z)
}
WritePovMesh<-function(mesh,col,OutFile="mesh.pov"){
    # expects a cubic mesh generated by cubicMesh and a vector of colours
    # GJ: 040406 Changed to orthographic camera
    if(dirname(OutFile)=="."){
	OutFile=file.path(RotDir,OutFile)
    }
    cat("OutFile =",OutFile,"\n")
    if(!file.create(OutFile)) stop(paste("Couldn\'t create file",Outfile))

    # Write out some header information
    cat(paste("//",OutFile),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=OutFile,append=F)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=OutFile,append=TRUE)
    MinAxes=c(0,0,0)
    MaxAxes=c(168,168,88)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)
    # Don't think I want to shift the origin
    NeuronCentre=c(0,0,0)

    # Set the camera
    CameraLoc=c(200,200,-200)
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(84,84,44),collapse=",")
    TCamera<-paste("camera { orthographic location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
    cat(TCamera,"\n\n",file=OutFile,append=TRUE)
    
    # AND FINALLY THE LIGHT
    LightPos<-CameraLoc
    # Decided I preferred the light to be at (0,0,SomeZVal)
    ##LightPos[c("X","Y")]<-LightPos["Z"]
    TLightPos<-paste(LightPos,collapse=",")
    TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
    cat(TLight,"\n\n",file=OutFile,append=TRUE)
    
    # This can be used as a way of moving the neurons for only part of the time
    #cat("#declare nrnClock=(1-(clock<0.5))*clock;\n",file=OutFile,append=TRUE)
    
    
    # DECIDED TO ADD SECOND LIGHT SOURCE TO BRIGHTEN THINGS UP
    cat("light_source { <0,0,0> color White } \n",file=OutFile,append=TRUE)
    cat("#declare myRadius=1;\n",file=OutFile,append=TRUE)

    textLines=paste("sphere{ <",mesh[,1],",",mesh[,2],",",mesh[,3],"> myRadius pigment{",col,"}}")
    textLines<<-textLines
    cat(textLines,file=OutFile,append=TRUE,sep="\n")
}


WriteBlendedPovMesh<-function(mesh,col,RegMesh,OutFile="BlendedMesh.pov"){
    # expects a cubic mesh generated by cubicMesh,
    # a vector of colours
    # a mesh from one of Torsten's registration files
    # (including Colours)
    # GJ: 040406 Changed to orthographic camera
    if(dirname(OutFile)=="."){
	OutFile=file.path(RotDir,OutFile)
    }
    cat("OutFile =",OutFile,"\n")
    if(!file.create(OutFile)) stop(paste("Couldn\'t create file",Outfile))

    # Write out some header information
    cat(paste("//",OutFile),paste("//",date()),"#include \"colors.inc\"",sep="\n",file=OutFile,append=F)
    # Now lets set the background Colour:
    cat("\n","background { color Black }","\n\n",file=OutFile,append=TRUE)
    MinAxes=c(0,0,0)
    MaxAxes=c(168,168,88)
    NeuronCentre<-apply(cbind(MinAxes,MaxAxes),1,mean)
    # Don't think I want to shift the origin
    NeuronCentre=c(0,0,0)

    # Set the camera
    CameraLoc=c(200,200,-200)
    TCameraLoc<-paste(CameraLoc,collapse=",")
    TCameraLookAt<-paste(c(84,84,44),collapse=",")
    TCamera<-paste("camera { orthographic location <",TCameraLoc,"> look_at <",TCameraLookAt,"> }",sep="")
    cat(TCamera,"\n\n",file=OutFile,append=TRUE)
    
    # AND FINALLY THE LIGHT
    LightPos<-CameraLoc
    # Decided I preferred the light to be at (0,0,SomeZVal)
    ##LightPos[c("X","Y")]<-LightPos["Z"]
    TLightPos<-paste(LightPos,collapse=",")
    TLight<-paste("light_source { <",TLightPos,"> color White }",sep="")
    cat(TLight,"\n\n",file=OutFile,append=TRUE)
    
    # This can be used as a way of moving the neurons for only part of the time
    #cat("#declare nrnClock=(1-(clock<0.5))*clock;\n",file=OutFile,append=TRUE)
    
    
    # DECIDED TO ADD SECOND LIGHT SOURCE TO BRIGHTEN THINGS UP
    cat("light_source { <0,0,0> color White } \n",file=OutFile,append=TRUE)
    cat("#declare myRadius=0.5;\n",file=OutFile,append=TRUE)
    cat("#declare mClock=clock;\n",file=OutFile,append=TRUE)

    dP=RegMesh[,1:3]-mesh
    
    
#     textLines=paste("sphere{ <",mesh[,1],"+mClock*",dP[,1],",",
# 	mesh[,2],"+mClock*",dP[,2],",",
# 	mesh[,3],"+mClock*",dP[,3],",",
# 	"> myRadius pigment{",col,"}",
# 	"translate <",dP[,1],",",dP[,2],",",dP[,3],">","}")
     textLines=paste("sphere{ <",mesh[,1],",",mesh[,2],",",mesh[,3],
 	"> myRadius translate mClock*<",dP[,1],",",dP[,2],",",dP[,3],">","pigment{",col,"}}")
    textLines<<-textLines
    cat(textLines,file=OutFile,append=TRUE,sep="\n")
}

RotWriteMesh<-function(mesh,col,OutFile="mesh.rot"){
    # expects a cubic mesh generated by cubicMesh and a vector of colours
    if(dirname(OutFile)=="."){
	OutFile=file.path(RotDir,OutFile)
    }
    write.table(cbind(mesh,col*-1),file=OutFile,row.names=F,col.names=F)
}

ScaleNeuronForPov<-function(ANeuron,Scale=c(1,-1,-1)){
	d =data.matrix(ANeuron$d[,c("X","Y","Z")])
	dd=t(t(d)*Scale)
	ANeuron$d[,c("X","Y","Z")]=dd
	ANeuron
}
