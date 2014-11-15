# Simple support functions for Mark Longair's landmarks format from
# the VIB/landmarks package 

ReadLongairPoints<-function(f){
	x=readLines(f)
	numbers=sub(".*\\[(.*)\\].*","\\1",x)
	rnames=sub("(\\\"){0,1}(.?*):\\s+\\[.*","\\2",x)
	rnames=sub("\\\"$","",rnames)
	tc=textConnection(numbers)
	points=read.csv(tc,col.names=c("X","Y","Z"),header=FALSE,row.names=rnames,blank.lines.skip = TRUE)
	close(tc)
	points
}

ReadNeuronsFromLongairTraces<-function(f,MergePaths=FALSE,...){
	neuronList=nat::read.neuron.fiji(f, simplify=FALSE)
	if(MergePaths) {
		if(length(neuronList)==1) return(neuronList[[1]])
		else MergeUnconnectedPathsToSingleNeuron(neuronList)
	}
	else neuronList
}
