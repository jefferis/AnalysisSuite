# a very simple R interface to Albert Cardona's trakem2 Java/ImageJ software
# This allows programmatic access to the VectorString3D class that implements polylines
# and Editions that implements approximate matching

# See TestNeuriteBlast.R for an example

# source("/Users/jefferis/projects/AnalysisSuite/R/Code/Startup.R")
# load("MyNeurons.rda") 

# x=GetLongestSegment(MyNeurons[[1]])
# xx=.jnew('ini/trakem2/vector/VectorString3D',x$X,x$Y,x$Z,FALSE)
# y=GetLongestSegment(MyNeurons[[2]])
# yy=.jnew('ini/trakem2/vector/VectorString3D',y$X,y$Y,y$Z,FALSE)
# .jcall('ini/trakem2/vector/VectorString3D','D','distance',xx,integer(1),yy,integer(1))

trakem2.init<-function(ImageJ="/Applications/ImageJ"){
	 library(rJava)
	.jinit(classpath=file.path(ImageJ,'ij.jar'),
		params=paste('-Dplugins.dir=',ImageJ))
	tem2=list.files(file.path(ImageJ,"plugins"),
		recurs=T,patt=glob2rx("TrakEM2_.jar"),full=T)
	if(length(tem2)==0) stop("Unable to find TrakEM2_.jar in ImageJ plugins dir")
	if(length(tem2)>1) warning(paste("Found multiple copies of TrakEM2_.jar in ImageJ plugins dir.\nUsing:",tem2[1],"\n"))
	.jaddClassPath(tem2[1])
	testLib=try(.jcall('ini/trakem2/utils/Utils','S','d2s',10.3,as.integer(2)))
	if(inherits(testLib,'try-error')) stop("Unable to load TrakEM2 java library")
}

trakem2.findBestMatch<-function(line1,line2, delta=0.5, skipEnds=FALSE, maxMut=5, minChunk=0.5){
	# line1 and line2 are collections of points defining a polyline
	vs1=trakem2.VectorString3D(line1)
	vs2=trakem2.VectorString3D(line2)	
	
	# resample to same delta
	trakem2.VectorString3D.resample(vs1,delta)
	trakem2.VectorString3D.resample(vs2,delta)
	
	matchResult=.jcall("ini/trakem2/vector/Compare",'[Ljava/lang/Object;',
	  'findBestMatch',vs1,vs2,delta,FALSE,as.integer(5),.jfloat(0.5))	
	# returns a list of 2 objects
	# the first can be used to find a Levenshtein edit distance
	# ((Editions)ob[0]).getDistance()
	.jcall(matchResult[[1]],'D','getDistance')
}

trakem2.VectorString3D<-function(x,y=NULL,z=NULL,closed=FALSE){
	xyz=xyz.coords(x,y,z)
	.jnew('ini/trakem2/vector/VectorString3D',xyz$x,xyz$y,xyz$z,closed)
}

trakem2.VectorString3D.resample<-function(vs,delta){
	.jcall(vs,'V','resample',delta)
}

trakem2.VectorString3D.getPoints<-function(vs){
	X=.jcall(vs,'[D','getPoints',as.integer(0))
	Y=.jcall(vs,'[D','getPoints',as.integer(1))
	Z=.jcall(vs,'[D','getPoints',as.integer(2))
	cbind(X,Y,Z)	
}

trakem2.VectorString3D.getVectors<-function(vs){
	X=.jcall(vs,'[D','getVectors',as.integer(0))
	Y=.jcall(vs,'[D','getVectors',as.integer(1))
	Z=.jcall(vs,'[D','getVectors',as.integer(2))
	cbind(X,Y,Z)	
}

trakem2.NeuronLongestSegmentDistance<-function(ANeuron,BNeuron,...){
	# calculates the Levenshtein edit distance between two segments of a pair of neurons
	trakem2.findBestMatch(GetLongestSegment(ANeuron),GetLongestSegment(BNeuron),...)
}
