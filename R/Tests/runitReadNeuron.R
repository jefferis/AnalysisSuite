# TestReadNeuron.R
# Just a little script to profile the contour reading functions
# 
# source(file.path(CodeDir,"TestReadNeuron.R"))
Rprof("TestReadNeuron.out")
f="/GD/AnalysingTraceFiles/TraceFilesWithMBData/VM2/traced VM2/SK17R2_2.clm2.asc"
AscData=ReadAscData(f)
for(i in 1:20){
	if(F){
		RawLHData=GetContourData(AscData,"LH")
		#OK Lets scale/reorient these data based on the info
		#in OrientInfo which was filled in when the axon data was read in
		RawLHData[,c("X","Y","Z")]=InvertAndReorderData(RawLHData[,c("X","Y","Z")],MyNeurons.new[[125]]$OrientInfo)
		RawMBData<-GetContourData(AscData,"MB")
		RawMBData[,c("X","Y","Z")]<-InvertAndReorderData(RawMBData[,c("X","Y","Z")],MyNeurons.new[[125]]$OrientInfo)
	} else {
		ReadNeuronFromAsc(f)
	}
}

Rprof(NULL)
#tFinish=Sys.time()-tStart
#cat("time =",tFinish,"\n")