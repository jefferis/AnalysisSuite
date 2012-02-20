# PDBFunctions.R

WriteNeuronToPDB<-function(ANeuron,filename=NULL,suffix="pdb",Force=F,MakeDir=T #,width=8,dp=3
){
	# Function to write out a neuron as a pdb file based on some example perl code from Madan Babu

    if(is.null(filename)) filename=paste(sub("(.*)\\.[^.]*$","\\1",ANeuron$InputFileName),sep=".",suffix)
    if(!Force && file.exists(filename) ){
		warning(paste(filename,"already exists; use Force=T to overwrite"))
		return()
    }
    if(!file.exists(dirname(filename))){
		# either bail
		if(!MakeDir){
		    warning(paste(dirname(filename),"does not exist; use MakeDir=T to overwrite"))
		    return()
		} else {
		    # or try to make a directory
		    if(!dir.create(dirname(filename))){
				warning(paste("Unable to create",dirname(filename)))
		    }
		}
    }
    if(!file.create(filename)){
		warning(paste("Unable to write to file",filename))
		return()
    }
    
    OutFile=file(filename,"w")

	cat("HEADER\nTITLE\nCOMPND\nSOURCE\n",file=OutFile)
	xyz=ANeuron$d[,c("X","Y","Z")]
	xyz$atm=c("O","N")
	PDBLines=sprintf("HETATM%5d  %s   MSE A 124    %8.3f%8.3f%8.3f",seq(xyz$X),xyz$atm,xyz$X,xyz$Y,xyz$Z)
	writeLines(PDBLines,con=OutFile)
    close(OutFile)
}
