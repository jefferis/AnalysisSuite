# functions to test OSUtils functions

# runTestFile(file.path(TestDir,"Core","runitOSUtils.R"))
require(RUnit)

test.abs2rel<-function(x)
{
	testans=abs2rel("/Volumes/JData/JPeople/Sebastian/images","/Volumes/JData/JPeople/")
	realans="Sebastian/images"
	checkEquals(testans,realans)

	testans=abs2rel("/Volumes/JData/JPeople/Sebastian/images","/Volumes/JData/JPeople")
	realans="Sebastian/images"
	checkEquals(testans,realans)
		
	checkException(abs2rel("/some/other/path","/Volumes/JData/JPeople/",
		StopIfNoCommonPath=TRUE),silent=TRUE)
	
	testans=abs2rel("/some/other/path","/Volumes/JData/JPeople/")
	realans="/some/other/path"
	checkEquals(testans,realans)	
}

test.fix.dir<-function(){
	testans = fix.dir("/some/path")
	realans = "/some/path/"
	checkEquals(testans,realans)

	testans = fix.dir("/some/path/")
	realans = "/some/path/"
	checkEquals(testans,realans)

	testans = fix.dir(c("/some/path/","/some/other/path"))
	realans = c("/some/path/","/some/other/path/")
	checkEquals(testans,realans)	
}

test.RunCmdForNewerInput<-function(){
	# make a test directory
	tf=replicate(5,tempfile())
	for (i in 1:2) cat("Hello",i,"!",file=tf[i])
	Sys.sleep(1.5)
	for (i in 3:4) cat("Hello",i,"!",file=tf[i])
	on.exit(unlink(tf[1:4]))

	# one older input file
	checkEquals(RunCmdForNewerInput(NULL,infiles=tf[1],outfile=tf[3]),FALSE)

	# one newer input file
	checkTrue(RunCmdForNewerInput(NULL,infiles=tf[4],outfile=tf[1]))

	# multiple older inputfiles
	checkEquals(FALSE,
		RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[2]),outfile=tf[3]))

	# one newer and one older input file
	checkTrue(RunCmdForNewerInput(NULL,infiles=c(tf[4],tf[1]),outfile=tf[2]))

	# single missing input file
	checkEquals(RunCmdForNewerInput(NULL,infiles=tf[5],outfile=tf[2]),FALSE)

	# empty input file vector
	checkEquals(RunCmdForNewerInput(NULL,infiles=character(0),
		outfile=tf[2]),FALSE)

	# one input missing, another present
	checkEquals(FALSE,
		RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[5]),outfile=tf[2]))

	# missing output file 
	checkTrue(RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[4]),outfile=tf[5]))

	# single older input, multiple newer outputs
	checkEquals(FALSE,
		RunCmdForNewerInput(NULL,infiles=tf[1],outfile=c(tf[3],tf[4])))

	# multiple older inputs, multiple newer outputs
	checkEquals(FALSE,
		RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[2]),outfile=c(tf[3],tf[4])))

	# single input, multiple outputs, of which one is older
	checkTrue(RunCmdForNewerInput(NULL,infiles=tf[3],outfile=c(tf[1],tf[4])))

	# one missing output file, older inputs
	checkTrue(RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[2]),outfile=c(tf[4],tf[5])))

	# multiple inputs, multiple outputs, one older
	checkTrue(RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[4]),outfile=c(tf[2],tf[3])))

	# multiple inputs, multiple outputs, one older
	checkTrue(RunCmdForNewerInput(NULL,infiles=c(tf[1],tf[4]),outfile=c(tf[2],tf[3]),Verbose=TRUE))
}
