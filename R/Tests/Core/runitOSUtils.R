# functions to test OSUtils functions

# runTestFile(file.path(TestDir,"Core","runitOSUtils.R"))
require(RUnit)

test.abs2rel<-function(x)
{
	testans=abs2rel("/Volumes/JData/JPeople/Sebastian/images","/Volumes/JData/JPeople/")
	realans="Sebastian/images"
	checkEquals(testans,realans)
	
	testans=abs2rel("/some/other/path","/Volumes/JData/JPeople/")
	realans="/some/other/path"
	checkEquals(testans,realans)
	
	checkException(abs2rel("/some/other/path","/Volumes/JData/JPeople/",Stop=TRUE))
	
	testans=abs2rel("/some/other/path","/Volumes/JData/JPeople/")
	realans="/some/other/path"
	checkEquals(testans,realans)	
}