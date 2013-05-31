# functions to test OSUtils functions

# runTestFile(file.path(TestDir,"Core","runitOSUtils.R"))
require(RUnit)

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
