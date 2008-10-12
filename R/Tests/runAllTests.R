# runAllTests.R
# run tests for all the directories of scripts/sources files

# source(file.path(RootDir,"R","Tests","runAllTests.R"))

require(RUnit)
TestDirs=file.path(RootDir,"R","Tests")
cat("Run all tests:\n")
MyTestSuites<-lapply(TestDirs,function(x) defineTestSuite(basename(x),x,testFileRegexp="^runit.+.[srR]{1}$") )	
TestResults=lapply(MyTestSuites,runTestSuite)
print(TestResults)
rm(TestDirs,MyTestSuites)