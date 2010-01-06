# runAllTests.R
# run tests for all the directories of scripts/sources files

# source(file.path(RootDir,"R","Tests","runAllTests.R"))

require(RUnit)
TestDirs=rownames(subset(file.info(dir(TestDir,full=TRUE)),isdir==TRUE))
cat("Run all tests:\n")
MyTestSuites<-lapply(TestDirs,function(x) defineTestSuite(basename(x),x,testFileRegexp="^runit.+.[srR]{1}$") )	
TestResults=runTestSuite(MyTestSuites)
print(TestResults)
printHTMLProtocol(TestResults,fileName=file.path(TestDir,"LatestTestResults.html"))
rm(TestDirs,MyTestSuites)