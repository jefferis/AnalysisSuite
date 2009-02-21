# Tests for reading Mark Longair's tracing format and points data

test.ReadLongairNeuron.SinglePath<-function(){
	# Test reading of Mark Longair's tracing format
	fieldsToCheck=c( "NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
	"NumSegs", "SegList", "d")
	
	result<-structure(list(NumPoints = 11L, StartPoint = 1, 
	    BranchPoints = numeric(0), EndPoints = c(1, 11), NumSegs = 1L, 
	    SegList = list(c(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)), d = structure(list(
	        PointNo = 1:11, Label = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 
	        2, 2), X = c(228.337419509888, 228.337419509888, 228.337419509888, 
	        227.788531482220, 227.788531482220, 227.239643454552, 
	        227.239643454552, 227.239643454552, 226.690755426884, 
	        226.690755426884, 226.690755426884), Y = c(92.2131886482239, 
	        92.2131886482239, 92.2131886482239, 92.2131886482239, 
	        92.2131886482239, 91.6643006205559, 91.1154125928879, 
	        91.1154125928879, 90.5665245652199, 90.0176365375519, 
	        89.4687485098839), Z = c(39, 40, 41, 42, 43, 44, 45, 
	        46, 47, 48, 49), radius = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 
	        1, 1), Parent = c(-1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)), .Names = c("PointNo", 
	    "Label", "X", "Y", "Z", "radius", "Parent"), row.names = c(NA, 
	    -11L), class = "data.frame")), .Names = c("NumPoints", 
	"StartPoint", "BranchPoints", "EndPoints", "NumSegs", "SegList", 
	"d"))
	
	result.new=ReadNeuronsFromLongairTraces(file.path(TestDir,"SinglePath.traces"))[[1]]
	checkEquals(result,result.new[fieldsToCheck],tol=1e-6)	
}

test.ReadLongairNeuron.PathJoining<-function(){
	# Test reading of Mark Longair's tracing format
	fieldsToCheck=c( "NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
	"NumSegs", "SegList", "d")
	
	result<-structure(list( NumPoints = 11L, StartPoint = 1, 
	    BranchPoints = numeric(0), EndPoints = c(1, 11), NumSegs = 1L, 
	    SegList = list(c(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)), d = structure(list(
	        PointNo = 1:11, Label = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 
	        2, 2), X = c(228.337419509888, 228.337419509888, 228.337419509888, 
	        227.788531482220, 227.788531482220, 227.239643454552, 
	        227.239643454552, 227.239643454552, 226.690755426884, 
	        226.690755426884, 226.690755426884), Y = c(92.2131886482239, 
	        92.2131886482239, 92.2131886482239, 92.2131886482239, 
	        92.2131886482239, 91.6643006205559, 91.1154125928879, 
	        91.1154125928879, 90.5665245652199, 90.0176365375519, 
	        89.4687485098839), Z = c(39, 40, 41, 42, 43, 44, 45, 
	        46, 47, 48, 49), radius = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 
	        1, 1), Parent = c(-1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)), .Names = c("PointNo", 
	    "Label", "X", "Y", "Z", "radius", "Parent"), row.names = c(NA, 
	    -11L), class = "data.frame")), .Names = c( "NumPoints", 
	"StartPoint", "BranchPoints", "EndPoints", "NumSegs", "SegList", 
	"d"))
	
	result.new=ReadNeuronsFromLongairTraces(file.path(TestDir,"SinglePath.traces"))[[1]]
	checkEquals(result,result.new[fieldsToCheck],tol=1e-6)	
}