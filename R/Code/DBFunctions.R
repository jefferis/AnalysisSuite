# DBFunctions.R
# Utility functions for database queries

TidyDBResult<-function(df,maxFactorLength=60,minUniqueFracToKeepChar=1.0,keepAsChar){
	# Function to tidy up data frame returned by a db query
	charCols=colnames(df)[sapply(df,is.character)]
	if(!missing(keepAsChar)) charCols=setdiff(charCols,keepAsChar)
	for (c in charCols){
		if(length(grep("date",c,ignore.case=TRUE))>0){
			# try and convert to a date
			d<-try(as.Date(df[[c]]),silent=TRUE)
			if(!inherits(d,"try-error")){
				df[[c]]<-d
				next
			} else {
				warning(paste("Failed to convert column",c,"to date format"))
			}
		}
		# FACTORS AND T/F
		df[[c]]=sub("\\N","",df[[c]],fixed=TRUE)
		maxstrlen=max(nchar(df[[c]]))
		if(maxstrlen==1){
			vals=toupper(unique(df[[c]]))
			# TRUE/FALSE values
			if(all(vals%in%c("T","F")))
				df[[c]]=(toupper(df[[c]])=="T")
			else
				df[[c]]=factor(df[[c]])
		} else if(length(unique(df[[c]]))/length(df[[c]]) < minUniqueFracToKeepChar){
			# we will have less than minUniqueFracToKeepChar factor levels
			# so go ahead and convert
			# (idea is that if there would be as many levels as rows we probably
			# don't want to make a factor)
			if(maxstrlen<=maxFactorLength)
				df[[c]]=factor(df[[c]])
		}
	}
	df
}