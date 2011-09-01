# Try to provide a simple online help system for functions outside of
# a package.  Idea is to extract information from roxygen comments
# in a manner similar to halp package

#' This function provides help for roxygen commented functions
#' that do not have conventional help (in Rd files)
hlp<-function(...){
	# try built in help first
	x=help(...)
	if(length(x)>0) return(eval(x))
	# if not, see if we can get somewhere with roxygen comment
	# get source code for function
	al=try(pairlist(...),silent=TRUE)
	if(inherits(al,'try-error')){
		warning('Unable to find function')
		return(invisible(NULL))
	}
	b=body(al[[1]])
	if(is.null(b)) {
		warning('No help or source code found for: ',as.character(al[[1]]))
		return(invisible(NULL))
	}
	# we have source code, let's parse it
	srccode=as.character(attr(b,'wholeSrcref'))
	# line that finishes the function definition
	functionenddefline=as.integer(attr(b,'srcref')[[1]])[1]
	# FIXME: still sensitive to comments after function def (but on same line)
	# or to cases where function name is 
	functiondefline=max(grep("function",srccode[1:functionenddefline],fixed=T))

	if(length(functiondefline)==0 || functiondefline==1) {
		warning('No help found in source code found for: ',as.character(al[[1]]))
		return(invisible(NULL))
	}
	lastline=length(srccode)
	precedingline=functiondefline-1
	# reverse because we want to work backwards from function definition line
	commentlines=rev(grep('^#',srccode[1:precedingline]))
	if(!any(commentlines==precedingline)){
		# preceding line is not a comment
		warning('No help found in source code found for: ',as.character(al[[1]]))
		return(invisible(NULL))
	}
	drcs=diff(commentlines)
	commentbreak=which(drcs < -1)[1]
	helptext=srccode[rev(commentlines[seq(commentbreak)])]
	helptext=sub("^#[']{0,1}[ ]{0,1}","",helptext)
	cat(paste(helptext,collapse="\n"))
	invisible(helptext)
}