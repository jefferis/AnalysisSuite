# Try to provide a simple online help system for functions outside of
# a package.  Idea is to extract information from roxygen comments
# in a manner similar to halp package

gjsource<-function (file, local = FALSE, echo = verbose, print.eval = echo, 
		verbose = getOption("verbose"), prompt.echo = getOption("prompt"), 
		max.deparse.length = 150, chdir = FALSE, encoding = getOption("encoding"), 
		continue.echo = getOption("continue"), skip.echo = 0, keep.source = getOption("keep.source")) 
{
	eval.with.vis <- function(expr, envir = parent.frame(), enclos = if (is.list(envir) || 
							is.pairlist(envir)) 
						parent.frame()
					else baseenv()) .Internal(eval.with.vis(expr, envir, enclos))
	envir <- if (local) 
				parent.frame()
			else .GlobalEnv
	have_encoding <- !missing(encoding) && encoding != "unknown"
	if (!missing(echo)) {
		if (!is.logical(echo)) 
			stop("'echo' must be logical")
		if (!echo && verbose) {
			warning("'verbose' is TRUE, 'echo' not; ... coercing 'echo <- TRUE'")
			echo <- TRUE
		}
	}
	if (verbose) {
		cat("'envir' chosen:")
		print(envir)
	}
	ofile <- file
	from_file <- FALSE
	srcfile <- NULL
	if (is.character(file)) {
		if (identical(encoding, "unknown")) {
			enc <- utils::localeToCharset()
			encoding <- enc[length(enc)]
		}
		else enc <- encoding
		if (length(enc) > 1L) {
			encoding <- NA
			owarn <- options("warn")
			options(warn = 2)
			for (e in enc) {
				if (is.na(e)) 
					next
				zz <- file(file, encoding = e)
				res <- tryCatch(readLines(zz), error = identity)
				close(zz)
				if (!inherits(res, "error")) {
					encoding <- e
					break
				}
			}
			options(owarn)
		}
		if (is.na(encoding)) 
			stop("unable to find a plausible encoding")
		if (verbose) 
			cat(gettextf("encoding = \"%s\" chosen", encoding), 
					"\n", sep = "")
		if (file == "") 
			file <- stdin()
		else {
			if (isTRUE(keep.source)) 
				srcfile <- srcfile(file, encoding = encoding)
			file <- file(file, "r", encoding = encoding)
			on.exit(close(file))
			from_file <- TRUE
			loc <- utils::localeToCharset()[1L]
			encoding <- if (have_encoding) 
						switch(loc, `UTF-8` = "UTF-8", `ISO8859-1` = "latin1", 
								"unknown")
					else "unknown"
		}
	}
	exprs <- .Internal(parse(file, n = -1, NULL, "?", srcfile, 
					encoding))
	Ne <- length(exprs)
	if (from_file) {
		close(file)
		on.exit()
	}
	if (verbose) 
		cat("--> parsed", Ne, "expressions; now eval(.)ing them:\n")
	if (Ne == 0) 
		return(invisible())
	if (chdir) {
		if (is.character(ofile)) {
			isURL <- length(grep("^(ftp|http|file)://", ofile)) > 
					0L
			if (isURL) 
				warning("'chdir = TRUE' makes no sense for a URL")
			if (!isURL && (path <- dirname(ofile)) != ".") {
				owd <- getwd()
				if (is.null(owd)) 
					warning("cannot 'chdir' as current directory is unknown")
				else on.exit(setwd(owd), add = TRUE)
				setwd(path)
			}
		}
		else {
			warning("'chdir = TRUE' makes no sense for a connection")
		}
	}
	if (echo) {
		sd <- "\""
		nos <- "[^\"]*"
		oddsd <- paste("^", nos, sd, "(", nos, sd, nos, sd, ")*", 
				nos, "$", sep = "")
	}
	srcrefs <- attr(exprs, "srcref")
	for (i in 1L:Ne) {
		if (verbose) 
			cat("\n>>>> eval(expression_nr.", i, ")\n\t\t =================\n")
		ei <- exprs[i]
		if (echo) {
			if (i > length(srcrefs) || is.null(srcref <- srcrefs[[i]])) {
				dep <- substr(paste(deparse(ei, control = c("showAttributes", 
												"useSource")), collapse = "\n"), 12, 1e+06)
				dep <- paste(prompt.echo, gsub("\n", paste("\n", 
										continue.echo, sep = ""), dep), sep = "")
				nd <- nchar(dep, "c") - 1
			}
			else {
				if (i == 1) 
					lastshown <- min(skip.echo, srcref[3L] - 1)
				dep <- getSrcLines(srcfile, lastshown + 1, srcref[3L])
				leading <- srcref[1L] - lastshown
				lastshown <- srcref[3L]
				while (length(dep) && length(grep("^[[:blank:]]*$", 
								dep[1L]))) {
					dep <- dep[-1L]
					leading <- leading - 1L
				}
				functionheader <- if(leading>1) dep [seq(leading-1)] else NULL
				if(!is.null(functionheader)) attr(exprs[i],"functionheader") = functionheader
				cat("functionheader:\n",paste(functionheader,collapse="\n"))
				dep <- paste(rep.int(c(prompt.echo, continue.echo), 
								c(leading, length(dep) - leading)), dep, sep = "", 
						collapse = "\n")
				nd <- nchar(dep, "c")
			}
			if (nd) {
				do.trunc <- nd > max.deparse.length
				dep <- substr(dep, 1L, if (do.trunc) 
									max.deparse.length
								else nd)
				cat("\n", dep, if (do.trunc) 
							paste(if (length(grep(sd, dep)) && length(grep(oddsd, 
															dep))) 
												" ...\" ..."
											else " ....", "[TRUNCATED] "), "\n", sep = "")
			}
		}
		yy <- eval.with.vis(ei, envir) # this is what actually runs the expression
		i.symbol <- mode(ei[[1L]]) == "name"
		if (!i.symbol) {
			curr.fun <- ei[[1L]][[1L]]
			if (verbose) {
				cat("curr.fun:")
				utils::str(curr.fun)
			}
		}
		if (verbose >= 2) {
			cat(".... mode(ei[[1L]])=", mode(ei[[1L]]), "; paste(curr.fun)=")
			utils::str(paste(curr.fun))
		}
		if (print.eval && yy$visible) {
			if (isS4(yy$value)) 
				methods::show(yy$value)
			else print(yy$value)
		}
		if (verbose) 
			cat(" .. after ", sQuote(deparse(ei, control = c("showAttributes", 
											"useSource"))), "\n", sep = "")
	}
#	invisible(yy)
}

gjhalp<-function (fun, use.pager = logical(0), character.only = FALSE, 
		...) 
{
	if (!is.logical(use.pager)) {
		stop("Argument 'use.pager' must be of type logical")
	}
	if (!is.logical(character.only)) {
		stop("Argument 'character.only' must be of type logical")
	}
	if (is.character(fun) || all(character.only)) {
		fun.name <- fun
	}
	else {
		fun.name <- deparse(substitute(fun))
	}
	fun.obj <- get(fun.name, mode = "function", ...)
	fun.source <- attr(fun.obj, "source")
	opt <- get("halp.options", envir = as.environment("package:halp"), 
			mode = "function", inherits = FALSE)
	if (is.null(opt("pattern"))) {
		stop("Invalid comment pattern for the 'halp' package, see '?halp.options'")
	}
	if (is.null(fun.source)) {
		stop("Source code for function '", fun.name, "' is not available\n", 
				sep = "")
	}
	fun.halp.raw <- gsub(opt("pattern"), "", fun.source[grep(opt("pattern"), 
							fun.source)])
	fun.halp.raw <- gsub("^@param (\\w+)"," + \\1:",fun.halp.raw)
	fun.halp.raw <- gsub("^@(\\w+)"," \\1:",fun.halp.raw)
	if (length(fun.halp.raw) > 0) {
		if ((length(use.pager) > 0 && all(use.pager)) || (length(use.pager) == 
					0 && length(fun.halp.raw) >= opt("pager.minLines"))) {
			files <- tempfile()
			if (is.null(files)) {
				stop("Could not acquire temporary filename for pager display")
			}
			file <- files[1]
			cat(paste(opt("pager.pre"), fun.halp.raw, collapse = "\n", 
							sep = ""), "\n", sep = "", file = file)
			file.show(file)
		}
		else {
			cat(paste(opt("screen.pre"), fun.halp.raw, collapse = "\n", 
							sep = ""), "\n", sep = "")
		}
	}
	
	else {
		warning("The are no halp comments for function '", fun.name, 
				"'", call. = FALSE)
	}
	invisible(fun.halp.raw)
}

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