
#' Subset a neuronlist returning either a new neuronlist or the names of chosen neurons
#'
#' EITHER use its attached dataframe as the basis of subset operation. Then use
#' rownames of the new dataframe to select neuronlist entries and return that
#' sublist
#' OR apply a function to every item in the list that returns TRUE/FALSE
#' to determine inclusion in output list
#'
#' * When ReturnList is F just return the indices into the list
#' * When INDICES are specified, then use a for loop to iterate over only those
#' members of the list. This is equivalent to nl[INDICES] but is much
#' faster for big lists when memory swapping occurs. Note that any indices not 
#' present in nl will be dropped with a warning
#'
#' @param nl a neuronlist
#' @param INDICES optional indices to subset neuronlist (faster for big lists)
#' @param ReturnList whether to return the selected neurons (when T) or just their names
#' @param ... either a function or column names in the attached dataframe
#' @export
#' @examples
#' #Apply a 3d search function to the first 100 neurons in the neuronlist dataset
#' subset(dps[1:100],function(x) {length(subset(x,s3d))>0},ReturnList=F)
#' #The same but using INDICES, which is up to 100x faster when neuronlist is large
#' subset(dps,function(x) {length(subset(x,s3d))>0},INDICES=names(dps)[1:100])
subset.neuronlist<-function(nl, ..., INDICES=NULL, ReturnList=is.null(INDICES)){
	arglist=try(pairlist(...),silent=TRUE)
	if(!inherits(arglist,"try-error") && is.function(arglist[[1]])){
		# we are going to apply a function to every element in neuronlist 
		# and expect a return value
		if(length(arglist)>1) stop("I don't know how to handle optional function args.",
			" Use an anonymous function instead")
		if(is.null(INDICES)){
			snl=sapply(nl,arglist[[1]])
			if(ReturnList) return(nl[snl])
			else return(names(nl)[snl])
		} else {
			if(inherits(INDICES,"character")){
				snl=logical(length(INDICES))
				names(snl)=INDICES
			} else if(inherits(INDICES,"logical")){
				snl=logical(sum(INDICES))
				names(snl)=names(nl)[INDICES]
			} else if(inherits(INDICES,"integer")){
				snl=logical(length(INDICES))
				names(snl)=names(nl)[INDICES]
			}
			# trim this list of indices down in case any are not present
			missing_names=setdiff(names(snl),names(nl))
			if(length(missing_names)>0){
				if(length(missing_names)>10)
					warning("Dropping ",length(missing_names)," indices, which are not present in neuronlist")
				else warning("Dropping indices: ",paste(missing_names,collapse=", "),"\n, which are not present in neuronlist")
				snl=snl[intersect(names(snl),names(nl))]
			}
			if(ReturnList) {
				newlist=list()
				for (n in names(snl)){
					include=arglist[[1]](nl[[n]])
					if(include) newlist[[n]]=nl[[n]]
				}
				return(newlist)
			}
			else{
				for (n in names(snl)){
					snl[n]=arglist[[1]](nl[[n]])
				}
				return(names(which(snl)))
			} 
		}
	} else {
		df=attr(nl,'df')
		sdf=subset(df,...)
	}
	if(ReturnList) nl[rownames(sdf)]
	else return(rownames(sdf))
}


#' Evaluate expression in the context of dataframe attached to a neuronlist.
#' @param data A neuronlist object
#' @param expr The expression to evaluate
#' @param ... Ignored
with.neuronlist<-function(data,expr,...) {
  eval(substitute(expr), attr(data,'df'), enclos = parent.frame())
}

#' Read one or more neurons from file to a neuronlist in memory
#'
#' This function will cope with the same set of file formats offered by
#' read.neuron.
#' If the paths argument specifies a (single) directory then all files in that 
#' directory will be read unless an optional regex pattern is also specified.
#' neuronnames must specify a unique set of names that will be used as the 
#' names of the neurons in the resultant neuronlist. If neuronnames is a
#' a function then this will be applied to the path to each neuron. The default
#' value is the function basename which results in each neuron being named for 
#' the input file from which it was read.
#' The optional dataframe (df) detailing each neuron should have rownames that
#' match the names of each neuron. It would also make sense if the same
#' key was present in a column of the data frame.  If the dataframe contains
#' more rows than neurons, the superfluous rows are dropped with a warning.
#' If the dataframe is missing rows for some neurons an error is generated.
#' If SortOnUpdate is TRUE then updating an existing neuronlist should result
#' in a new neuronlist with ordering identical to reading all neurons from scratch.
#' @param paths Paths to neuron input files (or directory containing neurons)
#' @param pattern If paths is a directory, regex that file names must match.
#' @param neuronnames Character vector or function that specified neuron names
#' @param df Optional data frame containing information about each neuron
#' @param OmitFailures Omit failures (when TRUE) or leave an NA value in the list
#' @param SortOnUpdate Sort the neuronlist when update adds new neurons 
#' @return neuronlist object containing the neurons
#' @export
#' @seealso \code{\link{read.neuron}}
#' @examples
read.neurons<-function(paths, pattern=NULL, neuronnames=basename, nl=NULL,
	df=NULL, OmitFailures=TRUE, SortOnUpdate=FALSE, ...){
	if(!is.character(paths)) stop("Expects a character vector of filenames")
	
	if(length(paths)==1 && file.info(paths)$isdir)
		paths=dir(paths,pattern=pattern,full=TRUE)
	
	if(is.function(neuronnames))
		nn=neuronnames(paths)
	else
		nn=neuronnames
	duplicateNames=nn[duplicated(nn)]
	if(length(duplicateNames)) {
		stop("Neurons cannot have duplicate names: ",
				paste(duplicateNames,collapse=" "))
	}
	all_names=nn
	names(paths)=nn
	# Handle updates of an existing neuronlist
	if(!is.null(nl)) {
		if(!is.neuronlist(nl)) stop("nl must be a neuronlist")
		new_neurons=setdiff(nn,names(nl))
		old_neurons_we_can_see=intersect(names(nl),nn)
		old_paths=paths[old_neurons_we_can_see]
		
		new_md5s=md5sum(old_paths)
		old_md5s=sapply(nl,"[[","InputFileMD5")
		names(old_md5s)=names(nl)
		old_md5s=old_md5s[old_neurons_we_can_see]
		stopifnot(length(old_md5s)==length(new_md5s))
		modified_neurons=old_neurons_we_can_see[new_md5s!=old_md5s]
		# now just select the paths that need to be (re)loaded
		nn=c(modified_neurons,new_neurons)
		# no paths to load => existing list is up to date
		if(!length(nn)) return(nl)
		message("There are ",length(modified_neurons)," modified neurons",
				" and ",length(new_neurons),'new neurons')
		paths=paths[nn]
	} else nl=neuronlist()
	# Look after the attached dataframe
	if(!is.null(df)){
		matching_rows=intersect(nn,rownames(df))
		if(length(matching_rows)){
			missing_rows=setdiff(nn,matching_rows)
			if(length(missing_rows))
				stop("Some neurons are not recorded in dataframe: ",
						paste(missing_rows,collapse=" "))
			missing_neurons=setdiff(matching_rows,nn)
			if(length(missing_neurons))
				warning(length(missing_neurons), 
						" rows in dataframe do not have a matching neuron.")
		} else {
			stop("Dataframe rownames do not match neuron names.")
		}
	}
	# Actually read in the neurons
	for(n in names(paths)){
		f=paths[n]
		x=try(read.neuron(f))
		if(inherits(x,'try-error')){
			if(OmitFailures) x=NULL
			else x=NA
		}
		nl[[n]]=x
	}
	if(SortOnUpdate) {
		names_missing_from_all_names=setdiff(names(nl),all_names)
		if(length(names_missing_from_all_names)){
			warning("Cannot SortOnUpdate when supplied paths do not include all neurons: ",
					paste(names_missing_from_all_names,collapse=' '))
		} else {
			# nb names_we_have will be ordered like all_names
			names_we_have=intersect(all_names,names(nl))
			# resort if required
			if(!isTRUE(all.equal(names_we_have,names(nl))))
				nl=nl[names_we_have]
		}
	}
	# nb only keep dataframe rows for neurons that were successfully read in
	attr(nl,'df')=df[names(nl),]
	nl
}

#' Write neurons from a neuronlist object to individual files
#' 
#' NB using INDICES to subset a large neuron list can be much faster
#' @param nl neuronlist object
#' @param dir directory to write neurons
#' @param subdir String naming field in neuron that specifies a subdirectory
#'   OR expression to evaluate in the context of neuronlist's df attribute
#' @param INDICES names of neurons in neuronlist to write
#' @param ... Additional arguments passed to write.neuron
#' @return 
#' @author jefferis
#' @export
#' @seealso \code{\link{write.neuron}}
#' @examples
#' \dontrun{
#' write.neurons(MyNeurons,'/path/to/some/dir',subdir='CellType')
#' write.neurons(MyNeurons,'/path/to/some/dir',subdir=file.path(PNType,Glomerulus,Sex))
#' }
write.neuronlist<-function(nl,dir,subdir=NULL,INDICES=names(nl),...){
  if(!file.exists(dir)) dir.create(dir)
  # Construct subdirectory structure based on 
  df=attr(nl,'df')
  ee=substitute(subdir)
  subdirs=NULL
  if(is.call(ee) && !is.null(df)){
    df=df[INDICES,]
    subdirs=file.path(dir,eval(ee,df,parent.frame()))
    names(subdirs)=INDICES
  }
  for(nn in INDICES){
    n=nl[[nn]]
    thisdir=dir
    if(is.null(subdirs)){
      propval=n[[subdir]]
      if(!is.null(propval)) thisdir=file.path(dir,propval)
    } else {
      thisdir=subdirs[nn]
    }
    if(!file.exists(thisdir)) dir.create(thisdir,recursive=TRUE)
    write.neuron(n,dir=thisdir,...)
  }
}

#' Arithmetic for neuron coordinates applied to neuronlists
#'
#' If x is one number or 4-vector, multiply xyz and diameter by that
#' If x is a 3-vector, multiply xyz only
#' TODO Figure out how to document arithemtic functions in one go
#' @param x a neuronlist
#' @param y (a numeric vector to multiply coords in neuronlist members)
#' @return modified neuronlist
#' @export
#' @rdname neuronlist.arithmetic
#' @examples
#' mn2<-MyNeurons[1:10]*2
`*.neuronlist` <- function(x,y) {
	# TODO look into S3 generics for this functionality
	nlapply(x,`*`,y)
}

`+.neuronlist` <- function(x,y) nlapply(x,`+`,y)
`-.neuronlist` <- function(x,y) x+(-y)
`/.neuronlist` <- function(x,y) x*(1/y)

