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
