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
  .Deprecated("nat::write.neurons")
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
