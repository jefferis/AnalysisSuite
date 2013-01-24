# code to provide simple IO to neuro ML 
# see e.g. https://raw.github.com/openworm/CElegansNeuroML/master/CElegans/generatedNeuroML2/AIAL.nml

#' Convert a dotprops neuron representation into NeuroML format
#'
#' @ details Note that this representation currently produces neurons without any 
#' connections betweeng segments that will consequently be regarded as invalid
#' by most software, including neuroConstruct.
#' @references see http://www.neuroml.org/ for details
#' @param x Dotprops object
#' @param f Path to neuroML ouptut file
#' @param id NeuroML string id for cell
#' @param neurite_diam,soma_diam Default diams for neurites and soma, respectively
#' @param notes Comment block for NeuroML
#' @param ... Additional arguments that could be used to specify neuroml 
#'  blocks FIXME not yet implemented
#' @return Last argument evaluated by brew function
#' @export
#' @seealso \code{\link{brew},\link{dotprops}}
dotprops2nml<-function(x,f,id,neurite_diam=1,soma_diam=4,notes=NULL,...){
	require(brew)
	neuroml_tpl=file.path(CodeDir,'IO','NeuroML.brew')
	seg_tpl=file.path(dirname(neuroml_tpl),'NeuroML_segment.brew')
	# save time by reading this to memory once
	seg_tpl_txt=paste(readLines(seg_tpl),collapse='\n')

	if(is(f,'connection')) con<-f
	else con=file(f,open='wt')
	on.exit(close(con))
	
	soma_id=0

	segmentParser<-function(btpl){
		tc=textConnection("segments",open='w')
		on.exit(close(tc))
		if(grepl('segments',btpl,fixed=TRUE)){
			for(i in seq(nrow(x$points))){
				segment_id=i
				p1=c(x$points[i,]-x$vect[i,]/2,neurite_diam)
				p2=c(x$points[i,]+x$vect[i,]/2,neurite_diam)
				brew(text=seg_tpl_txt,output=tc)
			}
		} else if(grepl('soma',btpl,fixed=TRUE)){
			segment_id=0
			p1=p2=c(unlist(x$soma),soma_diam)
			brew(text=seg_tpl_txt,output=tc)
		}
		
		paste(segments,collapse='\n')
	}

	brew(neuroml_tpl,con,tplParser=segmentParser)
}