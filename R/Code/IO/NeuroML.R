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
neuroml_tpl_txt<-'<?xml version="1.0" encoding="UTF-8"?>
<neuroml xmlns="http://www.neuroml.org/schema/neuroml2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.neuroml.org/schema/neuroml2  http://neuroml.svn.sourceforge.net/viewvc/neuroml/NeuroML2/Schemas/NeuroML2/NeuroML_v2alpha.xsd" id="<%=id%>">

    <cell id="<%=id%>">

        <notes>
        <%=notes%>
        </notes>

        <morphology id="morphology_<%=id%>">

            <%%soma%%>
            <%%segments%%>
            <segmentGroup id="Soma">
                <member segment="<%=soma_id%>"/>
            </segmentGroup>

            <segmentGroup id="Neurites">
                <%%neuritesegments%%>
            </segmentGroup>

            <segmentGroup id="all">
                <include segmentGroup="Soma"/>
                <include segmentGroup="Neurites"/>
            </segmentGroup>

            <segmentGroup id="soma_group">
                <include segmentGroup="Soma"/>
            </segmentGroup>
            
        </morphology>

    </cell>
    
</neuroml>
'

seg_tpl_txt<-'<segment id="%d" name="seg_%d">
    <proximal x="%g" y="%g" z="%g" diameter="%g"/>
    <distal x="%g" y="%g" z="%g" diameter="%g"/>
</segment>
'
	# indent by correct amount
	seg_tpl_txt<-gsub('\n','\n            ',seg_tpl_txt,perl=F)
	require(brew)
	
	if(is(f,'connection')) con<-f
	else con=file(f,open='wt')
	on.exit(close(con))
	
	soma_id=0

	segmentParser<-function(btpl){
		if(btpl=='segments'){
			ids=seq.int(len=nrow(x$points))
			p1=x$points-x$vect/2
			p2=x$points+x$vect/2
			segments<-sprintf(seg_tpl_txt,
				ids,ids,
				p1[,1],p1[,2],p1[,3],neurite_diam,
				p2[,1],p2[,2],p2[,3],neurite_diam)
		} else if(btpl=='soma'){
			segment_id=0
			p1=p2=c(unlist(x$soma),soma_diam)
			segments<-sprintf(seg_tpl_txt,
				segment_id,segment_id,
				p1[1],p1[2],p1[3],soma_diam,
				p2[1],p2[2],p2[3],soma_diam)
		} else if(btpl=='neuritesegments'){
			return(paste('                <member segment="',
						seq.int(len=nrow(x$points)),
						'"/>',
						sep='',collapse='\n'))
		}
		paste(segments,collapse='\n')
	}

	brew(text=neuroml_tpl_txt,output=con,tplParser=segmentParser)
}
