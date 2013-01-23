# code to provide simple IO to neuro ML 
# see e.g. https://raw.github.com/openworm/CElegansNeuroML/master/CElegans/generatedNeuroML2/AIAL.nml

dotprops2nml<-function(x,f,id,neurite_diam=1,soma_diam=4,notes=NULL,...){
	require(brew)
	neuroml_tpl=file.path(CodeDir,'IO','NeuroML.brew')
	seg_tpl=file.path(dirname(neuroml_tpl),'NeuroML_segment.brew')

	if(is(f,'connection')) con<-f
	else con=file(f,open='wt')
	on.exit(close(con))
	
	soma_id=0

	segmentParser<-function(btpl){
		# just ignore incoming template for time being
		tc=textConnection("segments",open='w')
		on.exit(close(tc))
		if(grepl('segments',btpl)){
			for(i in seq(nrow(x$points))){
				segment_id=i
				p1=c(x$points[i,]-x$vect[i,]/2,neurite_diam)
				p2=c(x$points[i,]+x$vect[i,]/2,neurite_diam)
				brew(seg_tpl,tc)
			}
		} else if(grepl('soma',btpl)){
			segment_id=0
			p1=p2=c(unlist(x$soma),soma_diam)
			brew(seg_tpl,tc)
		}
		
		paste(segments,collapse='\n')
	}

	brew(neuroml_tpl,con,tplParser=segmentParser)
}