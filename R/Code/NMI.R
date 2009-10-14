# NMI.R
# source(file.path(CodeDir,"NMI.R"))
# function(s) to calculate the normalised mutual information (NMI)
# between a set of densities.  Expects a density array
plogp=function(x) x*log(x)

entropies<-function(d1,d2,breaks=seq(min(d1,d2),max(d1,d2),length=100)){
	# first precompute
	#c1=cut(as.vector(da[,,,1]),breaks,labels=F)
	#c2=cut(as.vector(da[,,,2]),breaks,labels=F)

	c1=cut(as.vector(d1),breaks)
	c2=cut(as.vector(d2),breaks)
	t=table(c1,c2)
	t=t/sum(t)
	HAB=-sum(plogp(t),na.rm=T)
	HA=-sum(plogp(rowSums(t,na.rm=T)),na.rm=T)
	HB=-sum(plogp(colSums(t,na.rm=T)),na.rm=T)
		
	return(c(HA,HB,HAB))
}

NMI<-function(...){
	l=entropies(...)
	return(( l[1]+l[2])/l[3])
}
MI<-function(...){
	l=entropies(...)
	return( l[1]+l[2]-l[3])
}
bothMI<-function(...){
	l=entropies(...)
	NMI=( l[1]+l[2])/l[3]
	MI=l[1]+l[2]-l[3]
	return(c(MI,NMI))
}

MIDist=function(...) 2/NMI(...) - 1

MIDistMatrix<-function(d){
	n=names(d)
	# could just do top triangle!
	sapply(n,function(g1) sapply(n,function(g2) MIDist(d[[g1]],d[[g2]])))
}

DensityStackEntropy<-function(d,nbreaks=100){
	# Takes a stack of densities (4D oe 2D) and for every voxel
	# calculates the entropy across the stack of different cell types
	# at that point
	
	# TOFIX - not sure this is really a sensible definition as is.
	
	# first make table for evey sub-stack
# 	levels=apply(d,2,function(x) cut(x,nbreaks))
# 	#dim(levels)=dim(d)
# 	ps=apply(levels,2,table)/nrow(d)
# 	plogptable=plogp(ps)
	levels=apply(d,2,function(x) cut(x,seq(0,max(x),len=nbreaks+1),include=T))
	ps=apply(levels,2,function(x) hist(x,plot=F,breaks=seq(0,max(x),len=nbreaks+1))$counts)/nrow(d)
	plogptable=plogp(ps)
	fullplogps=sapply(seq(ncol(d)),function(i) plogptable[levels[,i],i])
	rowSums(-fullplogps)
}
