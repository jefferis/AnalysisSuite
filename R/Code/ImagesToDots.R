# Functions to process images of neurons and generate more compact
# dot-based representations
# Based on Nick Masse's matlab code (esp extract_properties.m)

DotProperties<-function(points,k=20){
	npoints=nrow(points)
	if(npoints<k) stop("Too few points to calculate properties")
	if(ncol(points)!=3) stop("points must be a N x 3 matrix")
	
	alpha=rep(0,npoints)
	vect=matrix(0,ncol=3,nrow=npoints)

	nns=nn2(points,points,k=k)
	# transpose points to 3xN because 
	# R arithemtic of matric / vector operates column-wise
	pointst=t(points)
	for(i in 1:npoints){
		indNN=nns$nn.idx[i,]
		
		pt=pointst[,indNN]
		cpt=pt-rowMeans(pt)
		
		inertia=matrix(0,ncol=3,nrow=3)
		diag(inertia)=rowSums(cpt^2)
		inertia[1,2]<-inertia[2,1]<-sum(cpt[1,]*cpt[2,])
		inertia[1,3]<-inertia[3,1]<-sum(cpt[1,]*cpt[3,])
		inertia[2,3]<-inertia[3,2]<-sum(cpt[2,]*cpt[3,])
		
		v1d1=eigen(inertia,symmetric=TRUE)
		alpha[i]=(v1d1$values[1]-v1d1$values[2])/sum(v1d1$values)
		vect[i,]=v1d1$vectors[,1]
	}
	return(list(alpha=alpha,vect=vect))
}

# stop()
# 
# points=ReadAmiramesh("/GD/projects/Nick/FruCloneClustering/data/SAKW13-1_gj_dimension_reduced.am")