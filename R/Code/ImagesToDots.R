# Functions to process images of neurons and generate more compact
# dot-based representations
# Based on Nick Masse's matlab code (esp extract_properties.m)

DotProperties<-function(points,k=20){
	npoints=nrow(points)
	if(npoints<k) stop("Too few points to calculate properties")

	alpha=rep(0,npoints)
	vect=matrix(0,ncol=3,nrow=npoints)

	nns=nn2(points,points,k=k)
	for(i in 1:npoints){
		indNN=nns$nn.idx[i,]
		center_mass=colMeans(points[indNN,])
		
		inertia=colSums(
			(points[indNN,c(1:3,1:3,1:3)] - 
				matrix(center_mass,ncol=9,nrow=k,byrow=T)) *
			(points[indNN,c(1,1,1,2,2,2,3,3,3)] - 
				matrix(center_mass[c(1,1,1,2,2,2,3,3,3)],ncol=9,nrow=k,byrow=T))
		)
		
		v1d1=eigen(matrix(inertia,ncol=3,byrow=T),symmetric=TRUE)
		alpha[i]=(v1d1$values[1]-v1d1$values[2])/sum(v1d1$values)
		vect[i,]=v1d1$vectors[,1]
	}
	return(list(alpha=alpha,vect=vect))
}

# stop()
# 
# points=ReadAmiramesh("/GD/projects/Nick/FruCloneClustering/data/SAKW13-1_gj_dimension_reduced.am")