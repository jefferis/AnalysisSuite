# Simple example of clustering a set of Projection Neurons
# The DotProperties representation could undoubtedly be optimised
# and the WeightedNNBasedLinesetMatching.dotprops could be replaced
# with the more advanced algorithm from my nblast v2

# load traced neurons from Jefferis, Potter et al 2007
if(!exists("MyNeurons"))
	load(file.path(ObjDir,"MyNeurons.rda"))

# make a copy of the attached dataframe with information about the neurons for convenience
df=attr(MyNeurons,"df")

# convert to dot properties representation suitable for nblast
MyNeurons.dps=nlapply(MyNeurons,DotProperties)

# calculate all by all nblast distance matrix
dpns=plugindist.list(MyNeurons.dps,WeightedNNBasedLinesetMatching.dotprops)

# hierarchical clustering
# Ward's method tries to find compact spherical clusters
# I (GJ) have normslly found this effective
hcpns=hclust(dpns,meth="ward")

plot(hcpns)
# Plot labelling by glomerulus, 
# making the text a bit smaller so that is can be read on a laptop
plot(hcpns,labels=df$Glomerulus,cex=0.4)

