# ClusterFunctions.R

# a set of support functions associated with clustering

#RELEASE
#BEGINCOPYRIGHT
###############
# R Source Code to accompany the manuscript
#
# "Comprehensive Maps of Drosophila Higher Olfactory Centers: 
# Spatially Segregated Fruit and Pheromone Representation"
# Cell (2007), doi:10.1016/j.cell.2007.01.040
# by Gregory S.X.E. Jefferis*, Christopher J. Potter*
# Alexander M. Chan, Elizabeth C. Marin
# Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo
#
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# See flybrain.stanford.edu for further details
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################
#ENDMAINCOPYRIGHT

# source(file.path(CodeDir,"ClusterFunctions.R"))
sensClassClust<-function(distmat,k=4,method='ward',...){
	require(A2R)
	lhclust=hclust(distmat,method=method)
	lhclusters=cutree.order(lhclust,k=k)
	SensClass=CoutInfo$SensClass[match(rownames(as.matrix(distmat)),CoutInfo$Glomerulus)]
	PrettyHclustPlot(lhclust,fact.sup=factor(SensClass),k=k,...)
	lhclusters
}

sensClassClust<-function(distmat,k=4,method='ward',...){
	if(!inherits(distmat,"hclust")){
		lhclust=hclust(distmat,method=method)
	} else {
		lhclust=distmat
	}
	names=lhclust$labels
	require(A2R)
	lhclusters=cutree.order(lhclust,k=k)
	SensClass=CoutInfo$SensClass[match(names,CoutInfo$Glomerulus)]
	PrettyHclustPlot(lhclust,fact.sup=factor(SensClass),k=k,...)
	lhclusters
}

sensClassClustPlot<-function(lhclust.col,lhclust,lhc,k,
	selectedGroups=c("ab","at","pb"),legtext=c("unknown",selectedGroups),col,...){

	# 2006-05-30
	# Function more or less dedicated to plotting figure 6A1
	if(missing(lhc)) lhc=cutree.order(lhclust,k=k)
	n=length(lhc)
	
	# trim off v prefix where present 
	Glomeruli=sub("^v(.*)$","\\1",lhclust$labels)
	SensClass=as.character(CoutInfo$SensClass[match(Glomeruli,CoutInfo$Glomerulus)])
	SensClass[!SensClass%in%selectedGroups]=".unknown"
	SensClass=factor(SensClass)
	# Reorder according to order of clustering which goes bottom to top by 
	# default
	SensClass=rev(SensClass[lhclust$order])

	# Figure out colours
	if(missing(col)) {
		col=c("grey",rainbow(length(selectedGroups)))
	}
	# If there are no NAs, then we don't need the grey
	if(nlevels(SensClass)==length(selectedGroups)) { col=col[-1]; legtext=legtext[-1]}
	
	layout(matrix(1:2,ncol=2),wid=c(1,0.2))

	newparmar<-oldparmar<-par("mar")
	newparmar[4]=newparmar[4]*1.5
	par(mar=newparmar)
	plot(rev(lhclust.col),horiz=T,axes=F)
	# Legend for the sensillum classes
	legend("topleft", bty="n",
		#title="Sensillum Class",
		cex=1.1,inset=c(0,0.0),legtext,fill=col)

	# Get ready for new plot 
	ylims=par("usr")[3:4] # c(x1, x2, y1, y2)
	#c(bottom, left, top, right)
	par(mar=c(oldparmar[1],0,oldparmar[3],oldparmar[4]))
	image(x=0:1,y=1:n,xlim=c(-1,1),ylim=ylims,rbind(1:length(SensClass),NA), 
		col = col[as.integer(SensClass)], axes = F,xlab="",ylab="")

	# Make Thick White Borders around the set of boxes for each cluster
	tlhc=rev(table(lhc))
	ctlhc=cumsum(tlhc)
	rect(-1,ctlhc-tlhc+0.5,1,ctlhc+0.5,lwd=4,bor='white')
	#rect(-.8,ctlhc-tlhc+0.5,.8,ctlhc+0.5,lwd=4,bor=rev(rainbow(max(lhc))))
	# make thin white rectangles around each individual box
	rect(-.5,1:n-.5,0.5,1:n+0.5,bor='white')
	
	# Make Vertical lines in the colour corresponding to each LH cluster
	for(i in 1:5){
		lines(c(-.75,-.75),c( (ctlhc-tlhc+0.5)[i],ctlhc[i]+0.5),col=rev(rainbow(max(lhc)))[i],lwd=2)
	}
	
	rval=as.character(rev(SensClass))
	names(rval)=lhclust$labels[lhclust$order]
	return(rval)
}

"cutree.order"<-
function (hclu, k = NULL, h = NULL) 
{
# Function copied from A2R library
# by Romain Francois
# available at
# http://addictedtor.free.fr
	coupe <- cutree(hclu, k = k, h = h)
	coupe.or <- coupe[hclu$order]
	coupe.out <- rep(NA, length(coupe))
	j <- 1
	k <- coupe.or[1]
	for (i in 1:length(coupe)) {
		if (coupe.or[i] == k) 
			next
		else {
			coupe.out[which(coupe == k)] <- j
			j <- j + 1
			k <- coupe.or[i]
		}
	}
	coupe.out[is.na(coupe.out)] <- j
	names(coupe.out) <- names(coupe)
	coupe.out
}

"dendrapply.gj"<-
function (X, FUN, ...) 
{
	FUN <- match.fun(FUN)
	if (!inherits(X, "dendrogram")) 
		stop("'X' is not a dendrogram")
	Napply <- function(d) {
		r <- FUN(d, ...)
		if (!is.leaf(d)) {
			if (!is.list(r)) 
				r <- as.list(r)
			# NB gj changed
			# call Napply on the new list r rather than old list d
			for (j in seq(length = length(d))) r[[j]] <- Napply(r[[j]])
		}
		r
	}
	Napply(X)
}

.recursiveMergeCol<-function(d){
	
	.mergeColOnce<-function(x){
		edgecol=function(x) attr(x,"edgePar")$col
		if(!is.leaf(x)) {
			if( isTRUE(  edgecol(x[[1]])==edgecol(x[[2]])  ) ){
				attr(x, "edgePar")=list(col=edgecol(x[[1]]))
			}
		}
		x
	}

	d.old=d
	while(!identical((d.new<-dendrapply(d.old,.mergeColOnce)),d.old)){d.old=d.new}
	d.new
}

 
addCols<-function(d,cols){
	if(!inherits(d,"dendrogram")) d=as.dendrogram(d)
	.addCol <- function(n) {
		if(is.leaf(n)) {
			attr(n, "edgePar") <- list(col=cols[n])
			attr(n, "nodePar") <- list(lab.col=cols[n],pch=NA)
		}
		n
	}
	.removeUnbalancedCols<-function(x){
		edgecol=function(x) attr(x,"edgePar")$col
		if(!is.leaf(x)) {
			#cat(";",unlist(x),"\n")
			if( !isTRUE(  edgecol(x[[1]])==edgecol(x[[2]])  ) ){
				#cat(	unlist(attr(x[[1]], "edgePar")))
				#cat(	unlist(attr(x[[2]], "edgePar")))
				attr(x[[1]], "edgePar")=NULL
				attr(x[[2]], "edgePar")=NULL
				#cat("Setting to null")
			}
		}
		x
	}

	d1=dendrapply(d,.addCol)
	d1=dendrapply(d1,.recursiveMergeCol)
	dendrapply.gj(d1,.removeUnbalancedCols)
}
