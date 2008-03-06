# PrettyHclustPlot.R
# modified from A2Rplot.hclust by
# GJ 2006-04-03

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
# Addendum - see A2R R library for additional caopyright information
 
require(A2R)

sensClassClust<-function(distmat,k=4,method='ward',...){
	lhclust=hclust(distmat,method=method)
	lhclusters=cutree.order(lhclust,k=k)
	SensClass=factor(CoutInfo$SensClass[match(rownames(as.matrix(distmat)),CoutInfo$Glomerulus)])
	PrettyHclustPlot(lhclust,fact.sup=SensClass,k=k,...)
	lhclusters
}

sensClassClust2<-function(distmat,k=4,method='ward',...){
	lhclust=hclust(distmat,method=method)
	lhclusters=cutree.order(lhclust,k=k)
	SensClass=factor(CoutInfo$SensClass[match(rownames(as.matrix(distmat)),CoutInfo$Glomerulus)])
	plot(lhclust)
	
	lhclusters
}


PrettyHclustPlot<-
function (x, k = 2, col.up = "black", col.down = rainbow(k), 
	lty.up = 2, lty.down = 1, lwd.up = 1, lwd.down = 2, type = c("rectangle", 
		"triangle"), knot.pos = c("mean", "bary", "left", "right", 
		"random"), criteria, fact.sup, show.labels = TRUE, only.tree = FALSE, 
	main = paste("Colored Dendrogram (", k, " groups)"), boxes = TRUE, 
	members, FactorLabel=deparse(substitute(fact.sup)),...) 
{
	if (missing(members)) 
		members <- NULL
	opar <- par(no.readonly = TRUE)
	knot.pos <- match.arg(knot.pos)
	type <- match.arg(type)
	if (k < 2) 
		stop("k must be at least 2")
	._a2r_counter <<- 0
	._a2r_hclu <<- x
	._a2r_envir <<- environment()
	nn <- length(x$order) - 1
	._a2r_height_cut <<- mean(x$height[nn - k + 1:2])
	._a2r_group <<- 0
	n.indiv <- length(x$order)
	groups.o <- cutree.order(x, k = k)[x$order]
	bottom <- if (is.null(members)) 
		0
	else x$height[nn] * -0.2
	if (only.tree) {
		if (is.null(members)) 
			plot(0, type = "n", xlim = c(0.5, n.indiv + 0.5), 
				ylim = c(bottom, x$height[nn]), xaxs = "i", axes = FALSE, 
				xlab = "", ylab = "")
		else plot(0, type = "n", xlim = c(0.5, sum(members) + 
			0.5), ylim = c(bottom, x$height[nn]), xaxs = "i", 
			axes = FALSE, xlab = "", ylab = "")
		.rec.hclust(nn, col = col.up, lty = lty.up, lwd = lwd.up)
		if (boxes) {
			axis(2)
			box()
		}
		return(NULL)
	}
	matlayout <- matrix(c(2, 4, 6, 1, 3, 5), nc = 2, nr = 3)
	widths <- c(1, 9)
	heights <- c(8, 1, 1)
	if (!show.labels) {
		matlayout <- matrix(c(2, 4, 1, 3), nc = 2, nr = 2)
		widths <- c(1, 9)
		heights <- c(9, 1)
	}
	if (!missing(fact.sup)) {
		heights <- c(8, 1, 1)
	}
	if (missing(criteria) & missing(fact.sup)) {
		matlayout <- matrix(c(2, 4, 1, 3), nc = 2, nr = 2)
		widths <- c(1, 9)
		heights <- c(9, 1)
	}
	layout(matlayout, width = widths, height = heights)
	par(mar = c(0, 0, 3, 4))
	if (is.null(members)) 
		plot(0, type = "n", xlim = c(0.5, n.indiv + 0.5), ylim = c(bottom, 
			x$height[nn]), xaxs = "i", axes = FALSE, xlab = "", 
			ylab = "")
	else plot(0, type = "n", xlim = c(0.5, sum(members) + 0.5), 
		ylim = c(bottom, x$height[nn]), xaxs = "i", axes = FALSE, 
		xlab = "", ylab = "")
	.rec.hclust(nn, col = col.up, lty = lty.up, lwd = lwd.up)
	title(main)
	if (boxes) {
		box()
		axis(4)
	}
	if (!missing(criteria)) {
		par(mar = c(0, 0, 3, 0))
		plot(0, type = "n", xlim = range(criteria), ylim = c(0, 
			x$height[nn]), axes = FALSE, xlab = "", ylab = "")
		par(las = 2)
		n.crit <- length(criteria)
		heights.cut <- (tail(x$height, n.crit) + tail(x$height, 
			n.crit + 1)[-(n.crit + 1)])/2
		heights.cut <- rev(heights.cut)
		points(criteria, heights.cut, pch = 21, bg = "red", type = "o")
		points(criteria[k - 1], heights.cut[k - 1], pch = 21, 
			cex = 2, bg = "blue", xpd = NA)
		if (boxes) {
			axis(3)
			box()
		}
	}
	else {
		par(mar = c(0, 0, 3, 0))
		plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, 
			xlab = "", ylab = "")
	}
	if (show.labels) {
		# mar c(bottom, left, top, right) 
		par(mar = c(0, 0, 0, 4))
#		par(mar = c(0, 0, 0, 4))
		par(srt = 90)
		obs.labels <- substr(x$labels[x$order], 1, 6)
		if (is.null(members)) {
#			plot(0, type = "n", xlim = c(0.5, n.indiv + 0.5), 
				plot(0, type = "n", xlim = c(0.5, n.indiv + .7), 
				ylim = c(0, 1), xaxs = "i", axes = FALSE, xlab = "", 
				ylab = "")
			text(1:n.indiv, 0, obs.labels, pos = 4, col = col.down[groups.o])
		}
		else {
			plot(0, type = "n", xlim = c(0.5, sum(members) + 
				0.5), ylim = c(0, 1), xaxs = "i", axes = FALSE, 
				xlab = "", ylab = "")
			xo <- members[x$order]
			text(cumsum(xo) - xo/2, 0, obs.labels, pos = 4, col = col.down[groups.o])
		}
		par(srt = 0)
		if (boxes) {
			box()
		}
		par(mar = c(0, 0, 0, 0))
		plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
			axes = FALSE, xlab = "", ylab = "")
		text(0.5, 0.5, "Labels")
		if (boxes) {
			box()
		}
	}
	if (!missing(fact.sup)) {
		quali <- as.factor(fact.sup)[x$order]
		quanti <- as.numeric(quali)
		par(mar = c(1, 0, 0, 4))
		# par(mar = c(1, 0, 0, 4))
		n.levels <- length(levels(quali))
		plot(0, type = "n", xlim = c(0.5, n.indiv + 0.5), ylim = c(0, 
			n.levels), xaxs = "i", yaxs = "i", axes = FALSE, 
			xlab = "", ylab = "")
		rect(xleft = (1:n.indiv) - 0.5, xright = (1:n.indiv) + 
			0.5, ybottom = quanti - 1, ytop = quanti, col = col.down[groups.o])
		par(las = 1)
		axis(4, (1:n.levels) - 0.5, levels(quali), tick = FALSE)
		if (boxes) {
			box()
		}
		par(mar = c(1, 0, 0, 0))
		plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
			axes = FALSE, xlab = "", ylab = "")
		text(0.5, 0.5, FactorLabel)
		if (boxes) {
			box()
		}
	}
	par(opar)
}
