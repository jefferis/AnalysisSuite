# Some miscellaneous functions related to image analysis

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
cumnorm<-function(x,mu=0,sigma=1) 1/2 * ( 1 + erf((x - mu) / (sigma * sqrt(2))) )

FitCumulativeGaussian<-function(x,cumfreq,muest,sigmaest,...){
	# use est of median
	if(missing(muest)) muest=x[which.min(abs(cumfreq-0.5))]
	# see IQR {stats}
	if(missing(sigmaest)) sigmaest=( x[which.min(abs(cumfreq - 0.75))] - x[which.min(abs(cumfreq - 0.25))] ) / 1.349
	coeffs<-nls(cumfreq~cumnorm(x,mu,sigma),start=list (mu=muest,sigma=sigmaest), ...)
	return(coeffs$m$getPars())	
}

FitCumulativeGaussianToHistogram<-function(h,truncate=1.0,plot=T,...){
	xs=h$breaks[-length(h$breaks)]
	# find cumsum
	cumfreq=cumsum(h$counts)/sum(h$counts)
	x=xs;
	if(truncate!=1.0) {
		ids=seq(1,by=1,to=truncate*length(h$counts))
		if(length(ids)>10){
			x=x[ids]
			cumfreq=cumfreq[ids]
		} else warning("Ignoring truncate since too few points available")
	}
	params=FitCumulativeGaussian(x,cumfreq)
	if(plot){
		plot(xs,cumsum(h$counts)/sum(h$counts))
		lines(xs,cumnorm(xs,mu= params["mu"],sigma= params["sigma"]),type='l',col='red')
		if(truncate!=1.0) abline(v=x[length(x)])
	}
	params
}

MakeHistogramFromNrrd<-function(filename,...){
	tmp=tempfile()
	outfile=NrrdHisto(filename,outfile=tmp,...)
	if(outfile!=tmp) return(NULL)
	h=ReadHistogramFromNrrd(tmp)
	unlink(tmp)
	return(h)
}

optimalDownsamplingSigma<-function(downsampleby=2,sourcesigma=0.5,targetsigma=0.5,pixelSize=1)
{
	# Optimal downsampling - according to http://pacific.mpi-cbg.de/wiki/index.php/Downsample
	# A results of NaN means no downsampling required
	sigma.pixels=sqrt((targetsigma * downsampleby)^2 - sourcesigma^2)
	sigma=sigma.pixels * pixelSize
	sigma
}

findCDFCorner<-function(x,grad=1,scaleXTo1=TRUE)
{
	# take set of pixel intensities, find cumulative histogram
	n=1000
	e=ecdf(x)
	xmax=max(x)
	# easiest way to do this is adjust the target gradient
	# if xmax>>1 dy/dx will always be smaller than if we scaled xvals to max 1
	# so we want to divide target grad by xmax
	if(scaleXTo1) grad=grad/xmax
	xs=seq(from=0,to=xmax,len=n)
	grads=diff(e(xs))/diff(xs)
	# find smoothed derivative
	ls=lowess(xs[-n],grads,f=1/50,delta=1/1000*xmax)
	# find x val where derivative is closest to grad
	ls$x[which.min(abs(ls$y-grad))]
}

#' Calculate table of estimated background intensity parameters for histograms 
#'
#' By default only looks at the bottom 10% of image histogram (param truncate)
#' Suggests a default threshold of background.mu+2*background.sigma
#' By default the dataframe will have rownames of the name of the input file
#' @param images Paths to image histograms in nrrd format (or a directory)
#' @param filestems Character vector of filestems or function that can calculate them
#' @param truncate fraction of histogram to use for background calculations
#' @return data.frame with cols including background.mu, background.sigma
#' @export
#' @seealso \code{\link{NrrdHisto},\link{ReadHistogramFromNrrd}, 
#' \link{FitCumulativeGaussianToHistogram}, \link{UpdateOrCalculateBackgroundParams}}
CalculateBackgroundParams<-function(images,filestems=basename,truncate=0.1)
{
	if(length(images)==1 && file.info(images)$isdir)
		images=dir(images,full=TRUE)
	imagesdf=data.frame(HistogramFile=images,stringsAsFactors=F)
	if(is.function(filestems))
		rownames(imagesdf)=filestems(images)
	else
		rownames(imagesdf)=filestems
	imagesdf$filestem=rownames(imagesdf)
	imagesdf$background.mu=NA
	imagesdf$background.sigma=NA
	# imagesdf$min=NA
	imagesdf$max=NA

	for(r in rownames(imagesdf)){
		histogramfile=imagesdf[r,"HistogramFile"]
		if(file.exists(histogramfile)){
			imagesdf[r,"max"]=ReadNrrdHeader(histogramfile)$axismaxs
			params<-try(FitCumulativeGaussianToHistogram(
				ReadHistogramFromNrrd(histogramfile),truncate=truncate,plot=F))
			if(!inherits(params,"try-error")) {
				imagesdf[r,"background.mu"]=params["mu"]
				imagesdf[r,"background.sigma"]=params["sigma"]
				cat("+")				
			} else cat("file",basename(histogramfile),"has error",params,"\n")
		} else cat("missing",basename(histogramfile),"\n")
	}
	imagesdf$threshold=imagesdf$background.mu+imagesdf$background.sigma*2
	imagesdf$md5=md5sum(imagesdf$HistogramFile)
	imagesdf$mtime=file.info(imagesdf$HistogramFile)$mtime

	attr(imagesdf,"changed")=TRUE
	imagesdf
}

#' Update a table of calculated background params with new/modified histograms
#'
#' For details see CalculateBackgroundParams.
#' @param images Histograms in nrrd format (or a single directory)
#' @param imagesdf Existing dataframe to start from
#' @param ... Additional options passed to CalculateBackgroundParams
#' @return data.frame of background params
#' @export
#' @seealso \code{\link{CalculateBackgroundParams}}
UpdateOrCalculateBackgroundParams<-function(images,imagesdf,...)
{
	# with one argument just go ahead and calculate
	if(missing(imagesdf)) return(CalculateBackgroundParams(images))
	
	if(length(images)==1 && file.info(images)$isdir)
		images=dir(images,full=TRUE)
	
	if(is.character(imagesdf$mtime)){
		imagesdf$mtime=as.POSIXct(strptime(imagesdf$mtime, "%d/%m/%Y %H:%M"))
	} 
	
	# a flag to know if we changed anything
	attr(imagesdf,"changed")=FALSE
	# else just try to calculate new
	newImages<-setdiff(images,imagesdf$HistogramFile)
	# and changed images
	currentMD5s<-md5sum(imagesdf$HistogramFile)
	imagesToCalculate<-c(newImages,imagesdf$HistogramFile[imagesdf$md5!=currentMD5s])
	imagesToKeep<-setdiff(imagesdf$HistogramFile,imagesToCalculate)
	
	# nothing changed so just return what we were given
	if(length(imagesToCalculate)<1) return(imagesdf)
	
	# if we got this far, something has changed 
	newdf=CalculateBackgroundParams(imagesToCalculate,...)
	attr(newdf,"changed")=TRUE
	# if we had to (re)calculate everbody then we're done
	if(length(imagesToKeep)<1) return(newdf)
	
	# splice the calculated and recalculated images together
	newdf=rbind(newdf,subset(imagesdf,HistogramFile%in%imagesToKeep))
	# and sort by HistogramFile
	newdf=newdf[order(newdf$HistogramFile),]
	attr(newdf,"changed")=TRUE
	newdf
}
