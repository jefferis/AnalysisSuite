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