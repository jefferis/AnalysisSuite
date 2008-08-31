# jet.R

# function to produce a set of colours for n levels 
# according to the jet colour map - my version is 
# slightly different from the matlab version (last 1/8 where it
# goes from red -> brown is removed
# v2 accelerates the cold blue transition

# 2005-10-15
# New version of jet.colors using intrinsics.  Looks good and actually
# always returns the specified number of levels.
# nb colorRampPalette returns a function that takes one argument
# the number of levels of the palette
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))


jet.colors2<-function(n){
		w=round(n/6)
		red=c(rep(0,2*w),seq(0,1,len=2*w),rep(1,2*w))
		green=c(seq(0,1,len=2*w),rep(1,2*w),seq(1,0,len=2*w))
		blue=c(seq(0.5,1,len=w),rep(1,w),seq(1,0,len=2*w),rep(0,2*w))
		
		rgbvals=cbind(red,green,blue)
		
		x=apply(rgbvals,1,function(x) rgb(x[1],x[2],x[3]))
		# nb decided that I really wanted the colours to be unique
		# there is a glitch that sometimes produces repeated colours
		unique(x)
}

		
		