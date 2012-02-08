# MarchetteQhullInterface.R
###########################
# Code provided by David Marchette to call qhull
# and then read in results
# I have added qhull.va(matrx) which simply returns the
# surface area and volume of the hull
# Greg Jefferis 011106
# license as below!

# PS I compiled qhull 3.1 and that seems to work fine

# This software developed by David Marchette
# Naval Surface Warfare Center, Dahlgren
# Division, Code B10.  It may be freely used and 
# modified. All warranties, express or implied, are
# disclaimed. 

#RELEASENOCOPY

# 031117 - prevented STDERR from being reported

# get approx volume and surface area of hull
qhull.va<-function(matrx,silent=T,joggle=F){
    if (missing(matrx))
    {
	cat("Matrix argument required\n")
	return(NULL)
    }
		  
    ### data matrix ###
    x <- eval(matrx)
    d <- dim(x)
    n <- d[1]
    dimension <- d[2]
    if (sum(abs(x[!is.na(x)]) == Inf) > 0)
    {
      cat("Sorry, qhull can't handle Infs\n")
      return()
    }
    if (sum(is.na(x)) > 0)
    {
	cat("Sorry, qhull can't handle NAs\n")
	return()
    }
    dfile <- tempfile("unix")
    ofile <- tempfile("unix")
    y <- c(dimension,n)
    write(y, file = dfile)
    #write(t(x),file=dfile,append=T,ncol=dimension)
    write.table(x,file=dfile,append=T,col.names=F,row.names=F)
    
    command <- paste("cat",dfile,"| qhull",sep=" ")
    
    if(joggle) command <- paste(command,"QJ",sep=" ")
    command <- paste(command,"FS",sep=" ")
    if(silent) {
	command <- paste(command,"Pp",sep=" ")
    }
    
    command <- paste(command, "TO",ofile,sep=" ")

    stdout<-system(command, intern=T,ignore.stderr=T)

    q <- scan(ofile,quiet=TRUE,skip=1)
    if(length(q) == 0){
       return(NULL);
    }
   
    SurfArea<-q[2]
    Volume<-q[3]
    
    unlink(c(ofile,dfile))
    
    return(list(SurfArea=SurfArea,Volume=Volume))
}
    


qhull <- function(matrx,
				  joggle     = FALSE,
				  summary    = FALSE,
				  total      = TRUE,
				  verify	 = FALSE,
				  silent     = TRUE)
{
    if (missing(matrx))
    {
	cat("Matrix argument required\n")
	return(NULL)
    }
		  
    ### data matrix ###
    x <- eval(matrx)
    d <- dim(x)
    n <- d[1]
    dimension <- d[2]
    if (sum(abs(x[!is.na(x)]) == Inf) > 0)
    {
      cat("Sorry, qhull can't handle Infs\n")
      return()
    }
    if (sum(is.na(x)) > 0)
    {
	cat("Sorry, qhull can't handle NAs\n")
	return()
    }
    dfile <- tempfile("unix")
    ofile <- tempfile("unix")
    y <- c(dimension,n)
    write(y, file = dfile)
    #write(t(x),file=dfile,append=T,ncol=dimension)
    write.table(x,file=dfile,append=T,col.names=F,row.names=F)
    
    command <- paste("cat",dfile,"| qhull",sep=" ")
    if(joggle) command <- paste(command,"QJ",sep=" ")
    if(verify) command <- paste(command,"Tv",sep=" ")
    if(summary) command <- paste(command,"s",sep=" ")
    if(total) command <- paste(command,"FA",sep=" ")
    command <- paste(command,"Fx",sep=" ")
    if(silent) command <- paste(command,"Pp",sep=" ")
    command <- paste(command, "TO",ofile,sep=" ")

    system(command, FALSE)

    q <- scan(ofile,quiet=TRUE)
    if(length(q) == 0){
       return(NULL);
    }
    index <- q[2:length(q)]+1

    #clean up after yourself
    unlink(c(ofile,dfile))
    
    return(index)
}

qhull.hyperplanes <- function(matrx,
				  joggle     = FALSE,
				  summary    = FALSE,
				  total      = FALSE,
				  verify	 = FALSE,
				  silent     = TRUE)
{
    if (!missing(matrx))
    {
        ### data matrix ###
        x <- eval(matrx)
		d <- dim(x)
		n <- d[1]
		dimension <- d[2]
        if (sum(abs(x[!is.na(x)]) == Inf) > 0)
        {
          cat("Sorry, qhull can't handle Infs\n")
          return()
        }
		if (sum(is.na(x)) > 0)
		{
          cat("Sorry, qhull can't handle NAs\n")
          return()
		}
        dfile <- tempfile("unix")
		ofile <- tempfile("unix")
		y <- c(dimension,n)
        write(y, file = dfile)
		write(t(x),file=dfile,append=T,ncol=dimension)

		command <- paste("cat",dfile,"| qhull",sep=" ")

		if(joggle)
		   command <- paste(command,"QJ",sep=" ")

		if(verify)
		   command <- paste(command,"Tv",sep=" ")

		if(summary)
		   command <- paste(command,"s",sep=" ")

		if(total)
		   command <- paste(command,"FA",sep=" ")

		if(silent) 
		   command <- paste(command,"Pp",sep=" ")

		command <- paste(command,"n",sep=" ")

        command <- paste(command, "TO",ofile,sep=" ")

        system(command, FALSE)

		q <- scan(ofile,quiet=TRUE)
		if(length(q) == 0){
		   return(NULL);
		}
		index <- matrix(q[3:length(q)],byrow=T,ncol=(dimension+1))

		#clean up after yourself
		system(paste("/bin/rm",dfile,sep=" "))
		system(paste("/bin/rm",ofile,sep=" "))

		return(index)
    }
    else
      cat("Matrix argument required\n")
    NULL
}

qhull.interior <- function(y,hyperplanes,data)
{
   if(missing(hyperplanes)){
	  if(missing(data)){
	     cat("\nMust have either hyperplanes or data defined")
		 return(NULL)
	  }
      hyperplanes <- qhull.hyperplanes(data)
   }
   if(is.vector(y)){
      zz <- c(y,1)
	  ww <- hyperplanes %*% zz
	  sum(ww>1E-6) == 0
   }
   else {
	  zz <- cbind(y,rep(1,dim(y)[1]))
	  ww <- hyperplanes %*% t(zz)
	  aa <- apply(ww,2,">",1E-6)
	  apply(aa,2,sum) == 0
   }
}
