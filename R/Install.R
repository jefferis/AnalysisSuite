# Easy Installation script from github repository
# Also 
local({
  install_bioc<-function(x){
    source("http://bioconductor.org/biocLite.R")
    for (pkg in x) biocLite(pkg)
  }
  
  install_or_update_cran_package<-function(x){
    allpackages=.packages(all = TRUE)
    new=setdiff(x,allpackages)
    existing=intersect(x,allpackages)
    if(length(existing))
      update.packages(existing)
    if(length(new))
      install.packages(new)
  }
  
  gjanalysis_suite_install<-function(local_file=FALSE){
    message("Please choose a location in which to INSTALL AnalysisSuite")
    message("You will need to type a fake filename in the dialog box you are presented")
    flush.console()
    installpath=dirname(file.choose(new=TRUE))
    if(local_file){
      message("Please choose the location of downloaded AnalysisSuite zip file")
      flush.console()
      zip_file=file.choose()      
    } else {
      zip_file=paste(tempfile(),'-analysissuite.zip',sep="")
      message("Downloading AnalysisSuite from github")
      zip_url='https://github.com/jefferis/AnalysisSuite/zipball/master'
      request=GET(zip_url, config(ssl.verifypeer = FALSE))
      stop_for_status(request)
      writeBin(content(request), zip_file)
      on.exit(unlink(zip_file),add=TRUE)
      zf=unzip(zip_file,exdir=installpath)
      # this will include the name of the zip folder itself
      dirname(zf[1])
    }
  }
  
  set_option<-function(p){
    options(gjanalysissuite.startup=p)
  }
  
  add_path_to_rprofile<-function(p){
    rp=file.path(normalizePath('~',mustWork=TRUE),'.RProfile')
    optionline=paste("options(gjanalysissuite.startup='",p,"')",sep='')
    cat(file=rp,append=TRUE,optionline,sep="",'\n')
  }
  RDir<-dirname(attr(body(function() {}),'srcfile')$filename)
  if(RDir=='.'){
    # we are sourcing this from a connection (which happens during web install)
    installpath=gjanalysis_suite_install()
  } else {
    # file has already been unzipped (hopefully in the right place)
    # 
    installpath=dirname(RDir)
  }
  message("Installing required packages")
  install_or_update_cran_package(c("rgl",'RANN','igraph'))
  install_bioc('RBGL')
  startup=file.path(installpath,'R','Code','Startup.R')
  set_option(startup)
  add_path_to_rprofile(startup)
  message("********************* Install completed ********************")
  message("You can start AnalysisSuite (now and in the future) by doing")
  message("source(options()$gjanalysissuite.startup)")
  message("************************************************************")
})
