# This script is designed to be almost invariant and hosted at 
# an http URL so that it can be sourced in R.
# (R cannot directly source the https urls that are now required by github)

if(!require(httr)){
  message("Installing httr package to permit https downloads (from github)")
  install.packages("httr")
  if(!require(httr)) {
    browseURL("https://github.com/jefferis/AnalysisSuite/blob/master/README.md")
    stop("Manual install required. See section 'R Install', subsection 'Manual Install' of README.md")
  }
}

local({
  install_AS<-function(){
    req=GET('https://raw.github.com/jefferis/AnalysisSuite/master/R/Install.R')
    tc=textConnection(text_content(req))
    on.exit(close(tc))
    source(tc)
  }
  install_AS()
})