# Some useful extensions / changes to rgl defaults

#' Set up pan call back for current rgl device
#'
#' Copied verbatim from ?rgl.setMouseCallbacks for rgl version 0.92.892
#' Mouse button 2 is right and button 3 is middle (accessed by meta/alt key)
#' @param button Integer from 1 to 3 indicating mouse button
#' @export
#' @seealso \code{\link{rgl.setMouseCallbacks}}
#' @author Duncan Murdoch
#' @examples
#' \dontrun{
#'  open3d()
#'  pan3d(2)
#' }
pan3d <- function(button) {
  start <- list()
  begin <- function(x, y) {
    start$userMatrix <<- par3d("userMatrix")
    start$viewport <<- par3d("viewport")
    start$scale <<- par3d("scale")
    start$projection <<- rgl.projection()
    start$pos <<- rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4],
       0.5, projection=start$projection)
  }
  update <- function(x, y) {
    xlat <- (rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4],
        0.5, projection = start$projection) - start$pos)*start$scale
    mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
    par3d(userMatrix = start$userMatrix %*% t(mouseMatrix) )
  }
  rgl.setMouseCallbacks(button, begin, update)
  # cat("Callbacks set on button", button, "of rgl device",rgl.cur(),"\n")
}

#' Open customised rgl window
#'
#' Pan with right button (Ctrl+click), Zoom with middle (Alt/Meta+click)
#' Defaults to dark grey background and orthogonal projection (FOV=0)
#' @param bgcol Background colour
#' @param FOV Field of View
#' @param ... additional options passed to open3d
#' @return Current rgl device
#' @export
#' @seealso \code{\link{open3d},\link{pan3d}}
open3dgj<- function(bgcol='dark grey', FOV=0, ...){
  res=open3d(mouseMode=c("trackball","user","zoom"), FOV=FOV, ...)
  bg3d(col=bgcol)
  pan3d(2)
  res
}
