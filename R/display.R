#' @title Open New \code{\link{Image}} Display
#'
#' @description \code{newDisplay} creates a window to display \code{\link{Image}}
#'  objects.
#'
#' @param window_name A character string representing the name of the display
#'  window (default: "Display").
#'
#' @param height An integer representing the height in pixels of the display
#'  window.
#'
#' @param width An integer representing the width in pixels of the display
#'  window.
#'
#' @return This function does not return anything.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{Image}}, \code{\link{display}}, \code{\link{destroyDisplay}}
#'
#' @examples
#' # TODO
#' @export
newDisplay <- function(window_name = "Display", height = 480, width = 640) {
  invisible(`_newDisplay`(window_name, height, width))
}

#' @title Destroy \code{\link{Image}} Display
#'
#' @aliases destroyAllDisplays
#'
#' @description \code{destroyDisplay} closes a specific existing
#'  \code{\link{Image}} display window. \code{destroyAllDisplays} all existing
#'  \code{\link{Image}} display window.
#'
#' @usage destroyDisplay(window_name)
#' destroyAllDisplays()
#'
#' @param window_name A character string representing the name of the display
#'  window (default: "Display").
#'
#' @return This function does not return anything.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{Image}}, \code{\link{newDisplay}}, \code{\link{display}}
#'
#' @examples
#' # TODO
#' @export
destroyDisplay <- function(window_name = "Display") {
  invisible(`_destroyDisplay`(window_name))
}


#' @export
destroyAllDisplays <- function() {
  invisible(`_destroyAllDisplays`())
}
