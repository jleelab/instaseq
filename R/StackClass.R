#' @title An S4 Class Containing an OpenCV Image
#'
#' @name Stack-class
#'
#' @aliases Rcpp_Stack
#'
#' @docType class
#'
#' @description \code{Image} objects are the base objects of the \pkg{\link{Rvision}}
#'  package. They contain an \href{http://opencv.org/}{OpenCV} image that can
#'  originate from an image file, an array, a video file or a video stream.
#'  This image can be manipulated using the functions of \pkg{\link{Rvision}}.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{image}}, \code{\link{Video}}, \code{\link{Stream}}
#' @export
"Stack"


#' @title Create an Object of Class \code{Image}
#'
#' @description Function for creating \code{\link{Image}} objects from arrays
#'  and image files.
#'
#' @param ... When created from an image file, \code{image} takes one argument
#'  that is a character string indicating the path to the image file. When
#'  created from an array (e.g. a matrix), it takes this array as its single
#'  argument. An \code{Image} object can also be created without any argument,
#'  in which case it is empty and can be populated with an image later.
#'
#' @return An \code{\link{Image}} object.
#'
#' @note \code{Image} objects can be created from video files and video streams
#'  using the following functions: \code{\link{video}}, \code{\link{stream}}.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{Image}}
#'
#' @examples
#' # TODO
#'
#' @export
stack <- function(input, acq = 'ch.zplane') {
  if(acq == 'ch.zplane'){
    d<-gsub('\\D+','', basename(input))
    num<-gsub("(?<![0-9])0+", "", d, perl = TRUE)
    input<-input[order(as.numeric(num))]
    channels<-sub(".*c", "", basename(input))

    input<-input[order(channels)]
  }

  new(prickly::Stack, input)
}

#' @title Test for an Image Object
#'
#' @description Tests whether the object is of class \code{\link{Image}}
#'
#' @param object Any R object.
#'
#' @return A logical indicating whether the object is of class
#'  \code{\link{Image}} (TRUE) or not (FALSE).
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{Image}}, \code{\link{image}}
#'
#' @examples
#' # TODO
#'
#' @export
isStack <- function(object) {
  inherits(object, "Stack") & (tryCatch(object$ncol(), error = function(e) 0) > 0)
}


#' @title The Number of Rows/Columns/Channels of an Image
#'
#' @aliases nrow.Image ncol.Rcpp_Image ncol.Image nchan
#'
#' @description nrow, ncol and nchan return the number of rows, columns or
#'  channels present in an \code{\link{Image}} object.
#'
#' @usage nrow.Rcpp_Image(x)
#' ncol.Rcpp_Image(x)
#' nchan(x)
#'
#' @param x An \code{\link{Image}} object.
#'
#' @return A numeric value.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link[=dim.Rcpp_Image]{dim}} which returns \emph{all}
#'  dimensions.
#'
#' @examples
#' # TODO
#' @export
pop.Rcpp_Stack <- function(x) {
  cat('Hello!')
  if (!isStack(x))
    stop("This is not a Stack object.")
  x$display()
}
