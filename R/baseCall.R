#' @export
baseCall <- function(obj) {
  appDir <- system.file("GUI", "basecall", package = "prickly")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `prickly`.", call. = FALSE)
  }


  # get the name of the passed object
  object_to_change <- deparse(substitute(obj))

  # get the object from a given environment
  val <- get(object_to_change, envir = .GlobalEnv)
  # ?environment

  # Save the object as a reactive value
  assign("theStack", val, envir=globalenv())


  shiny::runApp(appDir, display.mode = "normal")
}
