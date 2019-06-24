shinyServer(function(input, output, session) {

  theActive <- reactiveVal("tracking")
  # thePath <- reactiveVal()
  # theFile <- reactiveVal()
  # theVideo <- reactiveVal()
  # theImage <- reactiveVal()
  # theBackground <- reactiveVal()
  # theMask <- reactiveVal(NULL)
  # theBlobs <- reactiveVal()
  # theBlobSizes <- {}
  # theSuccess <- reactiveVal(FALSE)
  #
  #newDisplay("trackR")
  theStack$render()
  #
  # source("SERVER/video.R", local = TRUE)
  #
  # source("SERVER/background.R", local = TRUE)
  #
  # source("SERVER/mask.R", local = TRUE)
  #
  # source("SERVER/blob.R", local = TRUE)
  #
  source("SERVER/planeSelect.R", local = TRUE)
  #
  output$videoSlider <- renderUI({
    if (TRUE) {

      sliderInput("stackPos", label = h3("z-position"), width = "100%", value = as.integer(theStack$numFiles()/4/2), min = 1,
                  max = theStack$numFiles()/4, step = 1)


    } else {
      sliderInput("videoPos", "Frame", width = "100%", value = 1, min = 1,
                  max = 1, step = 1)
    }
  })

  # Clean up
  session$onSessionEnded(function() {
    destroyAllDisplays()
  })
})
