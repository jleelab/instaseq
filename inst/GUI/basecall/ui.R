shinyUI(fluidPage(

  titlePanel("base call"),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),

  fluidRow(
    column(
      12,

      bsCollapse(

        id = "bsCollapse",
        multiple = TRUE,
        open = "videoPanel",

        # source("UI/video.R", local = TRUE)$value,
        #
        # source("UI/background.R", local = TRUE)$value,
        #
        # source("UI/mask.R", local = TRUE)$value,
        #
        # source("UI/blob.R", local = TRUE)$value,
        #
        source("UI/planeSelect.R", local = TRUE)$value

      ),

      bsCollapse(
        open = "controlPanel",
        bsCollapsePanel(
          title = NULL,
          value = "controlPanel",
          htmlOutput("videoSlider")
        )
      )
    )
  )

))
