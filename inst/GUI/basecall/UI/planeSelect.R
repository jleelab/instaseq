bsCollapsePanel(
  title = actionLink("toggleTracking", "Segmentation module"),
  value = "trackingPanel",
  sliderInput("range", "Stacks to process:",
              min = 1, max = theStack$numFiles()/4,
              value = c(1, theStack$numFiles()/4)),
  checkboxInput("showMIP", "Display max. intensity projection", FALSE),
  actionButton("computeSeg", "Run stack segmentation"),
  htmlOutput("trackingStatus")
)
