# Toggle panel
observeEvent(input$toggleTracking, {
  if ("trackingPanel" %in% input$bsCollapse) {
    updateCollapse(session, "bsCollapse", close = "trackingPanel")
  } else {
    updateCollapse(session, "bsCollapse", open = "trackingPanel",
                   close = c("backgroundPanel", "blobPanel", "maskPanel", "videoPanel"))
    theActive("tracking")
  }
})


observeEvent(input$stackPos, {
  theStack$updateDisp(input$stackPos-1, 0, 0)
})
