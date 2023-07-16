library(shiny)
library(profvis)

profvis({
  shiny::runApp(appDir = "inst/shiny", display.mode = "normal")
})
