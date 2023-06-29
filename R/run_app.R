
#' @export
run_app <- function() {
  shiny::shinyAppDir(system.file("shiny", package = "TrialExpertElicitation"))
}
