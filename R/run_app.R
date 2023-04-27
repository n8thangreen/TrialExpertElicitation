#' @export
run_app <- function() {
  shiny::runApp(appDir = system.file("shiny", package = "TrialExpertElicitation"))
}
