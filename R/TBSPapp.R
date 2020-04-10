#' Run the TBSignatureProfiler Shiny application.
#'
#' Use this function to run the TBSignatureProfiler application.
#'
#'
#' @return The Shiny application will open.
#'
#' @export
#'
#' @examples
#' # Upload data through the app
#' if (interactive()){
#'   TBSPapp()
#'   }
#'
TBSPapp <- function() {
  appDir <- system.file("shiny", package = "TBSignatureProfiler")
  shiny::runApp(appDir, display.mode = "normal")
}
