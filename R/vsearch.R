#' Set the VSEARCH executable
#'
#' @param executable full path to the VSEARCH executable on your computer.
#'
#' @export set_vsearch_executable
#'
set_vsearch_executable <- function(vsearch_executable){
  options(Rsearch.vsearch_executable = vsearch_executable)
  save(vsearch_executable, file = system.file("data/vsearch_executable.rds", package = "Rsearch"))
}


# Non-exported function to gracefully fail when vsearch_executable
# does not contains a proper command line
vsearch_available <- function(){
  chr <- NULL
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  ok <- try(chr <- system2(vsearch_executable, stdout = TRUE), silent = TRUE)
  if(length(grep("Error", ok[1])) > 0){
    stop("Cannot run ", vsearch_executable, " from R, use set_vsearch_executable()\nto set proper command to invoke vsearch")
    return(FALSE)
  } else {
    return(TRUE)
  }
}
