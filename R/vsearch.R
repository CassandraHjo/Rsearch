#' Set the VSEARCH executable
#'
#' @description Specify the valid command to invoke VSEARCH.
#'
#' @param vsearch_executable full path to the VSEARCH executable on your computer (text).
#'
#' @details Use this function to change the command used to invoke the external
#' software VSEARCH on this computer. When the \code{Rsearch} package is installed this
#' command is by default just \code{"vsearch"}.
#'
#' If you have a windows computer and have copied the binary \code{vsearch.exe}
#' to the folder C:/Documents/ on your computer, you update R with this
#' information by  \code{set_vsearch_executable("C:/Documents/vsearch")}.
#'
#' You may use the function \code{\link{vsearch}} to test if the command is valid.
#'
#' @return Nothing is returned, but the option \code{Rsearch.vsearch_executable}
#' is updated. The string is also saved to a file for later R sessions, i.e. you
#' only need to update this once (or if you change how you run/install VSEARCH).
#'
#' @seealso \code{\link{vsearch}}.
#'
#' @export
#'
set_vsearch_executable <- function(vsearch_executable){
  options(Rsearch.vsearch_executable = vsearch_executable)
  save(vsearch_executable, file = system.file("extdata/vsearch_executable.rds", package = "Rsearch"))
}


#' Test if VSEARCH can be executed
#'
#' @description Tests if the VSEARCH executable is a valid command.
#'
#' @details Use this function to test the command used to invoke the external
#' software VSEARCH on this computer.
#'
#' @seealso \code{\link{set_vsearch_executable}}.
#'
#' @export vsearch
#'
vsearch <- function(){
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  cat("The VSEARCH executable is:", vsearch_executable, "\n")
  ok <- vsearch_available(vsearch_executable)
  if(ok){
    cat("This is a valid command to invoke VSEARCH on this computer!\n")
  }
}


# Non-exported function to gracefully fail when vsearch_executable
# does not contain a proper command line
vsearch_available <- function(vsearch_executable){
  chr <- NULL
  ok <- try(chr <- system2(vsearch_executable, args = c("-h", "--quiet", ""), stdout = TRUE), silent = TRUE)
  if(length(grep("Error", ok[1])) > 0){
    stop("Cannot run ", vsearch_executable, " from R, use set_vsearch_executable() to set proper command to invoke vsearch")
    return(FALSE)
  } else {
    return(TRUE)
  }
}
