#' Check status of VSEARCH system call
#'
#' @description Checks the status attribute of the output from a system2 call
#' and stops execution with a detailed error message if the call failed.
#'
#' @param vsearch_output Character vector output from \code{system2()} with
#' \code{stdout = TRUE} and \code{stderr = TRUE}.
#' @param args Arguments used in the \code{system2()} call.
#'
#' @return Invisibly returns \code{NULL} if no error occurred. Otherwise, stops
#' with an informative error message.
#'
#' @noRd
check_vsearch_status <- function(vsearch_output, args) {
  status <- attr(vsearch_output, "status")

  if (!is.null(status)) {
    errmsg <- attr(vsearch_output, "errmsg")

    stop(
      paste0(
        "VSEARCH execution failed.\n",
        "Arguments used: ", paste(args, collapse = " "), "\n",
        "Exit status: ", status, "\n",
        if (!is.null(errmsg)) paste0("System message: ", errmsg, "\n") else "",
        "VSEARCH output:\n", paste(vsearch_output, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}
