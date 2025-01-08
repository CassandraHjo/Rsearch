#' Optimize trimming for best possible merging
#'
#' @description Optimizes trimming to get the best possible merging results. The function searches for the best parameters by looping through different parameter values for the trimming options.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads See Details.
#' @param minovlen The minimum overlap between the merged reads. Must be at least 5. Defaults to \code{10}.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}} and \code{\link{vs_fastx_trim_filt}} where the arguments to this function are described in detail.
#'
#' @returns x
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}}
#' @export
#'
vs_optimize_trim <- function(fastq_input,
                             reverse,
                             minovlen = 10,
                             minlen = 0,
                             threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


}
