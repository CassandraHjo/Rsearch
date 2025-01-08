#' Optimize trimming for best possible merging
#'
#' @description Optimizes trimming to get the best possible merging results. The function searches for the best parameters by looping through different parameter values for the trimming options.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads See Details.
#' @param minovlen The minimum overlap between the merged reads. Must be at least 5. Defaults to \code{10}.
#' @param truncqual The sequences are truncated starting from the first base with the
#' specified base quality score value or lower. Defaults to \code{20}.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}} and \code{\link{vs_fastx_trim_filt}} where the arguments to this function are described in detail.
#'
#' @returns Results table
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}}
#' @export
#'
vs_optimize_trim <- function(fastq_input,
                             reverse,
                             minovlen = 10,
                             truncqual = 20,
                             minlen = 0,
                             threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Create empty data frame for storing results
  res.df <- data.frame(
    stripright_R1 = numeric(0),
    stripright_R2 = numeric(0),
    num_merged = numeric(0),
    R1_file = character(),
    R2_file = character()
  )

  # Decide which stripping values to use
  stripright_R1_values <- c(0, 10, 20, 50, 80, 100)
  stripright_R2_values <- c(0, 10, 20, 50, 80, 100)

  # Looping
  for (R1.val in stripright_R1_values) {
    for (R2.val in stripright_R2_values){

      # Trim R1 reads
      trim_R1.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                       reverse = NULL,
                                       maxee_rate = NULL,
                                       minlen = minlen,
                                       truncqual = truncqual,
                                       truncee = NULL,
                                       stripright = R1.val,
                                       threads = threads)

      # Trim R2 reads
      trim_R2.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                       reverse = NULL,
                                       maxee_rate = NULL,
                                       minlen = minlen,
                                       truncqual = truncqual,
                                       truncee = NULL,
                                       stripright = R2.val,
                                       threads = threads)

      # Sync R1 and R2 files
      sync_R1 <- fastx_synchronize(file1 = trim_R1.df,
                                   file2 = trim_R2.df)
      sync_R2 <- attr(sync_R1, "sync_file2")

      # Merge R1 and R2 reads
      merge.df <- vs_merging_lengths(fastq_input = sync_R1,
                                     reverse = sync_R2,
                                     minovlen = minovlen,
                                     minlen = minlen,
                                     threads = threads)

      # Find number of merged read pairs
      num_merged <- sum(!is.na(merge.df$length_overlap))

      # Add results to table
      new_row <- data.frame(
        stripright_R1 = R1.val,
        stripright_R2 = R2.val,
        num_merged = num_merged,
        R1_file = basename(fastq_input),
        R2_file = basename(reverse)
      )

      res.df <- rbind(new_row, res.df)

    }
  }
  return(res.df)
}
