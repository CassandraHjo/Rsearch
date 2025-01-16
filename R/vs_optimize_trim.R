#' Optimize trimming for best possible merging
#'
#' @description Optimizes trimming to get the best possible merging results. The function searches for the best parameters by looping through different parameter values for the trimming options.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads See Details.
#' @param stripright_R1_range The stripright values to be tested. The number of
#' bases stripped from the right end of the forward reads. Defaults to \code{c(0, 10, 50)}.
#' The range must be given as a vector with numbers.
#' @param stripright_R2_range The stripright values to be tested. The number of
#' bases stripped from the right end of the reverse reads. Defaults to \code{c(0, 10, 50)}.
#' The range must be given as a vector with numbers.
#' @param minovlen The minimum overlap between the merged reads. Must be at least 5. Defaults to \code{10}.
#' @param truncqual The sequences are truncated starting from the first base with the
#' specified base quality score value or lower. Defaults to \code{1}.
#' @param maxee_rate Threshold for average expected error. Numeric value ranging form
#' \code{0.0} to \code{1.0}. Defaults to \code{1}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{1}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}},
#' \code{\link{vs_fastx_trim_filt}}, \code{\link{fastx_synchronize}}, and
#' \code{\link{vs_fastx_uniques}} where the arguments to this functions are
#' described in detail.
#'
#' The best possible merging option is measured by the sum of copy numbers after
#' dereplication of the merged reads bigger than 1
#' (singletons are removed when summing).
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{stripright_R1}: The stripright value used on the R1 reads in the trimming.
#'   \item \code{stripright_R2}: The stripright value used on the R2 reads in the trimming.
#'   \item \code{sum_size}: Sum of the copy numbers for the dereplicated sequences with copynumber above 1.
#'   \item \code{R1}: The name of the R1 table/file.
#'   \item \code{R2}: The name of the R2 table/file.
#' }
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' \code{\link{fastx_synchronize}}, \code{\link{vs_fastx_uniques}}
#'
#' @export
#'
vs_optimize_trim <- function(fastq_input,
                             reverse,
                             stripright_R1_range = c(10, 20, 50),
                             stripright_R2_range = c(10, 20, 50),
                             minovlen = 10,
                             truncqual = 1,
                             maxee_rate = 0.01,
                             minlen = 1,
                             threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


  # Create empty data frame for storing results
  res.df <- data.frame(
    stripright_R1 = numeric(0),
    stripright_R2 = numeric(0),
    sum_size = numeric(0),
    R1 = character(),
    R2 = character()
  )

  # Setting up progress bar
  pb = utils::txtProgressBar(min = 0,
                             max = length(stripright_R1_range) * length(stripright_R1_range),
                             initial = 0,
                             style = 3)

  # Create counting variable
  stepi <- 0

  # Looping
  for (R1.val in stripright_R1_range) {
    for (R2.val in stripright_R2_range){

      # Update counting variable and progress bar
      stepi <- stepi + 1
      utils::setTxtProgressBar(pb,stepi)

      # Trim R1 reads
      trim_R1.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                       reverse = NULL,
                                       maxee_rate = maxee_rate,
                                       minlen = minlen,
                                       truncqual = truncqual,
                                       stripright = R1.val,
                                       threads = threads)

      # Trim R2 reads
      trim_R2.df <- vs_fastx_trim_filt(fastx_input = reverse,
                                       reverse = NULL,
                                       maxee_rate = maxee_rate,
                                       minlen = minlen,
                                       truncqual = truncqual,
                                       stripright = R2.val,
                                       threads = threads)

      # Sync R1 and R2 files
      sync_R1 <- fastx_synchronize(file1 = trim_R1.df,
                                   file2 = trim_R2.df)
      sync_R2 <- attr(sync_R1, "reverse")

      # Merge R1 and R2 reads
      merge.df <- vs_fastq_mergepairs(fastq_input = sync_R1,
                                      reverse = sync_R2,
                                      minovlen = minovlen,
                                      output_format = "fasta",
                                      minlen = minlen,
                                      threads = threads)

      # Dereplicate merged reads
      derep.df <- vs_fastx_uniques(fastx_input = merge.df,
                                   output_format = "fasta",
                                   relabel_sha1 = TRUE)

      # Find number of reads with size > 1
      tbl <- derep.df |>
        dplyr::mutate(size = stringr::str_remove(Header, ".+;size=")) |>
        dplyr::mutate(size = as.numeric(size)) |>
        dplyr::filter(size > 1)

      # Add results to table
      new_row <- data.frame(
        stripright_R1 = R1.val,
        stripright_R2 = R2.val,
        sum_size = sum(tbl$size),
        R1 = as.character(substitute(fastq_input)),
        R2 = as.character(substitute(reverse))
      )

      res.df <- rbind(new_row, res.df)
    }
  }

  # Close progress bar
  close(pb)

  # Sort results table
  res.df <- res.df |>
    dplyr::arrange(stripright_R1)

  return(res.df)
}
