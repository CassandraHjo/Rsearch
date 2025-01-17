#' Optimize trimming for best possible merging
#'
#' @description Optimizes truncation based on base quality score (\code{--truncqual}
#' in \code{VSEARCH}) to get the best possible merging results. The function
#' searches for the best parameters by looping through different parameter values
#' for the \code{--truncqual} option.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads See Details.
#' @param minovlen The minimum overlap between the merged reads. Must be at least 5. Defaults to \code{10}.
#' @param truncqual_range The truncqual values to be tested. The sequences are
#' truncated starting from the first base with the specified base quality score
#' value or lower. Defaults to \code{1} to \code{20}. Must be given as a vector with numbers.
#' @param maxee_rate Threshold for average expected error. Numeric value ranging form
#' \code{0.0} to \code{1.0}. Defaults to \code{1}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{1}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' and \code{\link{vs_fastx_uniques}} where the arguments to this functions
#' are described in detail.
#'
#' The best possible merging option is measured by the sum of copy numbers after
#' dereplication of the merged reads bigger than 1
#' (singletons are removed when summing).
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{truncqual_value}: The truncqual value used in the trimming.
#'   \item \code{sum_size}: Sum of the copy numbers for the dereplicated sequences with copynumber above 1.
#'   \item \code{R1}: The name of the R1 table/file.
#'   \item \code{R2}: The name of the R2 table/file.
#' }
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' \code{\link{vs_fastx_uniques}}
#'
#' @export
#'
vs_optimize_truncqual <- function(fastq_input,
                                  reverse,
                                  minovlen = 10,
                                  truncqual_range = paste(1:20),
                                  maxee_rate = 0.01,
                                  minlen = 1,
                                  threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


  # Create empty data frame for storing results
  res.df <- data.frame(
    truncqual_value = numeric(0),
    sum_size = numeric(0),
    R1 = character(),
    R2 = character()
  )

  # Setting up progress bar
  pb = utils::txtProgressBar(min = 0,
                             max = length(truncqual_range),
                             initial = 0,
                             style = 3)

  # Create counting variable
  stepi <- 0

  # Looping through truncqual values

  for (val in truncqual_range) {

    # Update counting variable and progress bar
    stepi <- stepi + 1
    utils::setTxtProgressBar(pb,stepi)

    # Trim R1 and R2 reads together
    trim_R1.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                     reverse = reverse,
                                     maxee_rate = maxee_rate,
                                     minlen = minlen,
                                     truncqual = val,
                                     stripright = 0,
                                     threads = threads)

    # Merge R1 and R2 reads
    merge.df <- vs_fastq_mergepairs(fastq_input = trim_R1.df,
                                    reverse = attr(trim_R1.df, "reverse"),
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
      truncqual_value = val,
      sum_size = sum(tbl$size),
      R1 = as.character(substitute(fastq_input)),
      R2 = as.character(substitute(reverse))
    )

    res.df <- rbind(new_row, res.df)

  }
  # Close progress bar
  close(pb)

  # Sort results table
  res.df <- res.df |>
    dplyr::arrange(truncqual_value)

  return(res.df)
}
