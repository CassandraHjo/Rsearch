#' Optimize trimming for best possible merging
#'
#' @description Optimizes truncation, based on base quality score
#' (\code{--truncqual} in \code{VSEARCH}) to get the best possible merging
#' results. The function searches for the best parameters by looping through
#' different parameter values for the \code{--truncqual} option.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads. See
#' Details.
#' @param minovlen The minimum overlap between the merged reads. Must be at
#' least 5. Defaults to \code{10}.
#' @param truncqual_range The truncqual values to be tested. The sequences are
#' truncated starting from the first base with the specified base quality score
#' value or lower. Defaults to \code{1} to \code{20}. Must be given as a vector
#' with numbers.
#' @param min_size The minimum copy number (size) for a given read to be
#' included in the results. Defaults to \code{1}.
#' @param maxee_rate Threshold for average expected error. Numeric value ranging
#' form \code{0.0} to \code{1.0}. Defaults to \code{1}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be
#' retained. Defaults to \code{1}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}.
#' Defaults to \code{1}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}},
#' \code{\link{vs_fastx_trim_filt}}, and \code{\link{vs_fastx_uniques}} where
#' the arguments to this functions are described in detail.
#'
#' The best possible merging option is measured by the sum of copy numbers after
#' dereplication of the merged reads bigger than 1
#' (singletons are removed when summing).
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{truncqual_value}: The truncqual value used in the trimming.
#'   \item \code{sum_size}: Sum of the copy numbers for the dereplicated
#'   sequences with copynumber above the number specified by \code{min_size}.
#'   \item \code{R1_length}: The average length of R1-reads.
#'   \item \code{R2_length}: The average length of R2-reads.
#' }
#'
#' The data frame has an attribute \code{"sum_size_plot"} containing a
#' \code{\link{ggplot2}} object based on the returned data frame. In the plot
#' the \code{truncqual} values are plotted against the \code{sum_size} values,
#' with the optimal \code{truncqual} marked in red.
#'
#' The data frame also has an attribute \code{"read_lengths_plot"} containing a
#' \code{\link{ggplot2}} object based on the returned data frame. In the plot
#' the \code{truncqual} values are plotted against the mean read lengths
#' (\code{R1_length} and \code{R2_length}).
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' \code{\link{vs_fastx_uniques}}
#'
#' @examples
#' \dontrun{
#' # Read example FASTQ files
#' R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"), "R2_sample1_small.fq")
#'
#' # Run optimizing function
#' optimize.tbl <- vs_optimize_truncqual(fastq_input = R1.file,
#'                                       reverse = R2.file)
#'
#' # Display plots
#'
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_optimize_truncqual optimize_truncqual
#'
#' @export
#'
vs_optimize_truncqual <- function(fastq_input,
                                  reverse,
                                  minovlen = 10,
                                  truncqual_range = 1:20,
                                  min_size = 1,
                                  maxee_rate = 0.01,
                                  minlen = 1,
                                  threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


  # Create data frame for storing results
  res.df <- data.frame(
    truncqual_value = truncqual_range,
    sum_size = 0,
    R1_length = 0,
    R2_length = 0
  )

  # Setting up progress bar
  pb = utils::txtProgressBar(min = 0,
                             max = length(truncqual_range),
                             initial = 0,
                             style = 3)

  # Looping through truncqual values
  for (i in 1:length(truncqual_range)) {

    # Update progress bar
    utils::setTxtProgressBar(pb, i)

    # Trim R1 and R2 reads together
    trim_R1.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                     reverse = reverse,
                                     maxee_rate = maxee_rate,
                                     minlen = minlen,
                                     truncqual = truncqual_range[i],
                                     stripright = 0,
                                     threads = threads)
    trim_R2.df <- attr(trim_R1.df, "reverse")

    # Merge R1 and R2 reads
    merge.df <- vs_fastq_mergepairs(fastq_input = trim_R1.df,
                                    reverse = trim_R2.df,
                                    minovlen = minovlen,
                                    output_format = "fasta",
                                    minlen = minlen,
                                    threads = threads)

    # Dereplicate merged reads
    derep.df <- vs_fastx_uniques(fastx_input = merge.df,
                                 output_format = "fasta",
                                 relabel_sha1 = TRUE)

    # Find number of reads with size > min_size
    tbl <- derep.df |>
      dplyr::mutate(size = stringr::str_remove(Header, ".+;size=")) |>
      dplyr::mutate(size = as.numeric(size)) |>
      dplyr::filter(size > min_size)

    # Add results to table
    res.df$sum_size[i] = sum(tbl$size)
    res.df$R1_length[i] = round(mean(nchar(trim_R1.df$Sequence)), 2)
    res.df$R2_length[i] = round(mean(nchar(trim_R2.df$Sequence)), 2)

  }
  # Close progress bar
  close(pb)

  # Make plots

  # Plot truncqual value vs sum size. Optimal truncqual value is marked with red
  optimal_truncqual <- res.df$truncqual_value[which.max(res.df$sum_size)]

  p1 <- ggplot2::ggplot(res.df, ggplot2::aes(x = truncqual_value, y = sum_size)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_point(data = subset(res.df, truncqual_value == optimal_truncqual),
                        ggplot2::aes(x = truncqual_value, y = sum_size),
                        color = "red", size = 3, shape = 18) +
    ggplot2::labs(title = "Sum size vs. Trunqqual value",
                  x = "Truncqual value",
                  y = "Sum size")

  # Plot truncqual values vs read lengths
  p2 <- ggplot2::ggplot(res.df, ggplot2::aes(x = truncqual_value)) +
    ggplot2::geom_line(ggplot2::aes(y = R1_length, color = "R1")) +
    ggplot2::geom_point(ggplot2::aes(y = R1_length, color = "R1")) +
    ggplot2::geom_line(ggplot2::aes(y = R2_length, color = "R2")) +
    ggplot2::geom_point(ggplot2::aes(y = R2_length, color = "R2")) +
    ggplot2::labs(title = "R1 and R2 length vs. Truncqual value",
         x = "Truncqual value",
         y = "Mean read length",
         color = "Read")

  # Add plots as attributes
  attr(res.df, "sum_size_plot") <- p1
  attr(res.df, "read_lengths_plot") <- p2

  return(res.df)
}
