#' Optimize truncation for best possible merging
#'
#' @description Optimizes truncation, based on base quality score
#' (\code{--truncqual} in \code{VSEARCH}), to get the best possible merging
#' results. The function searches for the best \code{truncqual} value by looping
#' through different parameter values.
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
#' @param plot_title Title of the resulting output plot. Defaults to
#' \code{"Optimization of Read Merging Based on Truncqual Value"}.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}},
#' \code{\link{vs_fastx_trim_filt}}, and \code{\link{vs_fastx_uniques}} where
#' the arguments to this functions are described in detail.
#'
#' The best possible truncation option (\code{truncqual}) for merging is
#' measured by the proportion of merged sequences with a copy number above the
#' number specified by \code{min_size} after dereplication.
#'
#' \code{proportion_merged_high_quality_read_pairs} represents the proportion
#' of merged high-quality read-pairs relative to the highest observed value
#' across all tested \code{truncqual} values. This normalization allows
#' comparisons across different \code{truncqual} values. A value close to 1.0
#' indicates that the merging efficiency is near its maximum, while a lower
#' value suggests suboptimal merging conditions.
#'
#'
#' Changing \code{min_size} will affect the results. A low \code{min_size} will
#' include merged sequences with a lower copy number after dereplication, and a
#' higher \code{min_size} will filter out more reads and only count
#' high-frequency merged sequences.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{truncqual_value}: The tested truncqual value.
#'   \item \code{merged_high_quality_read_pairs}: Absoulute count of
#'   successfully merged sequence pairs with a copy number above \code{min_size}
#'   after dereplication.
#'   \item \code{proportion_merged_high_quality_read_pairs}: A relative metric,
#'   calculated as the number of merged high-quality read-pairs divided
#'   by the maximum observed merged read count.
#'   \item \code{R1_length}: The average length of R1-reads after trimming.
#'   \item \code{R2_length}: The average length of R2-reads after trimming.
#' }
#'
#' The data frame has an attribute \code{"plot"} containing a
#' \code{\link{ggplot2}} object based on the returned data frame. In the plot
#' the \code{truncqual} values are plotted against the
#' \code{proportion_merged_high_quality_read_pairs} values and the mean read
#' lengths (\code{R1_length} and \code{R2_length}). The optimal \code{truncqual}
#' value is marked by a red dashed line.
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
#' # Display plot
#' print(attr(optimize.tbl, "plot"))
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
                                  threads = 1,
                                  plot_title = "Optimization of Read Merging Based on Truncqual Value"){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


  # Create data frame for storing results
  res.df <- data.frame(
    truncqual_value = truncqual_range,
    merged_high_quality_read_pairs = 0,
    proportion_merged_high_quality_read_pairs = 0,
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
                                 relabel_sha1 = TRUE) |>
      dplyr::mutate(size = stringr::str_remove(Header, ".+;size=")) |>
      dplyr::mutate(size = as.numeric(size))

    # Find number of reads with size > min_size
    tbl <- derep.df |>
      dplyr::filter(size > min_size)

    # Add results to table
    res.df$merged_high_quality_read_pairs[i] = sum(tbl$size)
    res.df$R1_length[i] = round(mean(nchar(trim_R1.df$Sequence)), 2)
    res.df$R2_length[i] = round(mean(nchar(trim_R2.df$Sequence)), 2)

  }
  # Close progress bar
  close(pb)

  # Calculate the max merged_high_quality_read_pairs for normalization
  max_merged_high_quality_read_pairs <- max(res.df$merged_high_quality_read_pairs)

  # Calculate relative merged_high_quality_read_pairs
  res.df <- res.df |>
    dplyr::mutate(proportion_merged_high_quality_read_pairs = round(merged_high_quality_read_pairs / max_merged_high_quality_read_pairs, 4))

  # Find optimal truncqual value from res.df
  optimal_truncqual <- res.df$truncqual_value[which.max(res.df$proportion_merged_high_quality_read_pairs)]

  long.df <- res.df |>
    tidyr::pivot_longer(cols = c(proportion_merged_high_quality_read_pairs, R1_length, R2_length),
                        names_to = "metric",
                        values_to = "value") |>
    dplyr::mutate(facet = dplyr::case_when(
      metric == "proportion_merged_high_quality_read_pairs" ~ "Merged High-quality read-pairs",
      metric %in% c("R1_length", "R2_length") ~ "Read Lengths",
    ))

  # Define "pretty" labels
  label_mapping <- c(
    proportion_merged_high_quality_read_pairs = "Proportion Merged",
    R1_length = "Average R1 Length",
    R2_length = "Average R2 Length"
  )

  # Make plot
  p <- ggplot2::ggplot(long.df, ggplot2::aes(x = truncqual_value, y = value)) +
    ggplot2::geom_line(ggplot2::aes(color = metric)) +
    ggplot2::geom_point(ggplot2::aes(color = metric)) +
    ggplot2::facet_wrap(~ facet, scales = "free_y", ncol = 1) +
    ggplot2::geom_vline(xintercept = optimal_truncqual, color = "red", linetype = "dashed") +
    ggplot2::labs(title = plot_title,
                  x = "Truncqual value",
                  y = "Value",
                  color = "Measurement") +
    ggplot2::scale_color_manual(values = c("proportion_merged_high_quality_read_pairs" = "#6A5ACD",
                                           "R1_length" = "#008080",
                                           "R2_length" = "#FF7F50"),
                                labels = label_mapping)

  # Add plot as attribute
  attr(res.df, "plot") <- p

  return(res.df)
}
