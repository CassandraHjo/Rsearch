#' Optimize truncation with truncee_rate for optimal read merging
#'
#' @description \code{vs_optimize_truncee_rate} optimizes the truncation
#' parameter \code{truncee_rate} to achieve the best possible merging results.
#' The function iterates through a specified range of \code{truncee_rate} values
#' to identify the optimal value that maximizes the proportion of high-quality
#' merged read pairs.
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing (forward)
#' reads. See \emph{Details}.
#' @param reverse A FASTQ file path or FASTQ object containing (reverse) reads.
#' See \emph{Details}.
#' @param minovlen Minimum overlap between the merged reads. Must be at least 5.
#' Defaults to \code{10}.
#' @param truncee_rate_range A numeric vector of \code{truncee_rate} values to
#' test. Sequences are truncated starting from the first base with the specified
#' base quality score or lower. Defaults to
#' \code{c(0.0005, 0.001, seq(0.002, 0.014, by = 0.002)}. Provide as a numeric
#' vector.
#' @param minlen Minimum number of bases a sequence must have to be retained.
#' Defaults to \code{0}. See \emph{Details}.
#' @param min_size Minimum copy number (size) for a merged read to be
#' included in the results. Defaults to \code{2}.
#' @param maxee_rate Threshold for average expected error. Must range from
#' \code{0.0} to \code{1.0}. Defaults to \code{0.01}. See\emph{Details}.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param plot_title If \code{TRUE} (default), a summary title will be displayed
#' in the plot. Set to \code{FALSE} for no title.
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}},
#' \code{\link{vs_fastx_trim_filt}}, and \code{\link{vs_fastx_uniques}} where
#' the arguments to this functions are described in detail.
#'
#' The best possible truncation option (\code{truncee_rate}) for merging is
#' measured by the proportion of merged sequences with a copy number above the
#' number specified by \code{min_size} after dereplication.
#'
#' \code{merged_read_pairs} represents the proportion
#' of merged high-quality read-pairs relative to the highest observed value
#' across all tested \code{truncee_rate} values. This normalization allows
#' comparisons across different \code{truncee_rate} values. A value close to 1.0
#' indicates that the merging efficiency is near its maximum, while a lower
#' value suggests suboptimal merging conditions.
#'
#' Changing \code{min_size} will affect the results. A low \code{min_size} will
#' include merged sequences with a lower copy number after dereplication, and a
#' higher \code{min_size} will filter out more reads and only count
#' high-frequency merged sequences.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{truncee_rate_value}: Tested \code{truncee_rate} value.
#'   \item \code{merged_high_quality_read_pairs}: Absolute count of
#'   successfully merged sequence pairs with a copy number above \code{min_size}
#'   after dereplication.
#'   \item \code{merged_read_pairs}: A relative metric,
#'   calculated as the number of merged high-quality read-pairs divided
#'   by the maximum observed merged read count.
#'   \item \code{R1_length}: Average length of R1-reads after trimming.
#'   \item \code{R2_length}: Average length of R2-reads after trimming.
#' }
#'
#' The returned data frame has an attribute named \code{"plot"} containing a
#' \code{\link{ggplot2}} object based on the returned data frame. The plot
#' visualizes \code{truncee_rate} values against
#' \code{merged_read_pairs}, \code{R1_length}, and
#' \code{R2_length}, with the optimal \code{truncee_rate} value marked by a red
#' dashed line.
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' \code{\link{vs_fastx_uniques}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "R1_sample1_small.fq")
#' R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "R2_sample1_small.fq")
#'
#' # Run optimizing function
#' optimize.tbl <- vs_optimize_truncee_rate(fastq_input = R1.file,
#'                                       reverse = R2.file)
#'
#' # Display plot
#' print(attr(optimize.tbl, "plot"))
#'
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_optimize_truncee_rate optimize_truncee_rate
#'
#' @export
#'
vs_optimize_truncee_rate <- function(fastq_input,
                                     reverse,
                                     minovlen = 10,
                                     truncee_rate_range = c(0.0005,
                                                            0.001,
                                                            seq(0.002, 0.014,
                                                                by = 0.002)
                                     ),
                                     minlen = 1,
                                     min_size = 2,
                                     maxee_rate = 0.01,
                                     threads = 1,
                                     plot_title = TRUE){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)


  # Create data frame for storing results
  res.df <- data.frame(
    truncee_rate_value = truncee_rate_range,
    merged_read_pairs = 0,
    R1_length = 0,
    R2_length = 0
  )

  # Setting up progress bar
  pb = utils::txtProgressBar(min = 0,
                             max = length(truncee_rate_range),
                             initial = 0,
                             style = 3)

  # Looping through truncee_rate values
  for (i in 1:length(truncee_rate_range)) {

    # Update progress bar
    utils::setTxtProgressBar(pb, i)

    # Trim R1 and R2 reads together
    trim_R1.df <- vs_fastx_trim_filt(fastx_input = fastq_input,
                                     reverse = reverse,
                                     maxee_rate = maxee_rate,
                                     minlen = minlen,
                                     truncee_rate = truncee_rate_range[i],
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

    # Find number of dereplicated merged reads with size > min_size
    derep.df_filt <- derep.df |>
      dplyr::filter(size > min_size)

    # Add results to table
    res.df$merged_read_pairs[i] = sum(derep.df_filt$size)
    res.df$R1_length[i] = round(mean(nchar(trim_R1.df$Sequence)), 2)
    res.df$R2_length[i] = round(mean(nchar(trim_R2.df$Sequence)), 2)

  }
  # Close progress bar
  close(pb)

  # Find optimal truncee_rate value from res.df
  optimal_truncee_rate <- res.df$truncee_rate_value[which.max(res.df$merged_read_pairs)]

  long.df <- res.df |>
    tidyr::pivot_longer(cols = c(merged_read_pairs, R1_length, R2_length),
                        names_to = "metric",
                        values_to = "value") |>
    dplyr::mutate(facet = dplyr::case_when(
      metric == "merged_read_pairs" ~ "Merged read-pairs",
      metric %in% c("R1_length", "R2_length") ~ "Read Lengths",
    ))

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  # Make plot for merging proportion
  p1 <- ggplot2::ggplot(long.df[long.df$facet == "Merged read-pairs", ],
                        ggplot2::aes(x = truncee_rate_value, y = value, color = metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = optimal_truncee_rate, color = "red", linetype = "dashed") +
    ggplot2::labs(title = "Merged read-pairs",
                  x = "Truncee_rate value",
                  y = "Number of read-pairs",
                  color = "") +
    ggplot2::scale_color_manual(values = c("merged_read_pairs" = pal[2]),
                                labels = c(
                                  merged_read_pairs = "Number of read-pairs merged")) +
    ggplot2::theme_minimal() +
    # Remove x-axis because this is common with p2
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  # Make plot for read lengths
  p2 <- ggplot2::ggplot(long.df[long.df$facet == "Read Lengths", ],
                        ggplot2::aes(x = truncee_rate_value, y = value, color = metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = optimal_truncee_rate, color = "red", linetype = "dashed") +
    ggplot2::labs(title = "Read Lengths",
                  x = "truncee_rate",
                  y = "Length (bases)",
                  color = "") +
    ggplot2::scale_color_manual(values = c("R1_length" = pal[3],
                                           "R2_length" = pal[4]),
                                labels = c(R1_length = "Average R1 Length",
                                           R2_length = "Average R2 Length"
                                )) +
    ggplot2::theme_minimal()

  # Combine the two plots
  combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v")

  # Create the main title
  if (plot_title) {
    title <- paste(max(res.df$merged_read_pairs),
                   "read-pairs merged with truncee rate:",
                   optimal_truncee_rate,
                   "(total:",
                   num_readpairs,
                   ", size >",
                   min_size,
                   ")"
    )
  } else {
    title <- ""
  }

  # "Draw" the main title
  common_title <- cowplot::ggdraw() +
    cowplot::draw_label(title, size = 12, x = 0.01, hjust = 0)

  # Combine title and plot
  final_plot <- cowplot::plot_grid(common_title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

  # Add plot as attribute
  attr(res.df, "plot") <- final_plot

  return(res.df)
}
