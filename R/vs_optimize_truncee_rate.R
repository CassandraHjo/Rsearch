#' Optimize read truncation with truncee_rate
#'
#' @description \code{vs_optimize_truncee_rate} optimizes the truncation
#' parameter \code{truncee_rate} to achieve the best possible merging results.
#' The function iterates through a specified range of \code{truncee_rate} values
#' to identify the optimal value that maximizes the proportion of high-quality
#' merged read pairs.
#'
#' @param fastq_input (Required). A FASTQ file path, FASTQ tibble (forward
#' reads), or a paired-end tibble of class \code{"pe_df"}. See \emph{Details}.
#' @param reverse (Optional). A FASTQ file path or FASTQ tibble (reverse reads).
#' Optional if \code{fastq_input} is a \code{"pe_df"} object.
#' @param minovlen (Optional). Minimum overlap between the merged reads. Must be
#' at least 5. Defaults to \code{10}.
#' @param truncee_rate_range (Optional). A numeric vector of \code{truncee_rate}
#' values to test. Sequences are truncated so that their average expected error
#' per base is lower than the specified value. Defaults to \code{(0.002, 0.004,
#' 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024, 0.026,
#' 0.028, 0.030, 0.032, 0.034, 0.036, 0.038, 0.040)}.
#' @param minlen (Optional). Minimum number of bases a sequence must have to be
#' retained. Defaults to \code{0}. See \emph{Details}.
#' @param min_size (Optional). Minimum copy number (size) for a merged read to
#' be included in the results. Defaults to \code{2}.
#' @param maxee_rate (Optional). Threshold for average expected error. Must
#' range from \code{0.0} to \code{1.0}. Defaults to \code{0.01}. See
#' \emph{Details}.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param plot_title (Optional). If \code{TRUE} (default), a summary title will
#' be displayed in the plot. Set to \code{FALSE} for no title.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details
#' The function uses \code{\link{vs_fastq_mergepairs}},
#' \code{\link{vs_fastx_trim_filt}}, and \code{\link{vs_fastx_uniques}} where
#' the arguments to this functions are described in detail.
#'
#' If \code{fastq_input} has class \code{"pe_df"}, the reverse reads will be
#' automatically extracted from the \code{"reverse"} attribute unless
#' explicitly provided in the \code{reverse} argument.
#'
#' The best possible truncation option (\code{truncee_rate}) for merging is
#' measured by the number of merged read-pairs with a copy number above the
#' number specified by \code{min_size} after dereplication.
#'
#' Changing \code{min_size} will affect the results. A low \code{min_size} will
#' include merged sequences with a lower copy number after dereplication, and a
#' higher \code{min_size} will filter out more reads and only count
#' high-frequency merged sequences.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{truncee_rate_value}: Tested \code{truncee_rate} value.
#'   \item \code{merged_read_pairs}: Count of merged read-pairs with a copy
#'   number above \code{min_size} after dereplication.
#'   \item \code{R1_length}: Average length of R1-reads after trimming.
#'   \item \code{R2_length}: Average length of R2-reads after trimming.
#' }
#'
#' The returned data frame has an attribute named \code{"plot"} containing a
#' \code{\link[ggplot2]{ggplot2}} object based on the returned data frame. The
#' plot visualizes \code{truncee_rate} values against \code{merged_read_pairs},
#' \code{R1_length}, and \code{R2_length}, with the optimal \code{truncee_rate}
#' value marked by a red dashed line.
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}, \code{\link{vs_fastx_trim_filt}},
#' \code{\link{vs_fastx_uniques}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R1.fq")
#' R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R2.fq")
#'
#' # Run optimizing function
#' optimize.tbl <- vs_optimize_truncee_rate(fastq_input = R1.file,
#'                                          reverse = R2.file)
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
                                     reverse = NULL,
                                     minovlen = 10,
                                     truncee_rate_range = c(
                                       seq(0.002, 0.04,
                                           by = 0.002)),
                                     minlen = 1,
                                     min_size = 2,
                                     maxee_rate = 0.01,
                                     threads = 1,
                                     plot_title = TRUE,
                                     tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Handle pe_df class
  if (is_pe_df(fastq_input) && is.null(reverse)) {
    reverse <- attr(fastq_input, "reverse")
    if (is.null(reverse)) {
      stop("fastq_input has class 'pe_df' but no 'reverse' attribute found.")
    }
  }

  # Create data frame for storing results
  res.df <- data.frame(
    truncee_rate_value = truncee_rate_range,
    merged_read_pairs = 0,
    R1_length = 0,
    R2_length = 0
  )

  # Get the number of read pairs
  if (!is.character(fastq_input)) {
    num_readpairs <- nrow(fastq_input)
  } else {
    num_readpairs <- nrow(microseq::readFastq(fastq_input))
  }

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
                                     threads = threads,
                                     tmpdir = tmpdir)
    trim_R2.df <- attr(trim_R1.df, "reverse")

    # Merge R1 and R2 reads
    merge.df <- tryCatch({vs_fastq_mergepairs(fastq_input = trim_R1.df,
                                              reverse = trim_R2.df,
                                              minovlen = minovlen,
                                              output_format = "fasta",
                                              minlen = minlen,
                                              threads = threads,
                                              tmpdir = tmpdir)
    }, error = function(e) {
      res.df$merged_read_pairs[i] <<- 0
      res.df$R1_length[i] <<- if (nrow(trim_R1.df) > 0) round(mean(nchar(trim_R1.df$Sequence)), 2) else 0
      res.df$R2_length[i] <<- if (nrow(trim_R2.df) > 0) round(mean(nchar(trim_R2.df$Sequence)), 2) else 0
      return(NULL)
    })

    if (is.null(merge.df)) next

    # Dereplicate merged reads
    derep.df <- vs_fastx_uniques(fastx_input = merge.df,
                                 output_format = "fasta",
                                 relabel_sha1 = TRUE,
                                 tmpdir = tmpdir) |>
      dplyr::mutate(size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
      dplyr::mutate(size = as.numeric(size))

    # Find number of dereplicated merged reads with size > min_size
    derep.df_filt <- derep.df |>
      dplyr::filter(size >= min_size)

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
                                  merged_read_pairs = "Number of merged read-pairs")) +
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
