#' Display quality scores per position for FASTQ reads
#'
#' @description
#' Generates a plot displaying the quality scores for each
#' position in FASTQ reads.
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing (forward)
#' reads. See \emph{Details}.
#' @param reverse An optional FASTQ file path or FASTQ tibble containing reverse
#' reads. Defaults to \code{NULL}. See \emph{Details}.
#' @param quantile_lower The lower quantile threshold for the error bars in the
#' plot. Defaults to \code{0.25}.
#' @param quantile_upper The upper quantile threshold for the error bars in the
#' plot. Defaults to \code{0.75}.
#' @param plot_title The title of the plot. Defaults to
#' \code{"Per-position quality scores: median and mean"}. Set to \code{""} for
#' no title.
#' @param show_median If \code{TRUE} (default), a line representing the median
#' quality scores is added to the plot.
#' @param show_mean If \code{TRUE} (default), a line representing the mean
#' quality scores is added to the plot.
#'
#' @details
#' The mean and median quality scores for each base position over all reads in
#' the input file(s) (\code{fastq_input} (and \code{reverse}) are plotted as
#' curves. The vertical bars at each base indicate the interquartile range.
#'
#' \code{fastq_input} and \code{reverse} can either be file paths to FASTQ files
#' or FASTQ objects. FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}, see
#' \code{\link{readFastq}}.
#'
#' If \code{reverse} is provided, it is plotted together with the first plot in
#' its own panel. Note that the x-axis in this plot is reversed.
#'
#' The default error bars represent the interquartile range (25%-75%) in the
#' quality scores. Custom quantile ranges can be specified via
#' \code{quantile_lower} and \code{quantile_upper}. Additionally, the median and
#' mean quality lines may be turned off by
#' setting \code{show_median = FALSE} or \code{show_mean = FALSE}, respectively.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # Define inputs
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "small_R1.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R2.fq")
#'
#' # Generate and display quality plot with both median and mean lines
#' qual_plots <- plot_base_quality(fastq_input = fastq_input,
#'                                 reverse = reverse)
#' print(qual_plots)
#'
#' # Generate and display quality plot without the plot title
#' qual_plots_wo_title <- plot_base_quality(fastq_input = fastq_input,
#'                                          reverse = reverse,
#'                                          plot_title = "")
#' print(qual_plots_wo_title)
#'
#' # Generate a plot showing only the median quality line
#' qual_plots_median_only <- plot_base_quality(fastq_input = fastq_input,
#'                                             reverse = reverse,
#'                                             show_mean = FALSE)
#' print(qual_plots_median_only)
#'
#' # Generate a plot showing only the mean quality line
#' qual_plots_mean_only <- plot_base_quality(fastq_input = fastq_input,
#'                                           reverse = reverse,
#'                                           show_median = FALSE)
#' print(qual_plots_mean_only)
#' }
#'
#' @export
#'
#' @importFrom stats quantile median
#'
plot_base_quality <- function(fastq_input,
                              reverse = NULL,
                              quantile_lower = 0.25,
                              quantile_upper = 0.75,
                              plot_title = "Per-position quality scores: median and mean",
                              show_median = TRUE,
                              show_mean = TRUE) {

  # Handle input: file or tibble
  if (!is.character(fastq_input)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    fastq.tbl <- fastq_input
  } else {
    fastq.tbl <- microseq::readFastq(fastq_input)
  }

  # Handle reverse: file or tibble
  if (!is.null(reverse)){
    if (!is.character(reverse)){
      # Ensure required columns exist
      required_cols_rev <- c("Header", "Sequence", "Quality")
      if (!all(required_cols_rev %in% colnames(reverse))) {
        stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
      }
      reverse.tbl <- reverse
    } else {
      reverse.tbl <- microseq::readFastq(reverse)
    }
  }

  # Validate quantile range
  if (quantile_lower < 0 | quantile_lower > 1 |
      quantile_upper < 0 | quantile_upper > 1) {
    stop("Invalid values for quantile range. Choose values between 0 and 1.")
  }

  if (quantile_lower >= quantile_upper) {
    stop("Invalid quantile range: 'quantile_lower' must be smaller than 'quantile_upper'.")
  }

  # Make fastq.tbl plot

  # Convert quality symbols to numeric scores
  fastq.tbl$Q_scores <- lapply(fastq.tbl$Quality,
                               function(Q.seq) {Q.seq |>
                                   charToRaw() |>
                                   strtoi(16L) - 33
                               })

  # Find max length among all reads
  max_length <- max(sapply(fastq.tbl$Q_scores, length))

  # Pad shorter reads with NAs
  fastq.tbl$Q_scores <- lapply(fastq.tbl$Q_scores,
                               function(q) {
                                 c(q, rep(NA, max_length - length(q)))
                               })
  # Make quality score matrix
  fastq_qual_matrix <- do.call(rbind, fastq.tbl$Q_scores)

  # Calculate statistics
  median_quality_R1 <- apply(fastq_qual_matrix, 2, median, na.rm = TRUE)
  mean_quality_R1 <- apply(fastq_qual_matrix, 2, mean, na.rm = TRUE)
  quantiles_quality_R1 <- apply(fastq_qual_matrix,
                                2,
                                quantile,
                                probs = c(quantile_lower, quantile_upper),
                                na.rm = TRUE)

  # Make plotting data frame
  df_R1 <- data.frame(
    Position = 1:ncol(fastq_qual_matrix),
    MedianQuality = median_quality_R1,
    MeanQuality = mean_quality_R1,
    Lower = quantiles_quality_R1[1, ],
    Upper = quantiles_quality_R1[2, ])

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  # Plot error bars and labels
  R1.plot <- ggplot2::ggplot(df_R1, ggplot2::aes(x = Position)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                           width = 0.2, color = pal[2]) +
    ggplot2::labs(title = plot_title,
                  x = "Base position",
                  y = "Quality score",
                  color = "Statistic") +
    ggplot2::theme_minimal()

  # Create empty vector for storing color mapping
  color_mapping <- c()

  # Plot median line if specified
  if (show_median) {
    R1.plot <- R1.plot + ggplot2::geom_line(ggplot2::aes(y = MedianQuality, color = "Median"))
    color_mapping["Median"] <- pal[4]
  }

  # Plot mean line if specified
  if (show_mean) {
    R1.plot <- R1.plot + ggplot2::geom_line(ggplot2::aes(y = MeanQuality, color = "Mean"))
    color_mapping["Mean"] <- pal[3]
  }

  # Add correct colors
  if (length(color_mapping) > 0) {
    R1.plot <- R1.plot + ggplot2::scale_color_manual(values = color_mapping)
  }

  # Make reverse.tbl plot
  if (!is.null(reverse)){

    # Convert quality symbols to numeric scores
    reverse.tbl$Q_scores <- lapply(reverse.tbl$Quality,
                                   function(Q.seq) {Q.seq |>
                                       charToRaw() |>
                                       strtoi(16L) - 33
                                   })

    # Find max length among all reads
    max_length <- max(sapply(reverse.tbl$Q_scores, length))

    # Pad shorter reads with NAs
    reverse.tbl$Q_scores <- lapply(reverse.tbl$Q_scores,
                                   function(q) {
                                     c(q, rep(NA, max_length - length(q)))
                                   })

    # Make quality score matrix
    reverse_qual_matrix <- do.call(rbind, reverse.tbl$Q_scores)

    # Calculate statistics
    median_quality_R2 <- apply(reverse_qual_matrix, 2, median, na.rm = TRUE)
    mean_quality_R2 <- apply(reverse_qual_matrix, 2, mean, na.rm = TRUE)
    quantiles_quality_R2 <- apply(reverse_qual_matrix,
                                  2,
                                  quantile,
                                  probs = c(quantile_lower, quantile_upper),
                                  na.rm = TRUE)

    # Make plotting data frame
    df_R2 <- data.frame(
      Position = 1:ncol(reverse_qual_matrix),
      MedianQuality = median_quality_R2,
      MeanQuality = mean_quality_R2,
      Lower = quantiles_quality_R2[1, ],
      Upper = quantiles_quality_R2[2, ])

    # Plot error bars and labels
    R2.plot <- ggplot2::ggplot(df_R2, ggplot2::aes(x = Position)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                             width = 0.2, color = pal[2]) +
      ggplot2::labs(title = "R2 reads",
                    x = "Base position",
                    y = "Quality score",
                    color = "Statistic") +
      ggplot2::theme_minimal()

    # Create empty vector for storing color mapping
    color_mapping <- c()

    # Plot median line if specified
    if (show_median) {
      R2.plot <- R2.plot + ggplot2::geom_line(ggplot2::aes(y = MedianQuality, color = "Median"))
      color_mapping["Median"] <- pal[4]
    }

    # Plot mean line if specified
    if (show_mean) {
      R2.plot <- R2.plot + ggplot2::geom_line(ggplot2::aes(y = MeanQuality, color = "Mean"))
      color_mapping["Mean"] <- pal[3]
    }

    # Add correct colors
    if (length(color_mapping) > 0) {
      R2.plot <- R2.plot + ggplot2::scale_color_manual(values = color_mapping)
    }

    # Add new title to R1 plot
    R1.plot <- R1.plot +
      ggplot2::ggtitle("R1 reads")

    # Extract legend from R1 plot for grid plotting
    grobs <- ggplot2::ggplotGrob(R1.plot)$grobs
    legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

    # Remove legends from R1 and R2 plot
    R1.plot <- R1.plot + ggplot2::theme(legend.position = "none")
    R2.plot <- R2.plot + ggplot2::theme(legend.position = "none",
                                        axis.title.y = ggplot2::element_blank())

    # Create common title
    common_title <- cowplot::ggdraw() +
      cowplot::draw_label(plot_title, size = 14, x = 0.01, hjust = 0)

    # Combine the two plots
    combined_plot <- cowplot::plot_grid(R1.plot, R2.plot, ncol = 2)

    # Combine plots and title
    plot_with_title <- cowplot::plot_grid(common_title,
                                          combined_plot,
                                          ncol = 1,
                                          rel_heights = c(0.1, 1))

    # Combine plot with legend
    final_plot <- cowplot::plot_grid(plot_with_title,
                                     legend,
                                     rel_widths = c(1, 0.15))

    return(final_plot)
  }
  return(R1.plot)
}
