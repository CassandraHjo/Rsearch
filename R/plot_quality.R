#' Plot median quality scores per position for FASTQ reads
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing (forward)
#' reads. See \emph{Details}.
#' @param reverse An optional FASTQ file path or FASTQ tibble containing reverse
#' reads. Defaults to \code{NULL}. See \emph{Details}.
#' @param quantile_lower The lower quantile threshold for the error bars in the
#' plot. Defaults to \code{0.25}.
#' @param quantile_upper The upper quantile threshold for the error bars in the
#' plot. Defaults to \code{0.75}.
#'
#'
#' @details
#' Median quality scores for each position in the reads in the input files
#' (\code{fastq_input} and \code{reverse}) are plotted with ggplot2, displaying
#' the median quality score along with error bars indicating the selected
#' quantile range.
#'
#' \code{fastq_input} and \code{reverse} can either be file paths to FASTQ files
#' or FASTQ objects. FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' If \code{reverse} is provided, it is plotted together with the first plot in
#' its own panel. Note that the x-axis in this plot is reversed.
#'
#' The default error bars represent the interquartile range (25%-75%). Custom
#' quantile ranges can be specified via \code{quantile_lower} and
#' \code{quantile_upper}.
#'
#' @return A ggplot object with the quality plot(s).
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "R2_sample1_small.fq")
#'
#' # Plot
#' qual_plots <- plot_quality(fastq_input = fastq_input,
#'                            reverse = reverse)
#'
#' # Print plot
#' qual_plots
#' }
#'
#' @export
#'
#' @importFrom stats quantile median
#' @import patchwork
#'
plot_quality <- function(fastq_input,
                         reverse = NULL,
                         quantile_lower = 0.25,
                         quantile_upper = 0.75) {

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
  quantiles_quality_R1 <- apply(fastq_qual_matrix,
                                2,
                                quantile,
                                probs = c(quantile_lower, quantile_upper),
                                na.rm = TRUE)

  # Make plotting data frame
  df_R1 <- data.frame(
    Position = 1:ncol(fastq_qual_matrix),
    MedianQuality = median_quality_R1,
    Lower = quantiles_quality_R1[1, ],
    Upper = quantiles_quality_R1[2, ])

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  R1.plot <- ggplot2::ggplot(df_R1, ggplot2::aes(x = Position, y = MedianQuality)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                           width = 0.2, color = pal[2]) +
    ggplot2::geom_line(color = pal[4]) +
    ggplot2::geom_point(color = pal[3]) +
    ggplot2::labs(title = "Median quality score in each position",
                  x = "Base position",
                  y = "Quality score") +
    ggplot2::theme_minimal()

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
    quantiles_quality_R2 <- apply(reverse_qual_matrix,
                                  2,
                                  quantile,
                                  probs = c(quantile_lower, quantile_upper),
                                  na.rm = TRUE)

    # Make plotting data frame
    df_R2 <- data.frame(
      Position = 1:ncol(reverse_qual_matrix),
      MedianQuality = median_quality_R2,
      Lower = quantiles_quality_R2[1, ],
      Upper = quantiles_quality_R2[2, ])

    # Define color palette
    pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

    R2.plot <- ggplot2::ggplot(df_R2,
                               ggplot2::aes(x = Position, y = MedianQuality)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                             width = 0.2, color = pal[2]) +
      ggplot2::geom_line(color = pal[4]) +
      ggplot2::geom_point(color = pal[3]) +
      ggplot2::scale_x_reverse() +
      ggplot2::labs(title = "R2 reads",
                    x = "Base position",
                    y = "Quality score") +
      ggplot2::theme_minimal()

    # Add new title to R1 plot
    R1.plot <- R1.plot +
      ggplot2::ggtitle("R1 reads")

    # Combine using patchwork
    combined_plot <- R1.plot + R2.plot +
      patchwork::plot_annotation(title = "Median quality score in each position")

    return(combined_plot)
  }

  return(R1.plot)
}
