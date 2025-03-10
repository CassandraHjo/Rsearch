#' Plot read length vs. read quality
#'
#' @description
#' Generates a scatter plot visualizing the relationship between read length and
#' read quality. The y-axis can display either the mean quality score per read
#' or the expected error (EE) rate. Marginal histograms are included to show the
#' distribution of read lengths and quality metrics.
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing reads. See
#' \emph{Details}.
#' @param use_ee_rate If \code{TRUE}, the plot will display the expected error
#' rate (EE) on the y-axis instead of the mean quality score. Defaults to
#' \code{FALSE}.
#' @param plot_title If \code{TRUE} (default), a title will be displayed
#' in the plot. The title will either be "Read length vs Expected error rate
#' (EE) of read" or "Read length vs Average quality score of read", depending on
#' \code{use_ee_rate}. Set to \code{FALSE} for no title.
#'
#' @details
#' This function visualizes the relationship between read length and read
#' quality. The user can choose to plot either the
#' mean quality score per read or the expected error (EE) rate.
#'
#' \code{fastq_input} can either be a file path to a FASTQ file or a FASTQ
#' object. FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#' The EE rate is calculated as the mean of error probabilities per read, where
#' the error probability for each base is computed as \eqn{10^{(-Q/10)}} from
#' Phred scores. A lower EE rate indicates higher sequence quality, while a
#' higher EE rate suggests lower confidence in the read.
#'
#' Marginal histograms are added to display the distribution of read lengths
#' (top) and quality scores or EE rates (right).
#'
#' @return A ggplot2 object displaying the scatter plot with marginal histograms.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#'
#' # Generate and display scatter plot with mean quality score on y-axis
#' p1 <- plot_read_quality(fastq_input = fastq_input)
#' print(p1)
#'
#' # Generate and display scatter plot with mean quality score on y-axis
#' p2 <- plot_read_quality(fastq_input = fastq_input,
#'                         use_ee_rate = TRUE)
#' print(p2)
#' }
#'
#' @export
#'
plot_read_quality <- function(fastq_input,
                              use_ee_rate = FALSE,
                              plot_title = TRUE) {

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

  # Convert quality symbols to numeric scores
  fastq.tbl$Q_scores <- lapply(fastq.tbl$Quality,
                               function(Q.seq) {Q.seq |>
                                   charToRaw() |>
                                   strtoi(16L) - 33
                               })

  # Calculate mean quality score for each read
  fastq.tbl$Mean_Q_score <- sapply(fastq.tbl$Q_scores, mean)

  # Calculate expected error (EE) rate for each read
  fastq.tbl$EE_rate <- sapply(fastq.tbl$Q_scores,
                              function(Q) {
                                mean(10^(-Q/10))})

  # Add read length column
  fastq.tbl <- fastq.tbl |>
    dplyr::mutate(Length = nchar(Sequence))

  # Define color palette
  pal <- RColorBrewer::brewer.pal(5, "YlGnBu")

  # Choose y-axis variable based on user selection
  y_var <- ifelse(use_ee_rate, "EE_rate", "Mean_Q_score")
  y_label <- ifelse(use_ee_rate,
                    "Expected error rate (EE) of read",
                    "Average quality score of read")

  # Define plot title
  if (plot_title) {
    title <- paste("Read length vs", y_label)
  } else {
    title <- ""
  }

  # Plot scatter plot
  p1 <- ggplot2::ggplot(fastq.tbl,
                        ggplot2::aes(x = Length, y = .data[[y_var]])) +
    ggplot2::geom_point(alpha = 0.5, color = pal[2]) +
    ggplot2::labs(title = title,
                  x = "Read length (bases)",
                  y = y_label) +
    ggplot2::theme_minimal()

  # Add marginal histograms
  plot_with_marginal_plots <- ggExtra::ggMarginal(p1,
                                                  type = "histogram",
                                                  fill = pal[3],
                                                  col = pal[4])

  return(plot_with_marginal_plots)

}
