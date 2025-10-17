#' Plot distribution of expected error (EE) rate of reads
#'
#' @description
#' Generates a histogram visualizing the distribution of the expected error (EE)
#' rate for reads. The EE rate represents the cumulative probability of errors
#' in a read, calculated from Phred quality scores.
#'
#' @param fastq_input (Required). A FASTQ file path or FASTQ object containing
#' reads. See \emph{Details}.
#' @param n_bins (Optional). Number of bins used in the histogram. Defaults to
#' \code{30}, which is the default value in \code{ggplot2::geom_histogram()}.
#' @param plot_title (Optional). The title of the plot. Defaults to
#' \code{"Distribution of the expected error (EE) rate of reads"}. Set to
#' \code{""} for no title.
#'
#' @details
#' A histogram is plotted using ggplot2 to visualize the distribution of EE
#' rates. The user can adjust the number of bins in the histogram using the
#' \code{n_bins} parameter.
#'
#' \code{fastq_input} can either be a file path to a FASTQ file or a FASTQ
#' object. FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}, see \code{\link[microseq]{readFastq}}.
#'
#' The EE rate is calculated as the sum of error probabilities per read, where
#' the error probability for each base is computed as \eqn{10^{(-Q/10)}} from
#' Phred scores. A lower EE rate indicates higher sequence quality, while a
#' higher EE rate suggests lower confidence in the read.
#'
#' If \code{fastq_input} contains more than 10 000 reads, the function will
#' randomly select 10 000 rows for downstream calculations. This subsampling is
#' performed to reduce computation time and improve performance on large
#' datasets.
#'
#' @return A ggplot2 object displaying the histogram of EE rate distribution.
#'
#' @examples
#' # Define input file path
#' fastq_input <- system.file("extdata/small_R1.fq", package = "Rsearch")
#'
#' # Generate and display histogram
#' ee_plot <- plot_ee_rate_dist(fastq_input = fastq_input)
#' print(ee_plot)
#'
#' @export
#'
plot_ee_rate_dist <- function(fastq_input,
                              n_bins = 30,
                              plot_title = "Distribution of the expected error (EE) rate of reads") {

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

  # If it is more than 10 000 reads, take a random sample of 10 000 reads
  if (nrow(fastq.tbl) > 10000) {
    sample_indices <- sample(seq_len(nrow(fastq.tbl)), 10000)
    fastq.tbl <- fastq.tbl[sample_indices, ]
  }

  # Convert quality symbols to numeric scores
  fastq.tbl$Q_scores <- lapply(fastq.tbl$Quality,
                               function(Q.seq) {Q.seq |>
                                   charToRaw() |>
                                   strtoi(16L) - 33
                               })

  # Calculate expected error (EE) rate for each read
  fastq.tbl$EE_rate <- sapply(fastq.tbl$Q_scores,
                              function(Q) {
                                mean(10^(-Q/10))})

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  # Create histogram
  ee_plot <- ggplot2::ggplot(fastq.tbl,
                             ggplot2::aes(x = EE_rate)) +
    ggplot2::geom_histogram(bins = n_bins,
                            fill = pal[3],
                            color = pal[4],
                            boundary = 0) +
    ggplot2::labs(title = plot_title,
                  x = "EE-rate",
                  y = "Number of reads") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(face = "bold"))

  return(ee_plot)
}
