#' Plot read length vs. read quality in heatmap
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing reads. See
#' \emph{Details}.
#'
#' @details
#' A histogram with read lengths is plotted with ggplot2, displaying
#' the number of reads with different lengths.
#'
#' \code{fastq_input} can either be a file path to a FASTQ file or a FASTQ
#' object. FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#' Note that this function is most useful if you have reads of different
#' lengths, like Nanopore reads.
#'
#' @return A ggplot object displaying a histogram of read lengths.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#'
#' # Plot heatmap
#' heatmap <- plot_heatmap(fastq_input = fastq_input)
#' }
#'
#' @export
#'
plot_heatmap <- function(fastq_input) {

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

  fastq.tbl <- fastq.tbl |>
    dplyr::mutate(Length = nchar(Sequence))

  # Define color palette
  pal <- RColorBrewer::brewer.pal(5, "YlGnBu")

  heatmap <- ggplot2::ggplot(fastq.tbl,
                             ggplot2::aes(x = Length, y = Mean_Q_score)) +
    ggplot2::geom_bin_2d(binwidth = c(10, 1)) +
    ggplot2::scale_fill_gradient(low = pal[2],
                                 high = pal[5],
                                 name = "Number of reads") +
    ggplot2::labs(title = "Read length vs average read quality",
                  x = "Read length",
                  y = "Average quality score") +
    ggplot2::theme_minimal()

  return(heatmap)

}
