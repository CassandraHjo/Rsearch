#' Plot histogram of read lengths
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
#' Note that this function is only useful if you have reads of different lengths,
#' like Nanopore reads.
#'
#' @return A ggplot object displaying a histogram of read lengths.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#'
#' # Plot histogram
#' hist_plot <- plot_read_length_histogram(fastq_input = fastq_input)
#' }
#'
#' @export
#'
plot_read_length_histogram <- function(fastq_input) {

  # Handle input: file or tibble
  if (!is.character(fastq_input)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    data.tbl <- fastq_input
  } else {
    data.tbl <- microseq::readFastq(fastq_input)
  }

  # Make histogram

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  hist.plot <- data.tbl |>
    dplyr::mutate(Length = nchar(Sequence)) |>
    ggplot2::ggplot(ggplot2::aes(x = Length)) +
    ggplot2::geom_histogram(fill = pal[3], color = pal[4]) +
    ggplot2::labs(title = "Read lengths",
                  x = "Read length",
                  y = "Number of reads") +
    ggplot2::theme_minimal()

  return(hist.plot)
}
