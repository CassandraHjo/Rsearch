#' Plot average quality scores per position for FASTQ reads
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing (forward)
#' reads. See \emph{Details}.
#' @param reverse An optional FASTQ file path or FASTQ tibble containing reverse
#' reads. Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Average quality scores for each position in the reads in the input files
#' (\code{fastq_input} and \code{reverse}) are plotted with ggplot2, displaying
#' average quality scores with error bars indicating ± one standard deviation.
#'
#' \code{fastq_input} and \code{reverse} can either be file paths to FASTQ files
#' or FASTQ objects. FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' If \code{reverse} is provided, it is plotted together with the first plot in
#' its own panel. Note that the x-axis in this plot is reversed.
#'
#' Note that this function assumes that all FASTQ reads are of equal length. Quality
#' scores are combined into a matrix where each row is a read and each column is
#' a base position.
#'
#' @return A ggplot object for the forward reads. If reverse reads are provided,
#' the reverse plot is attached as an attribute (\code{"reverse"}) to the
#' returned object.
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
#' # Extract plots
#' R1_plot <- qual_plots
#' R2_plot <- attr(qual_plots, "reverse")
#' }
#'
#' @export
#'
#' @importFrom stats sd
#'
plot_quality <- function(fastq_input,
                         reverse = NULL) {

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

  # Make fastq.tbl plot

  # Convert quality symbols to numeric scores
  fastq.tbl$Q_scores <- lapply(fastq.tbl$Quality,
                               function(Q.seq) {Q.seq |>
                                   charToRaw() |>
                                   strtoi(16L) - 33
                               })
  # Make quality score matrix
  fastq_qual_matrix <- do.call(rbind, fastq.tbl$Q_scores)

  # Calculate statistics
  mean_quality_R1 <- colMeans(fastq_qual_matrix)
  sd_quality_R1 <- apply(fastq_qual_matrix, 2, stats::sd)

  # Make plotting data frame
  df_R1 <- data.frame(
    Position = 1:ncol(fastq_qual_matrix),
    MeanQuality = mean_quality_R1,
    Lower = mean_quality_R1 - sd_quality_R1,
    Upper = mean_quality_R1 + sd_quality_R1)

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  R1.plot <- ggplot2::ggplot(df_R1, ggplot2::aes(x = Position, y = MeanQuality)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                           width = 0.2, color = pal[2]) +
    ggplot2::geom_line(color = pal[4]) +
    ggplot2::geom_point(color = pal[3]) +
    ggplot2::labs(title = "Average quality score in each position",
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
    # Make quality score matrix
    reverse_qual_matrix <- do.call(rbind, reverse.tbl$Q_scores)

    # Calculate statistics
    mean_quality_R2 <- colMeans(reverse_qual_matrix)
    sd_quality_R2 <- apply(reverse_qual_matrix, 2, sd)

    # Make plotting data frame
    df_R2 <- data.frame(
      Position = 1:ncol(reverse_qual_matrix),
      MeanQuality = mean_quality_R2,
      Lower = mean_quality_R2 - sd_quality_R2,
      Upper = mean_quality_R2 + sd_quality_R2)

    # Define color palette
    pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

    R2.plot <- ggplot2::ggplot(df_R2,
                               ggplot2::aes(x = Position, y = MeanQuality)) +
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

    combined_plot <- gridExtra::grid.arrange(R1.plot,
                                             R2.plot,
                                             ncol = 2,
                                             top = grid::textGrob("Average quality score in each position",
                                                                  x = grid::unit(0, "npc"),
                                                                  just = "left",
                                                                  gp = grid::gpar(fontsize = 18))
    )
    return(invisible(combined_plot))
  }

  return(R1.plot)
}
