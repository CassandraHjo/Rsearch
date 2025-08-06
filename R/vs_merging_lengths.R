#' Length statistics after merging
#'
#' @description \code{vs_merging_lengths} computes length statistics for forward
#' reads, reverse reads, merged reads, and their overlaps before and after
#' merging.
#'
#' @param fastq_input (Required). A FASTQ file path, a FASTQ tibble (forward
#' reads), or a paired-end tibble of class \code{"pe_df"}. See \emph{Details}.
#' @param reverse (Optional). A FASTQ file path or FASTQ tibble containing
#' reverse reads. Optional if \code{fastq_input} is a \code{"pe_df"} object.
#' @param minovlen (Optional). Minimum overlap between the merged reads. Must be
#' at least 5. Defaults to \code{10}.
#' @param minlen (Optional). Minimum number of bases a sequence must have to be
#' retained. Defaults to \code{0}. See \emph{Details}.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param plot_title (Optional). If \code{TRUE} (default), a summary title will
#' be displayed in the plot. Set to \code{FALSE} for no title.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details The function uses \code{\link{vs_fastq_mergepairs}} where
#' the arguments to this function are described in detail.
#'
#' If \code{fastq_input} is an object of class \code{"pe_df"}, the reverse reads
#' are automatically extracted from its \code{"reverse"} attribute unless
#' explicitly provided via the \code{reverse} argument. This allows streamlined
#' input handling for paired-end tibbles created by
#' \code{\link{fastx_synchronize}} or \code{\link{vs_fastx_trim_filt}}.
#'
#' These length statistics are most typically used in order to tune the filter
#' and trimming of reads such that the merged reads are of high quality.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{length_1}: The length of the forward reads.
#'   \item \code{length_2}: The length of the reverse reads.
#'   \item \code{length_merged}: The length of the merged reads.
#'   \item \code{length_overlap}: The length of the overlap between the forward
#'   and reverse reads.
#' }
#'
#' In case of missing values for the latter two columns, it means that the
#' corresponding reads were not merged.
#'
#' The tibble includes additional attributes:
#' \describe{
#'   \item{\code{plot}}{A \code{\link[ggplot2]{ggplot2}} object visualizing the
#'   returned data frame.}
#'   \item{\code{statistics}}{Additional statistics returned from
#'   \code{\link{vs_fastq_mergepairs}}.}
#' }
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R1.fq")
#' R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R2.fq")
#'
#' # Run function
#' merging.tbl <- vs_merging_lengths(fastq_input = R1.file,
#'                                   reverse = R2.file)
#'
#' # Display plot
#' merging_stats_plot <- attr(merging.tbl, "plot")
#' print(merging_stats_plot)
#'
#' }
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_merging_lengths merging_lengths
#'
#' @export
#'
vs_merging_lengths <- function(fastq_input,
                               reverse = NULL,
                               minovlen = 10,
                               minlen = 0,
                               threads = 1,
                               plot_title = TRUE,
                               tmpdir = NULL) {
  # The forward reads
  if (!is.character(fastq_input)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    R1.tbl <- fastq_input
  } else {
    R1.tbl <- microseq::readFastq(fastq_input)
  }

  # Check for pe_df and extract reverse if needed
  if (is_pe_df(fastq_input) && is.null(reverse)) {
    reverse <- attr(fastq_input, "reverse")
    if (is.null(reverse)) {
      stop("fastq_input has class 'pe_df' but no 'reverse' attribute found.")
    }
  }

  # Read reverse reads
  if (is.null(reverse)) {
    stop("No reverse reads provided. Please supply reverse or use a 'pe_df' object.")
  }

  # The reverse reads
  if (!is.character(reverse)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(reverse))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    R2.tbl <- reverse
  } else {
    R2.tbl <- microseq::readFastq(reverse)
  }

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # The merged read lengths and overlap lengths
  merged.tbl <- vs_fastq_mergepairs(R1.tbl,
                                    R2.tbl,
                                    minovlen = minovlen,
                                    minlen = minlen,
                                    threads = threads,
                                    tmpdir = tmpdir)

  # The lengths
  res.tbl <- R1.tbl |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
    dplyr::mutate(length_1 = stringr::str_length(Sequence)) |>
    dplyr::select(tag, length_1)
  res.tbl <- R2.tbl |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/2$")) |>
    dplyr::mutate(length_2 = stringr::str_length(Sequence)) |>
    dplyr::select(tag, length_2) |>
    dplyr::full_join(res.tbl, by = "tag")
  res.tbl <- merged.tbl |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
    dplyr::mutate(length_merged = stringr::str_length(Sequence)) |>
    dplyr::select(tag, length_merged) |>
    dplyr::full_join(res.tbl, by = "tag") |>
    dplyr::mutate(length_overlap = length_1 + length_2 - length_merged) |>
    dplyr::relocate(read_id = tag, length_1, length_2, length_merged, length_overlap)

  attr(res.tbl, "statistics") <- attr(merged.tbl, "statistics")

  # Plotting

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  # Check if Length R1 or R2 has only one unique value
  unique_length_1 <- unique(res.tbl$length_1)
  unique_length_2 <- unique(res.tbl$length_2)

  plot_r1 <- if(length(unique_length_1) == 1) {
    res.tbl |>
      dplyr::filter(!is.na(length_1)) |>
      ggplot2::ggplot(ggplot2::aes(x = as.factor(length_1))) +
      ggplot2::geom_bar(fill = pal[3], color = pal[4], width = 0.2) +
      ggplot2::labs(title = "Length R1", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(face = "bold"))
  } else {
    res.tbl |>
      dplyr::filter(!is.na(length_1)) |>
      ggplot2::ggplot(ggplot2::aes(x = length_1)) +
      ggplot2::geom_histogram(binwidth = 1, fill = pal[3], color = pal[4]) +
      ggplot2::scale_x_continuous(limits = c(min(res.tbl$length_1) - 5, max(res.tbl$length_1) + 5)) +
      ggplot2::labs(title = "Length R1", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(face = "bold"))
  }

  plot_r2 <- if(length(unique_length_2) == 1) {
    res.tbl |>
      dplyr::filter(!is.na(length_2)) |>
      ggplot2::ggplot(ggplot2::aes(x = as.factor(length_2))) +
      ggplot2::geom_bar(fill = pal[3], color = pal[4], width = 0.2) +
      ggplot2::labs(title = "Length R2", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(face = "bold"))
  } else {
    res.tbl |>
      dplyr::filter(!is.na(length_2)) |>
      ggplot2::ggplot(ggplot2::aes(x = length_2)) +
      ggplot2::geom_histogram(binwidth = 1, fill = pal[3], color = pal[4]) +
      ggplot2::scale_x_continuous(limits = c(min(res.tbl$length_2) - 5, max(res.tbl$length_2) + 5)) +
      ggplot2::labs(title = "Length R2", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(face = "bold"))
  }

  # Create separate plots for merged reads and overlap
  p3 <- ggplot2::ggplot(dplyr::filter(res.tbl, !is.na(length_merged)),
                        ggplot2::aes(x = length_merged)) +
    ggplot2::geom_histogram(binwidth = 5, fill = pal[3], color = pal[4], na.rm = TRUE) +
    ggplot2::labs(title = "Length of merged reads", x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(face = "bold"))

  p4 <-  ggplot2::ggplot(dplyr::filter(res.tbl, !is.na(length_overlap)),
                         ggplot2::aes(x = length_overlap)) +
    ggplot2::geom_histogram(binwidth = 5, fill = pal[3], color = pal[4], na.rm = TRUE) +
    ggplot2::labs(title = "Length of overlap", x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(face = "bold"))

  # Arrange plots in a grid
  combined_plot <- cowplot::plot_grid(plot_r1, plot_r2, p3, p4, ncol = 2)

  # Define plot title
  if (plot_title) {
    title <- paste0("Merged ",
                    sum(!is.na(res.tbl$length_merged)),
                    " read pairs out of ",
                    nrow(res.tbl),
                    " (",
                    round(100 * sum(!is.na(res.tbl$length_merged)) / nrow(res.tbl)), "%)")
  } else {
    title <- ""
  }

  # "Draw" the plot title
  common_title <- cowplot::ggdraw() +
    cowplot::draw_label(title,
                        size = 14,
                        x = 0.01,
                        hjust = 0,
                        fontface = "bold")

  # "Draw" common x-axis label
  common_x <- cowplot::ggdraw() +
    cowplot::draw_label("Length (bases)",
                        size = 14,
                        x = 0.5,
                        hjust = 0.5,
                        fontface = "bold")

  # "Draw" common y-axis label
  common_y <- cowplot::ggdraw() +
    cowplot::draw_label("Number of reads",
                        size = 14,
                        angle = 90,
                        y = 0.5,
                        vjust = 0.5,
                        fontface = "bold")

  # Combine title, main plot and common x-axis label
  final_plot_no_y <- cowplot::plot_grid(common_title, combined_plot, common_x,
                                        ncol = 1, rel_heights = c(0.1, 1, 0.1))

  # Add y-axis title
  final_plot <- cowplot::plot_grid(common_y, final_plot_no_y, ncol = 2, rel_widths = c(0.1, 1))

  attr(res.tbl, "plot") <- final_plot

  return(res.tbl)
}
