#' Plot distribution of size values
#'
#' @description
#' Generates a bar plot representing the distribution of size values from a
#' FASTA or FASTQ file/object. The y-axis is log10-transformed to handle
#' variations in read count across different size values.
#'
#' @param fastx_input A FASTA/FASTQ file path or FASTA/FASTQ object containing
#' reads with size values embedded in the \code{Header} column. See
#' \emph{Details}.
#' @param input_format The format of the input file. Must be \code{"fasta"} or
#' \code{"fastq"} if \code{fastx_input} is a file path. Defaults to \code{NULL}.
#' @param cutoff An optional numeric value specifying a size threshold. Reads
#' with size greater than this value will be grouped into a single category
#' labeled \code{"> cutoff"} in the plot. Defaults to \code{NULL} (no cutoff
#' applied).
#' @param y_breaks A numeric vector specifying the breakpoints for the y-axis.
#' Defaults to \code{c(1, 10, 100, 1000, 2000, 3000, 10000)}.
#'
#' @details
#'
#' \code{fastx_input} can either be a file path to FASTA/FASTQ file or a
#' FASTA/FASTQ object. FASTA objects are tibbles that contain the
#' columns \code{Header} and \code{Sequence}. FASTQ objects are tibbles that
#' contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#' The \code{Header} column must contain the size values for each read.
#'
#' The \code{Header} column must contain size annotations formatted as
#' \code{;size=<int>}.
#'
#' The y-axis of the plot is log10-transformed to handle variations in read
#' counts across different size values. The breakpoints of the y-axis can be
#' modified using the \code{y_breaks} parameter.
#'
#' @return A ggplot2 object displaying a bar plot of size distribution.
#'
#' @examples
#' \dontrun{
#' # Define input file
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "derep_R1_sample1_small.fa")
#'
#' # Generate and display plot without cutoff
#' size_plot <- plot_size_dist(fastx_input = fastx_input,
#'                             input_format = "fasta")
#' print(size_plot)
#'
#' # Generate and display plot with a cutoff at size 100
#' size_plot <- plot_size_dist(fastx_input = fastx_input,
#'                             input_format = "fasta",
#'                             cutoff = 100)
#' print(size_plot)
#'
#' # Generate and display plot with custom y-axis breaks
#' size_plot <- plot_size_dist(fastx_input = fastx_input,
#'                             input_format = "fasta",
#'                             y_breaks = c(1, 50, 500, 5000))
#' print(size_plot)
#' }
#'
#' @export
#'
plot_size_dist <- function(fastx_input,
                           input_format = NULL,
                           cutoff = NULL,
                           y_breaks = c(1, 10, 100, 1000, 2000, 3000, 10000)) {


  # Handle input if tibble is provided
  if (!is.character(fastx_input)){ # If tibble
    required_cols <- c("Header", "Sequence")
    if (!all(required_cols %in% colnames(fastx_input))) {
      stop("FASTX object must contain columns: Header and Sequence")
    }
    fastx.tbl <- fastx_input
  } else {
    # Handle input if file path is provided
    if (!file.exists(fastx_input)) {
      stop("Cannot find input file: ", fastx_input)
    }

    if (is.null(input_format) || !(input_format %in% c("fasta", "fastq"))) {
      stop("Input format must be specified as 'fasta' or 'fastq' if input is a file path.")
    }

    fastx.tbl <- if (input_format == "fasta") {
      microseq::readFasta(fastx_input)
    } else {
      microseq::readFastq(fastx_input) |>
        dplyr::select(-Quality)
    }
  }

  # Extract size value from Header and clean Header
  fastx.tbl <- fastx.tbl |>
    dplyr::mutate(size = stringr::str_remove(Header, ".+;size=")) |>
    dplyr::mutate(size = as.integer(size)) |>
    dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

  # Apply cutoff if specified: values above cutoff become "> cutoff"
  if (!is.null(cutoff)) {
    fastx.tbl <- fastx.tbl |>
      dplyr::mutate(size = ifelse(size > cutoff, paste0("> ", cutoff), as.character(size)))
  }

  # Group by size and count reads
  size_dist.tbl <- fastx.tbl |>
    dplyr::group_by(size) |>
    dplyr::summarize(num_reads = dplyr::n()) |>
    dplyr::ungroup()

  # Convert size to a factor for correct ordering in the plot
  size_dist.tbl <- size_dist.tbl |>
    dplyr::mutate(size = factor(size,
                                levels = c(sort(as.numeric(unique(size_dist.tbl$size[size_dist.tbl$size != paste0("> ", cutoff)]))),
                                           paste0("> ",
                                                  cutoff
                                           )
                                )
    )
    )

  # Define color palette
  pal <- RColorBrewer::brewer.pal(4, "YlGnBu")

  # Create bar plot
  size_plot <- ggplot2::ggplot(size_dist.tbl,
                               ggplot2::aes(x = size, y = num_reads)) +
    ggplot2::geom_bar(stat = "identity", fill = pal[3], color = pal[4]) +
    ggplot2::labs(title = "Size distribution",
                  x = "Size",
                  y = "Number of reads") +
    ggplot2::scale_y_log10(breaks = y_breaks) +
    ggplot2::theme_minimal()

  return(size_plot)
}
