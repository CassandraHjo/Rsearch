#' Length statistics after merging
#'
#' @description Statistics of read lengths before and after merging.
#'
#' @param fastq_input A FASTQ file path or object containing (forward) reads.
#' @param reverse A FASTQ file path or object containing (reverse) reads.
#' @param minovlen The minimum overlap between the merged reads. Must be at least 5. Defaults to \code{10}.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details This function calculates the length of the forward reads, reverse reads,
#' the merged reads, and the overlap lengths. It uses \code{\link{vs_fastq_mergepairs}} where
#' the arguments to this function are described in detail.
#'
#' These length statistics are most typically used in order to tune the filter
#' and trimming of reads such that the merged reads are of high quality.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{length_1}: The length of the forward reads.
#'   \item \code{length_2}: The length of the reverse reads.
#'   \item \code{length_merged}: The length of the merged reads.
#'   \item \code{length_overlap}: The length of the overlap between the forward and reverse reads.
#' }
#'
#' In case of missing values for the latter two columns, it means that the
#' corresponding reads were not merged.
#'
#' The data frame has an attribute \code{plot} containing a grid plot based on the returned data frame.
#'
#' The data frame also has an attribute \code{statistics} containing the same
#' attribute returned from \code{\link{vs_fastq_mergepairs}}.
#'
#' @seealso \code{\link{vs_fastq_mergepairs}}
#'
#' @examples
#' \dontrun{
#' # Read example FASTQ files
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"), "R2_sample1_small.fq")
#'
#' # Execute merging
#' merging.tbl <- vs_merging_lengths(fastq_input = fastq_input,
#'                                   reverse = reverse,
#'                                   minovlen = 10,
#'                                   minlen = 0,
#'                                   threads = 1)
#'
#' # Extract plot
#' merging_stats_plot <- attr(merging.tbl, "plot")
#'
#' }
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_merging_lengths <- function(fastq_input,
                               reverse,
                               minovlen = 10,
                               minlen = 0,
                               threads = 1){
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

  # The merged read lengths and overlap lengths
  merged.tbl <- vs_fastq_mergepairs(R1.tbl,
                                    R2.tbl,
                                    minovlen = minovlen,
                                    minlen = minlen,
                                    threads = threads)

  # The lengths
  res.tbl <- R1.tbl |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
    dplyr::mutate(length_1 = stringr::str_length(Sequence)) |>
    dplyr::select(tag, length_1)
  res.tbl <- R2.tbl |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
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
  plot1 <- res.tbl |>
    tidyr::pivot_longer(-read_id, names_to = "type", values_to = "length") |>
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = length)) +
    ggplot2::facet_wrap(dplyr::vars(type), scales = "free")

  attr(res.tbl, "plot") <- plot1

  return(res.tbl)
}
