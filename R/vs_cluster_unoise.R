#' Denoising FASTA sequences
#'
#' @description \code{vs_cluster_unoise} performs denoising of FASTA sequences from a
#' given file or object using \code{VSEARCH}´s \code{cluster_unoise} method.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object containing
#' reads to denoise. See \emph{Details}.
#' @param otutabout (Optional). A character string specifying the name of the
#' output file in an OTU table format. If \code{NULL} (default), the output is
#' returned as a tibble in R. See \emph{Details}.
#' @param minsize (Optional). Minimum abundance of cluster centroids.
#' Defaults to \code{8}.
#' @param unoise_alpha (Optional). Alpha value for the UNOISE algorithm.
#' Defaults to \code{2}.
#' @param log_file (Optional). Name of the log file to capture messages from
#' \code{VSEARCH}. If \code{NULL} (default), no log file is created.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details
#' Sequences are denoised according to the UNOISE version 3 algorithm by Robert
#' Edgar, but without the de novo chimera removal step. In this algorithm,
#' clustering of sequences depends both on their similarity and their
#' abundances. The abundance ratio (skew) is the abundance of a new
#' sequence divided by the abundance of the centroid sequence. This skew must
#' not be larger than beta if the sequences should be clustered together. Beta
#' is calculated as 2 raised to the power of minus 1 minus alpha times the
#' sequence distance. The sequence distance used is the number of mismatches in
#' the alignment, ignoring gaps. This means that the abundance must be
#' exponentially lower as the distance increases from the centroid for a new
#' sequence to be included in the cluster.
#'
#' The argument \code{minsize} will affect the total number of clusters,
#' specifying the minimum copy number required for any centroid. A larger value
#' means (in general) fewer clusters.
#'
#' \code{fasta_input} can either be a file path to a FASTA file or a FASTA
#' object. FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}, see \code{\link{readFasta}}.
#'
#' The \code{Header} column \strong{must} contain the size (copy number) for
#' each read. The size information must have the format ";size=X",
#' where X is the count for the given sequence. This is obtained by running all
#' reads through \code{\link{vs_fastx_uniques}} with \code{sizeout = TRUE}.
#'
#' You may use reads for a single sample or all reads from all samples as input.
#' In the latter case the \code{Header} must also contain sample information
#' on the format ";sample=xxx" where "xxx" is a unique sample identifier text.
#' Again, this is obtained by using \code{\link{vs_fastx_uniques}} on the reads
#' for each sample prior to this step. Use the \code{sample = "xxx"} argument,
#' where "xxx" is replaced with some unique text for each sample.
#'
#' If \code{log_file} is \code{NULL} and \code{centroids} is specified,
#' clustering statistics from \code{VSEARCH} will not be captured.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A read count table with one row for each cluster and one column for
#' each sample. If \code{otutabout} is a text it is assumed to be a file name,
#' and the results are written to this file. If no such text is supplied (default),
#' it is returned as a tibble.
#'
#' The first two columns of this tibble lists the \code{Header} and \code{Sequence}
#' of the centroid sequences for each cluster.
#'
#' The clustering statistics are included as an attribute named
#' \code{"statistics"} with the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: Total number of nucleotides used as input for
#'   clustering.
#'   \item \code{min_length_input_seq}: Length of the shortest sequence used as
#'   input for clustering.
#'   \item \code{max_length_input_seq}: Length of the longest sequence used as
#'   input for clustering.
#'   \item \code{avg_length_input_seq}: Average length of the sequences used as
#'   input for clustering.
#'   \item \code{num_clusters}: Number of clusters generated.
#'   \item \code{min_size_cluster}: Size of the smallest cluster.
#'   \item \code{max_size_cluster}: Size of the largest cluster.
#'   \item \code{avg_size_cluster}: Average size of the clusters.
#'   \item \code{num_singletons}: Number of singletons after clustering.
#'   \item \code{input}: Name of the input file/object for the clustering.
#' }
#'
#' @examples
#' \dontrun{
#' # A small fasta file
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "small.fasta")
#'
#' # Denoise sequences and read counts
#' denoised.tbl <- vs_cluster_unoise(fasta_input = fasta_input)
#' head(denoised.tbl)
#'
#' # Extract clustering statistics
#' statistics <- attr(denoised.tbl, "statistics")
#'
#' # Cluster sequences and write results to a file
#' vs_cluster_unoise(fasta_input = fasta_input,
#'                   otutabout = "otutable.tsv")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_cluster_unoise unoise denoise
#'
#' @export
#'
vs_cluster_unoise <- function(fasta_input,
                              otutabout = NULL,
                              minsize = 8,
                              unoise_alpha = 2,
                              log_file = NULL,
                              threads = 1,
                              vsearch_options = NULL,
                              tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Create empty vector for collecting temporary files
  temp_files <- character()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0 && is.character(temp_files)) {
      existing_files <- temp_files[file.exists(temp_files)]
      if (length(existing_files) > 0) {
        file.remove(existing_files)
      }
    }
  }, add = TRUE)

  # Check if FASTA input is file or tibble
  if (!is.character(fasta_input)){
    temp_file <- tempfile(pattern = "input",
                          tmpdir = tmpdir,
                          fileext = ".fa")
    temp_files <- c(temp_files, temp_file)
    microseq::writeFasta(fasta_input, temp_file)
    fasta_file <- temp_file

    # Capture original name for statistics table later
    fasta_input_name <- as.character(substitute(fasta_input))
  } else {
    fasta_file <- fasta_input

    # Capture original name for statistics table later
    fasta_input_name <- basename(fasta_input)
  }

  # Check if input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)

  # Temporary files
  centrfile <- tempfile(pattern = "centroid",
                        tmpdir = tmpdir,
                        fileext = ".fa")
  otutabfile <- tempfile(pattern = "otutab",
                         tmpdir = tmpdir,
                         fileext = ".tsv")
  temp_files <- c(temp_files, centrfile, otutabfile)

  # Build argument string for command line
  args <- c("--cluster_unoise", shQuote(fasta_file),
            "--threads", threads,
            "--minsize", minsize,
            "--unoise_alpha", unoise_alpha,
            "--centroids", centrfile,
            "--otutabout", otutabfile)

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Read results and make otu table
  centr.tbl <- microseq::readFasta(centrfile)
  if(nrow(centr.tbl) > 0){
    statistics <- calculate_cluster_statistics(centr.tbl,
                                               fasta_file,
                                               fasta_input_name)
    centr.tbl <- centr.tbl |>
      dplyr::mutate(tag = stringr::word(Header, 1, sep = ";")) |>
      dplyr::distinct(tag, .keep_all = T) |>
      dplyr::select(tag, Sequence)
    otu.tbl <- suppressMessages(readr::read_tsv(otutabfile)) |>
      dplyr::rename(tag = `#OTU ID`)
    sizes <- otu.tbl |>
      dplyr::select(-tag) |>
      as.matrix() |>
      rowSums()
    otu.tbl <- otu.tbl |>
      dplyr::left_join(centr.tbl, by = "tag") |>
      dplyr::mutate(size = sizes) |>
      dplyr::arrange(desc(size)) |>
      dplyr::mutate(Header = stringr::str_c("ZOTU_", 1:dplyr::n(), ";size=", size)) |>
      dplyr::select(-tag, -size) |>
      dplyr::relocate(Header, Sequence)
    attr(otu.tbl, "statistics") <- statistics
  } else {
    warning("No clusters found, try to lower minsize")
    otu.tbl <- NULL
  }

  # The return
  if(is.character(otutabout)){
    readr::write_delim(otu.tbl, delim = "\t", file = otutabout)
    return(invisible(NULL))
  } else {
    return(otu.tbl)
  }
}

