#' Denoising FASTA sequences
#'
#' @description \code{vs_cluster_unoise} performs denoising of FASTA sequences
#' from a given file or object using \code{VSEARCH}´s \code{cluster_unoise}
#' method.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object containing
#' reads to denoise. See \emph{Details}.
#' @param otutabout (Optional). A character string specifying the name of the
#' output file in an OTU table format. If \code{NULL} (default), the output is
#' returned as a tibble in R. See \emph{Details}.
#' @param clusters (Optional). Output each cluster to a separate FASTA file
#' using the prefix string provided and a ticker (0, 1, 2, etc.) to construct
#' the path and the filenames. Defaults to \code{NULL}, meaning no such files
#' are created.
#' @param uc (Optional). A character string specifying the name of the output
#' file in a uclust-like format. If \code{NULL} (default), no output is written
#' to a file. If \code{TRUE}, the output is returned as a tibble. See
#' \emph{Details}.
#' @param minsize (Optional). Minimum abundance of cluster centroids.
#' Defaults to \code{8}.
#' @param unoise_alpha (Optional). Alpha value for the UNOISE algorithm.
#' Defaults to \code{2}.
#' @param relabel (Optional). Relabel sequences using the given prefix and a
#' ticker to construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 (Optional). If \code{TRUE} (default), relabel sequences
#' using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
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
#' \code{Sequence}, see \code{\link[microseq]{readFasta}}.
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
#' \code{otutabout} gives the option to output the results in an OTU
#' table format with tab-separated columns. When writing to a file, the first
#' line starts with the string "#OTU ID", followed by a tab-separated list of
#' all sample identifiers (formatted as "sample=X"). Each subsequent line,
#' corresponding to an OTU, begins with the OTU identifier and is followed by
#' tab-separated abundances for that OTU in each sample. If \code{otutabout} is
#' a character string, the output is written to the specified file. If
#' \code{otutabout} is \code{TRUE}, the function returns the OTU table as a
#' tibble, where the first column is named \code{otu_id} instead of "#OTU ID".
#'
#' \code{uc} gives the option to output the results in a uclust-like format.
#' This is a tab-separated uclust-like format with 10 columns and 3 different
#' types of entries (S, H or C). Each FASTA sequence i the input can either be a
#' cluster centroid (S) or a hit (H) assigned to a cluster. Cluster records (C)
#' summarize information (size, centroid label) for each cluster. See the '--uc'
#' section in the \code{VSEARCH} manual for more information.
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
#' and the results are written to this file. If no such text is supplied
#' (default), it is returned as a tibble.
#'
#' The first two columns of this tibble lists the \code{Header} and
#' \code{Sequence} of the centroid sequences for each cluster.
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
#' @aliases vs_cluster_unoise cluster_unoise unoise denoise
#'
#' @export
#'
vs_cluster_unoise <- function(fasta_input,
                              otutabout = NULL,
                              uc = NULL,
                              clusters = NULL,
                              minsize = 8,
                              unoise_alpha = 2,
                              relabel = NULL,
                              relabel_sha1 = FALSE,
                              log_file = NULL,
                              threads = 1,
                              vsearch_options = NULL,
                              tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Only one table output allowed
  if (isTRUE(otutabout) && isTRUE(uc)) {
    stop("Only one of 'otutabout' and 'uc' can be TRUE. Please specify only one of these options.")
  }

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

  path_uc <- NULL
  if (!is.null(uc)) {
    path_uc <- if(is.character(uc)) uc else tempfile(pattern = "uc", tmpdir = tmpdir, fileext = ".uc")
    if (isTRUE(uc)) temp_files <- c(temp_files, path_uc)
  }

  # Build argument string for command line
  args <- c("--cluster_unoise", shQuote(fasta_file),
            "--threads", threads,
            "--minsize", minsize,
            "--unoise_alpha", unoise_alpha,
            "--centroids", centrfile,
            "--otutabout", otutabfile)

  if (!is.null(relabel))         args <- c(args, "--relabel", relabel)
  if (relabel_sha1)              args <- c(args, "--relabel_sha1", "")
  if (!is.null(vsearch_options)) args <- c(args, vsearch_options)
  if (!is.null(log_file))        args <- c(args, "--log", log_file)
  if (!is.null(clusters))        args <- c(args, "--clusters", clusters)
  if (!is.null(path_uc))         args <- c(args, "--uc", path_uc)

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Read results and make otu table
  if (file.size(centrfile) == 0) {
    centr.tbl <- tibble::tibble(Header = character(), Sequence = character())
  } else {
    centr.tbl <- microseq::readFasta(centrfile)
  }

  if(nrow(centr.tbl) > 0){
    centr.tbl <- centr.tbl |>
      dplyr::mutate(tag = stringr::word(Header, 1, sep = ";")) |>
      dplyr::distinct(tag, .keep_all = T) |>
      dplyr::select(tag, Sequence)

    otu.tbl <- suppressMessages(readr::read_tsv(otutabfile)) |>
      dplyr::rename(tag = `#OTU ID`)

    sizes <- otu.tbl |>
      dplyr::select(-tag) |>
      dplyr::select(tidyselect::where(is.numeric)) |>
      as.matrix() |>
      rowSums()

    otu.tbl <- otu.tbl |>
      dplyr::left_join(centr.tbl, by = "tag") |>
      dplyr::mutate(size = sizes) |>
      dplyr::arrange(dplyr::desc(size)) |>
      dplyr::mutate(Header = stringr::str_c(tag, ";size=", size)) |>
      dplyr::select(-tag, -size) |>
      dplyr::relocate(Header, Sequence)

    statistics <- calculate_unoise_statistics(otu.tbl,
                                              fasta_file,
                                              fasta_input_name)
    attr(otu.tbl, "statistics") <- statistics
  } else {
    warning("No clusters found, try to lower minsize")
    otu.tbl <- NULL
  }

  if (isTRUE(uc)) {
    res <- suppressMessages(readr::read_delim(path_uc, delim = "\t", col_names = FALSE))
    attr(res, "statistics") <- statistics
    return(res)
  }

  if (is.character(otutabout) || is.character(uc)) {
    if (is.character(otutabout) && !is.null(otu.tbl)) {
      readr::write_delim(otu.tbl, delim = "\t", file = otutabout)
    }
    return(invisible(NULL))
  }

  return(otu.tbl)
}

#' Calculate UNOISE statistics
#'
#' @description \code{calculate_unoise_statistics} calculates important
#' statistics after running \code{vs_cluster_unoise}, including the number of
#' clusters, sequences, and nucleotides, as well as the lengths and sizes of the
#' sequences and clusters.
#'
#' @param otu_tbl Output tibble from clustering. Contains the columns: Header,
#' Sequence, and one column for each cluster with read counts.
#' @param fasta_file File path to FASTA containing the input sequences to the
#' clustering.
#' @param fasta_input_name Name of the file/object with the input sequences
#' that was used in the clustering.
#'
#' @return A tibble with clustering statistics, including:
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
#' @return A tibble with clustering statistics.
#'
#' @noRd
#'
calculate_unoise_statistics <- function(otu_tbl,
                                        fasta_file,
                                        fasta_input_name) {

  # Process clustering output
  otu_tbl <- otu_tbl |>
    dplyr::mutate(cluster_size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
    dplyr::mutate(cluster_size = as.numeric(cluster_size)) |>
    dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

  # Make tibble from input sequences to the clustering
  input.df <- microseq::readFasta(fasta_file)

  # Calculate statistics
  num_nucleotides <- sum(nchar(input.df$Sequence))
  min_length_input_seq <- min(nchar(input.df$Sequence))
  max_length_input_seq <- max(nchar(input.df$Sequence))
  avg_length_input_seq <- mean(nchar(input.df$Sequence))
  num_clusters <- nrow(otu_tbl)
  min_size_cluster <- min(otu_tbl$cluster_size)
  max_size_cluster <- max(otu_tbl$cluster_size)
  avg_size_cluster <- round(mean(otu_tbl$cluster_size), 1)
  num_singletons <- sum(otu_tbl$cluster_size == 1)

  # Create table
  result_table <- tibble::tibble(
    num_nucleotides = num_nucleotides,
    min_length_input_seq = min_length_input_seq,
    max_length_input_seq = max_length_input_seq,
    avg_length_input_seq = avg_length_input_seq,
    num_clusters = num_clusters,
    min_size_cluster = min_size_cluster,
    max_size_cluster = max_size_cluster,
    avg_size_cluster = avg_size_cluster,
    num_singletons = num_singletons,
    input = fasta_input_name
  )

  return(result_table)
}

