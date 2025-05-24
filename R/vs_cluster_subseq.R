#' Cluster FASTA sequences
#'
#' @description \code{vs_cluster_subseq} clusters FASTA sequences from a given
#' file or object using \code{VSEARCH}´s \code{cluster_fast} method and 100%
#' identity. The function automatically sorts sequences by decreasing length
#' before clustering.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object containing
#' reads to cluster. See \emph{Details}.
#' @param centroids (Optional). A character string specifying the name of the
#' FASTA output file for the cluster centroid sequences. If \code{NULL}
#' (default), no output is written to a file and the centroid sequences are
#' returned as a FASTA object. See \emph{Details}.
#' @param strand (Optional). Specifies which strand to consider when comparing
#' sequences. Can be either \code{"plus"} (default) or \code{"both"}.
#' @param sizein (Optional). If \code{TRUE} (default), abundance annotations
#' present in sequence headers are taken into account.
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Defaults to \code{0}, which eliminates wrapping.
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
#' After merging/dereplication some sequences may be sub-sequences of longer
#' sequences. This function will cluster such sequences at 100% identity
#' (terminal gaps ignored), and keep the longest in each cluster as the
#' centroid.
#'
#' \code{fasta_input} can either be a file path to a FASTA file or a FASTA
#' object. FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}, see \code{\link{readFasta}}.
#'
#' If \code{sizein = TRUE} (default) the FASTA headers must contain text
#' matching the regular expression \code{"size=[0-9]+"} indicating the copy
#' number (=size) of each input sequence. This is then summed for each cluster
#' and added to the output. This text is typically added by de-replication, see
#' \code{\link{vs_fastx_uniques}}.
#'
#' The number of distinct sequences in each cluster is output as \code{members}.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{centroids} is specified the centroid sequences are written to the
#' specified file, and no tibble is returned.
#'
#' If \code{centroids} is not specified, a FASTA object
#' is returned. This is a \code{tibble} with columns \code{Header} and
#' \code{Sequence}, and also the additional column(s) \code{members} and, if
#' \code{sizein = TRUE}, \code{size}.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                                    "small.fasta")
#'
#' # De-replicating
#' derep.tbl <- vs_fastx_uniques(fasta_input, output_format = "fasta")
#'
#' # Clustering subsequences
#' cluster.tbl <- vs_cluster_subseq(fasta_input = derep.tbl)
#'
#' # Cluster sequences and write centroids to a file
#' vs_cluster_subseq(fasta_input = derep.tbl,
#'                   centroids = "distinct.fa")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_cluster_length
#'
#' @export
#'
vs_cluster_subseq <- function(fasta_input,
                              centroids = NULL,
                              strand = "plus",
                              sizein = TRUE,
                              fasta_width = 0,
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

  ##########################################
  # The arguments needed to build the args

  # Check if fasta_input is a file or a tibble
  if (!is.character(fasta_input)){
    temp_file <- tempfile(pattern = "input",
                          tmpdir = tmpdir,
                          fileext = ".fa")
    temp_files <- c(temp_files, temp_file)
    microseq::writeFasta(fasta_input, temp_file)
    fasta_file <- temp_file
  } else {
    fasta_file <- fasta_input
  }

  # Check if input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  ################################
  # Building the command line

  # The temporary UC-file
  uc_file <- tempfile(pattern = "uc",
                      tmpdir = tmpdir,
                      fileext = ".txt")
  temp_files <- c(temp_files, uc_file)

  # Build argument string for command line
  args <- c("--cluster_fast", shQuote(fasta_file),
            "--id", 1.0,
            "--threads", threads,
            "--strand", strand,
            "--uc", uc_file)

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Read uc_file and wrangle
  uc.tbl <- suppressMessages(readr::read_tsv(uc_file,
                                             col_names = c("type", "cluster",
                                                           "length", "identity",
                                                           "strand", "star1",
                                                           "star2", "star3",
                                                           "member", "centroid")
  )) |>
    dplyr::filter(type != "C") |>
    dplyr::select(centroid, member) |>
    dplyr::mutate(centroid = ifelse(centroid == "*", member, centroid)) |>
    dplyr::left_join(microseq::readFasta(fasta_file),
                     by = c("centroid" = "Header"))
  if (sizein) {
    out.tbl <- uc.tbl |>
      dplyr::mutate(size = stringr::str_extract(member, "size=[0-9]+")) |>
      dplyr::mutate(size = as.numeric(stringr::str_remove(size, "size="))) |>
      dplyr::group_by(centroid) |>
      dplyr::summarise(Sequence = Sequence[1],
                       members = dplyr::n(),
                       size = sum(size)) |>
      dplyr::ungroup() |>
      dplyr::mutate(centroid = stringr::str_replace(centroid,
                                                    "size=[0-9]+",
                                                    paste0("size=", size))) |>
      dplyr::mutate(centroid = paste0(centroid, ";members=", members))
  } else {
    out.tbl <- uc.tbl |>
      dplyr::group_by(centroid) |>
      dplyr::summarise(Sequence = Sequence[1], members = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(centroid = paste0(centroid, ";members=", members))
  }
  out.tbl <- out.tbl |>
    dplyr::rename(Header = centroid) |>
    dplyr::relocate(Header) |>
    dplyr::relocate(Sequence, .after = dplyr::last_col())

  # Determine return output
  if (!is.null(centroids)) {
    microseq::writeFasta(out.tbl, out.file = centroids, width = fasta_width)
    return(invisible(NULL)) # No return if centroids is specified
  } else {
    return(out.tbl)
  }
}
