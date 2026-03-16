#' Cluster FASTA sequences
#'
#' @description \code{vs_cluster_size} clusters FASTA sequences from a given
#' file or object using \code{VSEARCH}´s \code{cluster_size} method. The
#' function automatically sorts sequences by decreasing abundance before
#' clustering.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object containing
#' reads to cluster. See \emph{Details}.
#' @param centroids (Optional). A character string specifying the name of the
#' FASTA output file for the cluster centroid sequences. If \code{NULL}
#' (default), no output is written to a file and the centroid sequences are
#' returned as a FASTA object. See \emph{Details}.
#' @param otutabout (Optional). A character string specifying the name of the
#' output file in an OTU table format. If \code{NULL} (default), no output is
#' written to a file. If \code{TRUE}, the output is returned as a tibble. See
#' \emph{Details}.
#' @param clusters (Optional). Output each cluster to a separate FASTA file
#' using the prefix string provided and a ticker (0, 1, 2, etc.) to construct
#' the path and the filenames. Defaults to \code{NULL}, meaning no such files
#' are created.
#' @param uc (Optional). A character string specifying the name of the output
#' file in a uclust-like format. If \code{NULL} (default), no output is written
#' to a file. If \code{TRUE}, the output is returned as a tibble. See
#' \emph{Details}.
#' @param size_column (Optional). If \code{TRUE}, a column with the size of each
#' centroid is added to the centroid output tibble.
#' @param id (Optional). Pairwise identity threshold for sequence to be added to
#' a cluster. Defaults to \code{0.97}. See \emph{Details}.
#' @param strand (Optional). Specifies which strand to consider when comparing
#' sequences. Can be either \code{"plus"} (default) or \code{"both"}.
#' @param sizein (Optional). If \code{TRUE} (default), abundance annotations
#' present in sequence headers are taken into account.
#' @param sizeout (Optional). If \code{TRUE} (default), abundance annotations
#' are added to FASTA headers.
#' @param relabel (Optional). Relabel sequences using the given prefix and a
#' ticker to construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 (Optional). If \code{TRUE} (default), relabel sequences
#' using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Defaults to \code{0}, which eliminates wrapping.
#' @param sample (Optional). Add the given sample identifier string to sequence
#' headers. For instance, if the given string is "ABC", the text ";sample=ABC"
#' will be added to the header. This option is only applicable when the output
#' format is FASTA (\code{centroids}). If \code{NULL} (default), no identifier
#' is added.
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
#' Sequences are clustered based on the pairwise identity threshold specified by
#' \code{id}. Sequences are sorted by decreasing abundance before clustering.
#' The centroid of each cluster is the first sequence added to the cluster.
#'
#' \code{fasta_input} can either be a file path to a FASTA file or a FASTA
#' object. FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}, see \code{\link[microseq]{readFasta}}.
#'
#' If neither \code{centroids}, \code{otutabout}, nor \code{uc} is specified
#' (default), the function returns the centroid sequences as a FASTA object with
#' an additional column \code{otu_id}. This column contains the identifier
#' extracted from each sequence header.
#'
#' If \code{centroids} is specified, centroid sequences are written to the
#' specified file in FASTA format.
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
#' \code{id} is a value between 0 and 1 that defines the minimum pairwise
#' identity required for a sequence to be added to a cluster. A sequence is not
#' added to a cluster if its pairwise identity with the centroid is below the
#' \code{id} threshold.
#' Pairwise identity is calculated as the number of matching columns divided by
#' the alignment length minus terminal gaps.
#'
#' If \code{log_file} is \code{NULL} and \code{centroids} is specified,
#' clustering statistics from \code{VSEARCH} will not be captured.
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
#' If \code{otutabout} is \code{TRUE}, an OTU table is returned as a tibble.
#' If \code{otutabout} is a character string, the output is written to the file,
#' and no tibble is returned.
#'
#' If neither \code{centroids} nor \code{otutabout} is specified, a FASTA object
#' with the centroid sequences and additional column \code{otu_id} is returned.
#' The clustering statistics are included as an attribute named
#' \code{"statistics"}.
#'
#' The \code{"statistics"} attribute of the returned tibble (when
#' \code{centroids} is \code{NULL}) is a tibble with the following columns:
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
#' # Define arguments
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                                    "small.fasta")
#' centroids <- NULL
#'
#' # Cluster sequences and return a FASTA tibble
#' cluster_seqs <- vs_cluster_size(fasta_input = fasta_input,
#'                                 centroids = centroids)
#'
#' # Extract clustering statistics
#' statistics <- attr(cluster_seqs, "statistics")
#'
#' # Cluster sequences and write centroids to a file
#' vs_cluster_size(fasta_input = fasta_input,
#'                 centroids = "centroids_sequences.fa")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_cluster_size vs_cluster cluster_size cluster
#'
#' @export
#'
vs_cluster_size <- function(fasta_input,
                            centroids = NULL,
                            otutabout = NULL,
                            uc = NULL,
                            clusters = NULL,
                            size_column = FALSE,
                            id = 0.97,
                            strand = "plus",
                            sizein = TRUE,
                            sizeout = TRUE,
                            relabel = NULL,
                            relabel_sha1 = FALSE,
                            fasta_width = 0,
                            sample = NULL,
                            log_file = NULL,
                            threads = 1,
                            vsearch_options = NULL,
                            tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  # Validation: only one output table is allowed
  if (isTRUE(otutabout) && isTRUE(uc)) {
    stop("Only one table can be returned at a time ('otutabout' and 'uc' can not both be TRUE")
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

  # Determine output file based on user input

  # CENTROIDS: If centroids is specified, use that as output file. If centroids is NULL, use a temporary file.
  path_centroids <- NULL
  if (is.character(centroids) ) {
    path_centroids <- centroids
  } else if (is.null(centroids) && !isTRUE(otutabout) && !isTRUE(uc)) {
    path_centroids <- tempfile(pattern = "centroids",
                               tmpdir = tmpdir,
                               fileext = ".fa")
    temp_files <- c(temp_files, path_centroids)
  }

  # OTU table: If otutabout is specified, use that as output file. If otutabout is TRUE, use temp file.
  path_otutab <- NULL
  if (!is.null(otutabout)) {
    path_otutab <- ifelse(is.character(otutabout), otutabout, tempfile(pattern = "otutab",
                                                                       tmpdir = tmpdir,
                                                                       fileext = ".tsv"))
    if (isTRUE(otutabout)) temp_files <- c(temp_files, path_otutab)
  }

  # UC: If uc is specified, use that as output file. If uc is TRUE, use temp file.
  path_uc <- NULL
  if (!is.null(uc)) {
    path_uc <- ifelse(is.character(uc), uc, tempfile(pattern = "uc",
                                                     tmpdir = tmpdir,
                                                     fileext = ".txt"))
    if (isTRUE(uc)) temp_files <- c(temp_files, path_uc)
  }

  # Build argument string for command line
  args <- c("--cluster_size", shQuote(fasta_file),
            "--id", id,
            "--threads", threads,
            "--strand", strand,
            "--fasta_width", fasta_width)

  if (!is.null(path_centroids))  args <- c(args, "--centroids", path_centroids)
  if (!is.null(path_otutab))     args <- c(args, "--otutabout", path_otutab)
  if (!is.null(path_uc))         args <- c(args, "--uc", path_uc)

  if (!is.null(clusters))        args <- c(args, "--clusters", clusters)
  if (sizein)                    args <- c(args, "--sizein", "")
  if (sizeout)                   args <- c(args, "--sizeout", "")
  if (!is.null(relabel))         args <- c(args, "--relabel", relabel)
  if (relabel_sha1)              args <- c(args, "--relabel_sha1", "")
  if (!is.null(vsearch_options)) args <- c(args, vsearch_options)
  if (!is.null(sample))          args <- c(args, "--sample", sample)
  if (!is.null(log_file))        args <- c(args, "--log", log_file)

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Determine return output
  if (isTRUE(otutabout)) {
    df <- suppressMessages(readr::read_delim(path_otutab, delim = "\t"))
    colnames(df)[1] <- "otu_id"
    return(df)
  }

  if (isTRUE(uc)) {
    df <- suppressMessages(readr::read_delim(path_uc, delim = "\t", col_names = FALSE))
    return(df)
  }

  if (!is.null(centroids) || is.character(otutabout) || is.character(uc)) {
    return(invisible(NULL))
  }

  if (file.exists(path_centroids) && file.info(path_centroids)$size > 0) {
    centroids_fasta <- microseq::readFasta(path_centroids) |>
      dplyr::mutate(otu_id = stringr::str_extract(Header, "^[^;]+"))
  } else {
    centroids_fasta <- tibble::tibble(Header = character(), Sequence = character())
    warning("No centroid sequences were returned by VSEARCH. Check input quality or parameters.")
  }

  # Håndtering av size_column
  if (size_column && nrow(centroids_fasta) > 0) {
    centroids_fasta <- centroids_fasta |>
      dplyr::mutate(centroid_size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
      dplyr::mutate(centroid_size = as.numeric(centroid_size))
  }

  # Beregning av statistikk (med 0-tabell fallback)
  if (nrow(centroids_fasta) > 0) {
    statistics <- calculate_cluster_statistics(centroids_fasta,
                                               fasta_file,
                                               fasta_input_name)
  } else {
    statistics <- tibble::tibble(num_nucleotides = 0,
                                 min_length_input_seq = 0,
                                 max_length_input_seq = 0,
                                 avg_length_input_seq = 0,
                                 num_clusters = 0,
                                 min_size_cluster = 0,
                                 max_size_cluster = 0,
                                 avg_size_cluster = 0,
                                 num_singletons = 0,
                                 input = fasta_input_name)
  }

  attr(centroids_fasta, "statistics") <- statistics
  return(centroids_fasta)
}

#' Calculate clustering statistics
#'
#' @description \code{calculate_cluster_statistics} calculates important
#' clustering statistics after running \code{vs_cluster_size}, including the
#' number of clusters, sequences, and nucleotides, as well as the lengths and
#' sizes of the sequences and clusters.
#'
#' @param centroids_fasta Output tibble from clustering with centroids.
#' Contains the columns: Header, Sequence, and centroid_size (if
#' \code{size_column} is specified in \code{vs_cluster_size}).
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
calculate_cluster_statistics <- function(centroids_fasta,
                                         fasta_file,
                                         fasta_input_name) {

  # Process clustering output
  if (!"centroid_size" %in% colnames(centroids_fasta)) {
    centroids_fasta <- centroids_fasta |>
      dplyr::mutate(centroid_size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
      dplyr::mutate(centroid_size = as.numeric(centroid_size)) |>
      dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))
  }
  # Make tibble from input sequences to the clustering
  input.df <- microseq::readFasta(fasta_file)

  # Calculate statistics
  num_nucleotides <- sum(nchar(input.df$Sequence))
  min_length_input_seq <- min(nchar(input.df$Sequence))
  max_length_input_seq <- max(nchar(input.df$Sequence))
  avg_length_input_seq <- mean(nchar(input.df$Sequence))
  num_clusters <- nrow(centroids_fasta)
  min_size_cluster <- min(centroids_fasta$centroid_size)
  max_size_cluster <- max(centroids_fasta$centroid_size)
  avg_size_cluster <- round(mean(centroids_fasta$centroid_size), 1)
  num_singletons <- sum(centroids_fasta$centroid_size == 1)

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

