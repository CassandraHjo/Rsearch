#' Cluster FASTA sequences
#'
#' @description Clusters FASTA sequences in the given file or object.
#'
#' @param fasta_input A FASTA file path or a FASTA object with reads to cluster.
#' See Details.
#' @param centroids Name of the FASTA output file for the cluster centroid
#' sequences. If \code{NULL} (default) no output will be written to file.
#' See Details.
#' @param id The pairwise identity threshold for sequence to be added to
#' cluster. Defaults to \code{0.97}. See Details.
#' @param strand \code{"plus"} (default) or \code{"both"}. When comparing
#' sequences only check the \code{plus} strand or \code{both} strands.
#' @param sizein Decides if abundance annotations present in sequence headers
#' should be taken into account. Defaults to \code{TRUE}.
#' @param sizeout Decides if abundance annotations should be added to
#' FASTA headers. Defaults to \code{TRUE}.
#' @param relabel Relabel sequences using the given prefix and a ticker to
#' construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 Relabel sequences using the SHA1 message digest
#' algorithm. Defaults to \code{FALSE}.
#' @param threads The number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param fasta_width The number of characters in the width of sequences in the
#' output FASTA file. Defaults to \code{0}. See Details.
#' @param log_file Name of the log file to capture messages from \code{VSEARCH}.
#' If \code{NULL}, no log file is created. Defaults to \code{NULL}.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See Details.
#'
#' @details FASTA sequences in the input file are clustered,
#' using \code{VSEARCH}´s \code{cluster_size}.The function will automatically
#' sort by decreasing sequence abundance beforehand.
#'
#' \code{fasta_input} can either be a FASTA file or object. FASTA objects are
#' tibbles that contain the columns \code{Header} and \code{Sequence}.
#'
#' The centroids in \code{centroids} are the sequences that seeded the clusters
#' (i.e. the first sequence of the cluster). If \code{centroids} is specified,
#' the centroid sequences are output to this file in FASTA format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object.
#'
#' \code{id} is a value between 0 and 1, and describes the the minimum pairwise
#' identity with the centroid for sequence to be added to cluster.
#' The sequence is not added if pairwise identity is bellow \code{id}.
#' The pairwise identity is defined as the number of
#' (matching columns) / (alignment length - terminal gaps).
#'
#' FASTA files produced by \code{VSEARCH} are wrapped
#' (sequences are written on lines of integer nucleotides).\code{fasta_width} is
#' by default set to zero to eliminate the wrapping.
#'
#' \code{vsearch_options} can be used to pass additional arguments to
#' \code{VSEARCH}, that are not implemented in \code{Rsearch}. See the
#' \code{VSEARCH} manual for additional arguments, and how to use them.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{centroids} is unspecified, a FASTA object containing the centroid
#' sequences is returned. If \code{centroids} is specified, results are written
#' to file, and nothing is returned.
#'
#' When a FASTA object is returned, the statistics from the clustering,
#' \code{statistics}, is an attribute of the centroids tibble.
#' The statistics tibble has the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: The total number of nucleotides used as input
#'   for clustering.
#'   \item \code{min_length_input_seq}: The length of the shortest sequence used
#'   as input for clustering.
#'   \item \code{max_length_input_seq}: The length of the longest sequence used
#'   as input for clustering.
#'   \item \code{avg_length_input_seq}: The average length of the sequences used
#'   as input for clustering.
#'   \item \code{num_clusters}: The number of clusters generated.
#'   \item \code{min_size_cluster}: The size of the smallest cluster.
#'   \item \code{max_size_cluster}: The size of the largest cluster.
#'   \item \code{avg_size_cluster}: The average size of the clusters.
#'   \item \code{num_singletons}: The number of singletons after clustering.
#'   \item \code{input}: The name of the input file/object for the clustering.
#' }
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                                    "R1_sample1_small.fa")
#' centroids <- NULL
#'
#' # Cluster sequences, and return fasta tibble
#' cluster_seqs <- vs_cluster_size(fasta_input = fasta_input,
#'                                 centroids = centroids)
#'
#' # Extract clustering statistics
#' statistics <- attr(cluster_seqs, "statistics")
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
                            id = 0.97,
                            strand = "plus",
                            sizein = TRUE,
                            sizeout = TRUE,
                            relabel = NULL,
                            relabel_sha1 = FALSE,
                            threads = 1,
                            fasta_width = 0,
                            log_file = NULL,
                            vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
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
    temp_file <- tempfile(pattern = "input", fileext = ".fa")
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

  # Check is input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)

  # Determine centroids file
  if (is.null(centroids)) {
    outfile <- tempfile(pattern = "centroids", fileext = ".fa")
    temp_files <- c(temp_files, outfile)
  } else {
    outfile <- centroids
  }

  # Build argument string for command line
  args <- c("--cluster_size", shQuote(fasta_file),
            "--id", id,
            "--threads", 1,
            "--strand", strand,
            "--fasta_width", fasta_width,
            "--centroids", outfile)

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  # Add relabeling arguments if specified
  if (!is.null(relabel)){
    args <- c(args, "--relabel", relabel)
  }

  if (relabel_sha1){
    args <- c(args, "--relabel_sha1", "")
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(centroids)) {

    # Read output into FASTA object (tbl)
    centroids_fasta <- microseq::readFasta(outfile) |>
      dplyr::mutate(centroid_size = stringr::str_remove(Header, ".+;size=")) |>
      dplyr::mutate(centroid_size = as.numeric(centroid_size)) |>
      dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

    # Create statistics tibble
    statistics <- calculate_cluster_statistics(centroids_fasta,
                                               fasta_file,
                                               fasta_input_name)

    # Add additional tables as attributes to the primary table
    attr(centroids_fasta, "statistics") <- statistics
  }

  # Return results
  if (is.null(centroids)) { # Return tibble
    return(centroids_fasta)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}

#' Calculate clustering statistics
#'
#' @description Calculates important clustering statistics after running
#' \code{vs_cluster_size()}, including the number of clusters, sequences, and
#' nucleotides, and the lengths and sizes of the sequences and clusters.
#'
#' @param centroids_fasta The output tibble from clustering with the centroids.
#' Contains the columns: Header, Sequence, and centroid_size.
#' @param fasta_file The FASTA file containing the input sequences to the
#' clustering.
#' @param fasta_input_name The name of the file/object with the input sequences
#' that was used in the clustering.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: The total number of nucleotides used as input
#'   for clustering.
#'   \item \code{min_length_input_seq}: The length of the shortest sequence used
#'   as input for clustering.
#'   \item \code{max_length_input_seq}: The length of the longest sequence used
#'   as input for clustering.
#'   \item \code{avg_length_input_seq}: The average length of the sequences used
#'   as input for clustering.
#'   \item \code{num_clusters}: The number of clusters generated.
#'   \item \code{min_size_cluster}: The size of the smallest cluster.
#'   \item \code{max_size_cluster}: The size of the largest cluster.
#'   \item \code{avg_size_cluster}: The average size of the clusters.
#'   \item \code{num_singletons}: The number of singletons after clustering.
#'   \item \code{input}: The name of the input file/object for the clustering.
#' }
#'
#' @return A tibble with clustering statistics.
#'
#' @noRd
#'
calculate_cluster_statistics <- function(centroids_fasta,
                                         fasta_file,
                                         fasta_input_name) {

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
  result_table <- data.frame(
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

