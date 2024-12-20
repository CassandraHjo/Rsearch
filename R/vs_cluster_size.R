#' Cluster FASTA sequences
#'
#' @description Cluster FASTA sequences in the given file or object.
#'
#' @param fasta_input A FASTA file path or a FASTA object with reads to cluster. See Details.
#' @param centroids Name of the FASTA output file for the cluster centroid sequences. If \code{NULL} (default) no output will be written to file. See Details.
#' @param id The pairwise identity threshold for sequence to be added to cluster. Defaults to \code{0.97}. See Details.
#' @param strand \code{"plus"} (default) or \code{"both"}. When comparing sequences only check the \code{plus} strand or \code{both} strands.
#' @param sizein Decides if abundance annotations present in sequence headers should be taken into account. Defaults to \code{TRUE}.
#' @param sizeout Decides if abundance annotations should be added to FASTA headers. Defaults to \code{TRUE}.
#' @param relabel Relabel sequences using the given prefix and a ticker to construct new headers. Defaults to \code{"OTU"}.
#' @param threads The number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#' @param fasta_width The number of characters in the width of sequences in the output FASTA file. Defaults to \code{0}. See Details.
#'
#' @details Sequences in the input file are clustered, using \code{vsearch}´s \code{cluster_size}.
#' The function will automatically sort by decreasing sequence abundance beforehand.
#'
#' \code{fasta_input} can either be a FASTA file or object. FASTA objects are tibbles that contain the columns \code{Header} and \code{Sequence}.
#'
#' The centroids in \code{centroids} are the sequences that seeded the clusters (i.e. the first sequence of the cluster).
#' If \code{centroids} is specified, the remaining sequences after quality filtering are output to this file in FASTA format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object.
#'
#' \code{id} is a value between 0 and 1, and describes the the minimum pairwise identity with the centroid for sequence to be added to cluster.
#' The sequence is not added if pairwise identity is bellow \code{id}. The pairwise identity is defined as the number of (matching columns) / (alignment length - terminal gaps).
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' @importFrom magrittr %>%
#'
#' @return Tibble or \code{NULL}.
#'
#' If \code{centroids} is not specified, a FASTA object containing the centroid sequences is returned. If \code{centroids} is specified, results are written to file, and nothing is returned.
#'
#' When a FASTA object is returned, the statistics from the clustering, \code{statistics}, is an attribute, called \code{"statistics"} of the centroids tibble.
#' This tibble contains clustering statistics, including statistics about input sequences, number of clusters and their sizes.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fa")
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
#' @export
#'
vs_cluster_size <- function(fasta_input,
                            centroids = NULL,
                            id = 0.97,
                            strand = "plus",
                            sizein = TRUE,
                            sizeout = TRUE,
                            relabel = "OTU",
                            threads = 1,
                            fasta_width = 0){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0) {
      file.remove(temp_files)
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
  args <- c("--cluster_size", fasta_file,
            "--id", id,
            "--threads", 1,
            "--strand", strand,
            "--fasta_width", fasta_width,
            "--relabel", relabel,
            "--centroids", outfile)

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(centroids)) {

    # Read output into FASTA object (tbl)
    centroids_fasta <- microseq::readFasta(outfile) %>%
      dplyr::mutate(centroid_size = stringr::str_remove(Header, ".+;size=")) %>%
      dplyr::mutate(centroid_size = as.numeric(centroid_size)) %>%
      dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

    # Output statistics in table
    statistics <- parse_cluster_statistics(vsearch_output, fasta_input_name)

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

#' Parse clustering statistics from string to tibble
#'
#' @description This function transforms the output from \code{vsearch} when running \code{vs_cluster_size()} into a tibble.
#'
#' @param output A string of output from clustering sequences with \code{vsearch}.
#' @param input_file The name of the file/object with sequences used in the clustering
#'
#' @return A tibble with clustering metrics, including the number of nucleotides, sequences, clusters, and the lengths and sizes of the sequences and clusters.
#'
#' @noRd
#'
parse_cluster_statistics <- function(output, input_file) {

  # Extract values from output
  nucleotides <- as.numeric(stringr::str_extract(stringr::str_subset(output, " nt in"), "\\d+(?= nt )"))
  sequences <- as.numeric(stringr::str_extract(stringr::str_subset(output, "nt in"), "(?<=nt in )\\d+"))
  min_len_seq <- as.numeric(stringr::str_extract(stringr::str_subset(output, "nt in"), "(?<=min )\\d+"))
  max_len_seq <- as.numeric(stringr::str_extract(stringr::str_subset(output, "nt in"), "(?<=max )\\d+"))
  avg_len_seq <- as.numeric(stringr::str_extract(stringr::str_subset(output, "nt in"), "(?<=avg )\\d+"))

  num_clusters <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Clusters: "), "(?<=Clusters: )\\d+"))
  min_size_clusters <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Clusters: "), "(?<=min )\\d+"))
  max_size_clusters <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Clusters: "), "(?<=max )\\d+"))
  avg_size_clusters <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Clusters: "), "(?<=avg )\\d+\\.?\\d*"))

  num_singletons <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Singletons: "), "(?<=Singletons: )\\d+"))

  # Create table
  result_table <- data.frame(
    Tot_nucleotides = nucleotides,
    Num_sequences = sequences,
    Min_sequence_length = min_len_seq,
    Max_sequence_length = max_len_seq,
    Avg_sequence_length = avg_len_seq,
    Num_clusters = num_clusters,
    Min_size_clusters = min_size_clusters,
    Max_size_clusters = max_size_clusters,
    Avg_size_clusters = avg_size_clusters,
    Num_singletons = num_singletons,
    Input_file = basename(input_file)
  )

  return(result_table)
}

