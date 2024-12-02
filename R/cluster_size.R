#' Clusterize FASTA sequences
#'
#' @description Clustering the FASTA sequences in the given file or object.
#'
#' @param fasta_input a file with FASTA sequences or a FASTA object, see Details.
#' @param centroids file name for the FASTA file for output cluster centroid sequences. See Details.
#' @param id pairwise identity threshold for sequence to be added to cluster. See Details.
#' @param strand plus or both. When comparing sequences only check the plus strand or both strands.
#' @param sizein decides if abundance annotations present in sequence headers should be taken into account. True by default.
#' @param sizeout decides if abundance annotations should be added to FASTA headers.
#' @param relabel relabel sequences using the given prefix and a ticker to construct new headers.
#' @param threads number of computational threads to be used by vsearch.
#' @param fasta_width number of characters in the width of sequences in the output FASTA file. See Details.
#'
#' @details Sequences in the input file are clustered, using vsearch´s cluster_size.
#' The function will automatically sort by decreasing sequence abundance beforehand.
#'
#' \code{fasta_input} can either be a file with FASTA sequences or a FASTA object. The FASTA object needs to be a tibble
#' with columns \code{Header} and \code{Sequence}.
#'
#' The centroids in \code{centroids} are the sequences that seeded the clusters (i.e. the first sequence of the cluster).
#' If \code{centroids} is specified, the remaining sequences after quality filtering are output to this file in FASTA-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object, i.e. a tibble with
#' columns \code{Header}, \code{centroid_size} and \code{Sequence}.
#'
#' \code{id} is a value between 0 and 1, and describes the the minimum pairwise identity with the centroid for sequence to be added to cluster.
#' The sequence is not added if pairwise identity is bellow \code{id}. The pairwise identity is defined as the number of (matching columns) / (alignment length - terminal gaps).
#'
#' FASTA files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' @importFrom magrittr %>%
#'
#' @return SOMETHING....
#' @references https://github.com/torognes/vsearch
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

  # Check if FASTA input is file or tibble
  if (!is.character(fasta_input)){
    temp_file <- tempfile(pattern = "input", fileext = ".fa")
    microseq::writeFasta(fasta_input, temp_file)
    fasta_file <- temp_file
  } else {
    fasta_file <- fasta_input
  }

  # Check is input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)

  # Determine centroids file
  if (is.null(centroids)) {
    message("No filename for centroids file. No centroids file will be created.")
    outfile <- tempfile(pattern = "centroids", fileext = ".fa")
  } else {
    # Validate centroids file extention
    validate_fasta_file(centroids)

    message("Writing filtered sequences to file: ", centroids)
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

  # Read output into FASTA object (tbl)
  centroids_fasta <- microseq::readFasta(outfile) %>%
    dplyr::mutate(centroid_size = stringr::str_remove(Header, ".+;size=")) %>%
    dplyr::mutate(centroid_size = as.numeric(centroid_size))

  # LEGGE INN STATISTIKK HER!

  # Remove temp file for input if necessary
  if (!is.character(fasta_input)) {
    file.remove(fasta_file)
  }

  # Remove temp file for output if necessary
  if (is.null(centroids)) {
    file.remove(outfile)
  }

  return(list(centroids_fasta = centroids_fasta))
}
