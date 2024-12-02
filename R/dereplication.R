#' Dereplicate sequences
#'
#' @description Dereplication of sequences in FASTA file or object.
#'
#' @param fasta_input a FASTA file with reads or a FASTA object, see Detalis.
#' @param output name of the FASTA-file with the output or NULL, see Details.
#' @param fasta_width number of characters in the width of sequences in the output FASTA file. See Detalis.
#' @param minuniquesize minimum abundance value post-dereplication for sequence not to be discarded.
#' @param strand plus or both. When comparing sequences only check the plus strand or both strands.
#' @param sizein decides if abundance annotations present in sequence headers should be taken into account. True by default.
#' @param sizeout decides if abundance annotations should be added to FASTA headers.
#' @param relabel_sha1 relabel sequences using the SHA1 message digest algorithm.
#' @param relabel relabel sequences using the given prefix and a ticker to construct new headers
#'
#' @details Identical sequences in the input file/object are merged, using vsearch.
#' Identical sequences are defined as sequences with the same length and the same string of nucleotides.
#'
#' \code{fasta_input} can either be a FASTA file with reads or a FASTA object. The FASTA object needs to be a tibble
#' with columns \code{Header} and \code{Sequence}.
#'
#' If \code{output} is specified, the remaining sequences after quality filtering are output to this file in FASTA-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object, i.e. a tibble with
#' columns \code{Header} and \code{Sequence}.
#'
#' FASTA files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' @return a tibble with FASTA sequences with columns \code{Header} and \code{Sequence}.
#' @export
#'
vs_derep_fulllength <- function(fasta_input,
                                output = NULL,
                                minuniquesize = 1,
                                strand = "plus",
                                sizein = TRUE,
                                sizeout = TRUE,
                                relabel_sha1 = FALSE,
                                relabel = NULL,
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

  # Determine output file
  if (is.null(output)) {
    message("No filename for output file. No output file will be created.")
    outfile <- tempfile(pattern = "derep", fileext = ".fa")
  } else {
    # Validate output file extention
    validate_fasta_file(output)

    message("Writing filtered sequences to file: ", output)
    outfile <- output
  }

  # Build argument string for command line
  args <- c("--derep_fulllength", fasta_file,
            "--threads", 1,
            "--minuniquesize", minuniquesize,
            "--strand", strand,
            "--fasta_width", fasta_width,
            "--output", outfile)

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  if (relabel_sha1) {
    args <- c(args, "--relabel_sha1", "")
  }

  if (!is.null(relabel)) {
    args <- c(args, "--relabel", relabel)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Read output into FASTA object (tbl)
  derep_fasta <- microseq::readFasta(outfile)

  # Remove temp file for input if necessary
  if (!is.character(fasta_input)) {
    file.remove(fasta_file)
  }

  # Remove temp file for output if necessary
  if (is.null(output)) {
    file.remove(outfile)
  }

  return(list(derep_fasta = derep_fasta))
}
