#' Dereplicate sequences
#'
#' @description Dereplication of sequences in FASTA/FASTQ file or object.
#'
#' @param fastx_input A FASTQ/FASTA file path or object. See Details.
#' @param fastx_output Name of the output file for dereplicated reads from \code{fastx_input}. File can be in either FASTA or FASTQ format, depending on \code{output_format}. If \code{NULL} (default) no sequences will be written to file. See Details.
#' @param input_format Format of input file or object \code{fastx_input}: \code{"fasta"} or \code{"fastq"} (default).
#' @param output_format The desired output format for file/tibble: \code{"fasta"} or \code{"fastq"} (default).
#' @param minuniquesize The minimum abundance value post-dereplication for a sequence not to be discarded. Defaults to \code{1}.
#' @param strand \code{"plus"} (default) or \code{"both"}. When comparing sequences only check the plus strand or both strands.
#' @param sizein Decides if abundance annotations present in sequence headers in \code{fastx_input} should be taken into account. Defaults to \code{TRUE}.
#' @param sizeout Decides if abundance annotations should be added to headers in output table or output file (\code{fastx_output}). Defaults to \code{TRUE}.
#' @param relabel_sha1 Relabel sequences using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param relabel Relabel sequences using the given prefix and a ticker to construct new headers. Defaults to \code{NULL}.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param fastq_qout_max If \code{TRUE}, the quality score will be the highest (best) quality score observed in each position. Defaults to \code{FALSE}.
#'
#' @details The reads in the input file/object (\code{fastx_input}) are dereplicated by merging identical sequences, using \code{vsearch}.
#' Identical sequences are defined as sequences with the same length and the same string of nucleotides (case insensitive, T and U are considered the same).
#'
#' \code{fastx_input} can either be a FASTA/FASTQ file or object with reads. FASTA objects are tibbles that contain the columns \code{Header} and \code{Sequence}. FASTQ objects are tibbles that contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' If \code{fastx_output} is specified, the dereplicated sequences are output to this file in format given by \code{output_format}.
#' If unspecified (\code{NULL}) the result is returned as a FASTA/FASTQ object, depending on \code{output_format}.
#'
#' FASTA files produced by\code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' By default, the quality scores in FASTQ output files will correspond to the average error probability of the nucleotides in the each position.
#' If \code{fastq_qout_max = TRUE}, the quality score will be the highest (best) quality score observed in each position.
#'
#' @return Tibble or \code{NULL}.
#'
#' If \code{fastx_output} is not specified, a tibble containing the dereplicated reads is returned. If \code{fastx_output} is specified nothing is returned.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' fastx_output <- NULL
#' input_format <- "fastq"
#' output_format <- "fastq"
#'
#' # Dereplicate sequences, with tibble as output
#' derep_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
#'                              fastx_output = fastx_output,
#'                              input_format = input_format,
#'                              output_format = output_format)
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_fastx_uniques vs_fastq_uniques vs_fasta_uniques vs_fastx_dereplication
#'
#' @export
#'
vs_fastx_uniques <- function(fastx_input,
                             fastx_output = NULL,
                             input_format = "fastq",
                             output_format = "fastq",
                             minuniquesize = 1,
                             strand = "plus",
                             sizein = TRUE,
                             sizeout = TRUE,
                             relabel_sha1 = FALSE,
                             relabel = NULL,
                             fasta_width = 0,
                             fastq_qout_max = FALSE){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate input_format
  if (!input_format %in% c("fasta", "fastq")) {
    stop("Invalid input_format. Choose from fasta or fastq.")
  }

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
  }

  if (input_format == "fasta" && output_format == "fastq") {
    stop("Invalid output_format when input_format is 'fasta'")
  }

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  # Create empty vector for collecting temporary files
  temp_files <- character()

  # Set up cleanup of temporary files
  on.exit({
    existing_files <- temp_files[file.exists(temp_files)]
    if (length(existing_files) > 0) {
      file.remove(existing_files)
    }
  }, add = TRUE)

  # Handle input_format = "fasta"
  if (input_format == "fasta") {
    # Handle input: file or tibble
    if (!is.character(fastx_input)){
      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }
      temp_file <- tempfile(pattern = "input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file)
      microseq::writeFasta(fastx_input, temp_file)
      input_file <- temp_file
    } else {
      input_file <- fastx_input
    }
  }

  # Handle input_format = "fastq"
  if (input_format == "fastq") {
    # Handle input: file or tibble
    if (!is.character(fastx_input)){
      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }
      temp_file <- tempfile(pattern = "input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file)
      microseq::writeFastq(fastx_input, temp_file)
      input_file <- temp_file
    } else {
      input_file <- fastx_input
    }
  }


    # Handle output_format = "fasta"
  if (output_format == "fasta"){
    if (is.null(fastx_output)) {
      output_file <- tempfile(pattern = "derep", fileext = ".fa")
      temp_files <- c(temp_files, output_file)
    } else {
      output_file <- fastx_output
    }
  }

    # Handle output_format = "fastq"
    if (output_format == "fastq"){
    if (is.null(fastx_output)) {
      output_file <- tempfile(pattern = "derep", fileext = ".fq")
      temp_files <- c(temp_files, output_file)
    } else {
      output_file <- fastx_output
    }
  }

  # Check is input file exists at given path
  if (!file.exists(input_file)) stop("Cannot find input file: ", input_file)

  # Normalize file paths
  input_file <- normalizePath(input_file) |>
    shQuote()

  # Build argument string for command line
  args <- c("--fastx_uniques", input_file,
            "--threads", 1,
            "--minuniquesize", minuniquesize,
            "--strand", strand)

  if (output_format == "fasta") {
    args <- c(args,
              "--fasta_width", fasta_width,
              "--fastaout", output_file)
  }

  if (output_format == "fastq") {
    args <- c(args, "--fastqout", output_file)
  }

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

  if (fastq_qout_max) {
    args <- c(args, "--fastq_qout_max", "")
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(fastx_output)) {
    if (output_format == "fasta"){
      # Read output into FASTA object
      derep_tbl <- microseq::readFasta(output_file)
    }
    if (output_format == "fastq"){
      # Read output into FASTQ object
      derep_tbl <- microseq::readFastq(output_file)
    }
  }

  # Return results
  if (is.null(fastx_output)) { # Return tibble
    return(derep_tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
