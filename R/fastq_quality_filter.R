#' Quality filtering
#'
#' @description Quality filtering of sequences in FASTQ file or object
#'
#' @param fastq_input a FASTQ file with reads or a FASTQ object, see Detalis.
#' @param fastaout name of the FASTA-file with the output or NULL, see Details.
#' @param fastq_maxee_rate threshold for average expected error. Value ranging form 0.0 to 1.0. See Details.
#' @param fasta_width number of characters in the width of sequences in the output FASTA file. See Detalis.
#' @param threads number of computational threads to be used by vsearch.
#' @param log_file name of the log file with messages from running vsearch or NULL, see Details.
#'
#' @details The reads in the input FASTQ-file (\code{fastq_file}) are filtered based on ..., using vsearch.
#'
#' \code{fastq_input} can either be a FASTQ file with reads or a FASTQ object. The FASTQ object needs to be a tibble
#' with columns \code{Header}, \code{Sequence} and \code{Quality} (like the one outputed from \code{vs_fastq_mergepairs()}).
#'
#' If \code{fastaout} is specified, the remaining sequences after quality filtering are output to this file in FASTA-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object, i.e. a tibble with
#' columns \code{Header} and \code{Sequence}.
#'
#' Sequences with an average expected error greater than the specified \code{fastq_maxee_rate} are discarded.
#' For a given sequence, the average expected error is the sum of error probabilities for all the positions in the sequence,
#' divided by the length of the sequence.
#'
#' FASTA files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' #' If \code{log_file} is specified, the messages are output to this file.
#' If unspecified (\code{NULL}) no log file is written.
#'
#' @return a tibble with FASTA sequences with columns \code{Header} and \code{Sequence}.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_filter <- function(fastq_input,
                            fastaout = NULL,
                            fastq_maxee_rate,
                            fasta_width = 0,
                            threads = 1,
                            log_file = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Check if FASTQ input is file or tibble
  if (!is.character(fastq_input)){
    temp_file <- tempfile(pattern = "merged", fileext = ".fq")
    microseq::writeFastq(fastq_input, temp_file)
    fastq_file <- temp_file
  } else {
    fastq_file <- fastq_input
  }

  # Check is input file exists at given path
  if (!file.exists(fastq_file)) stop("Cannot find input file: ", fastq_file)

  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)

  # Determine output file
  if (is.null(fastaout)) {
    message("No filename for output file. No output file will be created.")
    outfile <- tempfile(pattern = "filtered", fileext = ".fa")
  } else {
    # Validate output file extention
    validate_fasta_file(fastaout)

    message("Writing filtered sequences to file: ", fastaout)
    outfile <- fastaout
  }

  # Build argument string for command line
  args <- c("--fastq_filter", fastq_file,
            "--threads", threads,
            "--fastq_maxee_rate", fastq_maxee_rate,
            "--fasta_width", fasta_width,
            "--fastaout", outfile)

  # Add log file if specified by user
  if (!is.null(log_file)) {
    validate_log_file(log_file)
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Read output into FASTA object (tbl)
  filt_fasta <- microseq::readFasta(outfile)

  # Read output from log file if log file was specified
  if (!is.null(log_file)) {
    vsearch_output <- readLines(log_file)
  }

  # Output statistics in table
  statistics <- parse_filter_statistics(vsearch_output, fastq_file)

  # Remove temp file for input if necessary
  if (!is.character(fastq_input)) {
    file.remove(fastq_file)
  }

  # Remove temp file for output if necessary
  if (is.null(fastaout)) {
    file.remove(outfile)
  }

  return(list(statistics = statistics, filt_fasta = filt_fasta))
}

#' Validate outfile extention
#'
#' @param fastaout
#'
#' @noRd
validate_fasta_file <- function(fastaout) {
  allowed_extensions <- c("\\.fa$", "\\.fasta$")
  if (!any(stringr::str_ends(fastaout, allowed_extensions))) {
    stop("The output file needs one of the following extentions: .fa, .fasta.")
  }
}

#' Parse statistics from filtering sequences to tibble
#'
#' @param output string of output from running vs_fastq_filter
#'
#' @return table with filtering metrics
#' @noRd
parse_filter_statistics <- function(output, fastq_file) {

  # Find line with statistics
  stats_line <- stringr::str_subset(output, "sequences kept")

  # Extract number of kept sequences
  kept <- as.numeric(stringr::str_extract(stats_line, "^\\d+"))

  # Extract number of truncated sequences
  truncated <- as.numeric(stringr::str_extract(stats_line, "(?<=of which )\\d+(?= truncated)"))

  # Extract number of discarded sequences
  discarded <- as.numeric(stringr::str_extract(stats_line, "(?<=, )\\d+(?= sequences discarded)"))

  # Create table
  result_table <- data.frame(
    Kept_Sequences = kept,
    Truncated_Sequences = truncated,
    Discarded_Sequences = discarded,
    fastq_file = basename(fastq_file)
  )

  return(result_table)
}
