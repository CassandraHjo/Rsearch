#' Trimming of reads
#'
#' @description Trimming sequences in FASTQ file or object
#'
#' @param fastq_input a FASTQ-file or a FASTQ object with forward reads (R1), see Detalis.
#' @param reverse a FASTQ-file or a FASTQ object with reverse reads (R2), see Detalis.
#' @param stripright_R1 number of bases stripped from the left end of the forward reads.
#' @param stripleft_R1 number of bases stripped from the right end of the forward reads.
#' @param stripright_R2 number of bases stripped from the left end of the reverse reads.
#' @param stripleft_R2 number of bases stripped from the right end of the reverse reads.
#' @param minlen minimum number of bases in input sequences.
#' @param threads number of computational threads to be used by vsearch.
#' @param fastqout_R1 name of the FASTQ-file with the output from trimming R1 reads or NULL, see Details.
#' @param fastqout_R2 name of the FASTQ-file with the output from trimming R2 reads or NULL, see Details.
#'
#' @details The reads in the input FASTQ-file (\code{fastq_input}) are trimmed based on the specified number of bases for each en of the read, using vsearch.
#'
#' \code{fastq_input} can either be a FASTQ file with reads or a FASTQ object. The FASTQ object needs to be a tibble
#' with columns \code{Header}, \code{Sequence} and \code{Quality} (like the one outputted from \code{vs_fastq_mergepairs()}).

#' \code{reverse} is an optional argument, for when you want to trim two files at a time (e.g. forward and reverse reads).
#' It can either be a FASTQ file with reads, a FASTQ object or \code{NULL}. The FASTQ object needs to be a tibble
#' with columns \code{Header}, \code{Sequence} and \code{Quality}.
#' If unspecified (\code{NULL}), only \code{fastq_input} is trimmed.
#'
#' If \code{fastqout_R1} is specified, the remaining sequences after trimming are output to this file in FASTQ-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTQ-object, i.e. a tibble with
#' columns \code{Header} and \code{Sequence}.
#'
#' #' If \code{fastqout_R2} is specified, the remaining sequences after trimming are output to this file in FASTQ-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTQ-object, i.e. a tibble with
#' columns \code{Header} and \code{Sequence}.
#'
#' @return A tibble, \code{trimmed_R1}, containing trimmed FASTQ sequences from the forward reads, with columns \code{Header}, \code{Sequence} and \code{Quality}.
#'
#' If \code{reverse} is specified, the resulting tibble (\code{trimmed_R2}) containing trimmed FASTQ sequences from the reverse reads, with columns \code{Header}, \code{Sequence} and \code{Quality} is an attribute to the primary table (\code{trimmed_R1}).
#' This table can be accessed by running \code{attributes(trimmed_R1)$trimmed_R2} or \code{attr(trimmed_R1, "trimmed_R2")}.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_trim <- function(fastq_input,
                          reverse = NULL,
                          stripright_R1 = 0,
                          stripleft_R1 = 0,
                          stripright_R2 = 0,
                          stripleft_R2 = 0,
                          minlen = 1,
                          threads = 1,
                          fastqout_R1 = NULL,
                          fastqout_R2 = NULL){


  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Check if FASTQ R1 input is file or tibble
  if (!is.character(fastq_input)){
    temp_file_R1 <- tempfile(pattern = "R1_reads", fileext = ".fq")
    temp_files <- c(temp_files, temp_file_R1)
    microseq::writeFastq(fastq_input, temp_file_R1)
    R1_file <- temp_file_R1
  } else {
    # Check is R1 input file exists at given path
    if (!file.exists(fastq_input)) stop("Cannot find input file: ", fastq_input)
    R1_file <- fastq_input
  }

  # Normalize file path
  R1_file <- normalizePath(R1_file)

  # Determine output file
  if (is.null(fastqout_R1)) {
    message("No filename for R1 output file. No output file will be created.")
    outfile <- tempfile(pattern = "trimmed_R1", fileext = ".fq")
    temp_files <- c(temp_files, outfile)
  } else {
    message("Writing trimmed R1 sequences to file: ", fastqout_R1)
    outfile <- fastqout_R1
  }

  # Build argument string for command line to trim R1 reads
  args <- c("--fastq_filter", R1_file,
            "--fastq_stripright", stripright_R1,
            "--fastq_stripleft", stripleft_R1,
            "--fastq_minlen", minlen,
            "--threads", threads,
            "--fastqout", outfile)

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Read output into fastq object (tbl)
  trimmed_R1 <- microseq::readFastq(outfile)

  # # Remove temp file if necessary
  # if (is.null(fastqout_R1)) {
  #   file.remove(outfile)
  # }

  # Check if FASTQ R2 input is given, and of which type (file or tibble)
  if (!is.null(reverse)){
    if (!is.character(reverse)){
      temp_file_R2 <- tempfile(pattern = "trimmed_R2", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_R2)
      microseq::writeFastq(reverse, temp_file_R2)
      R2_file <- temp_file_R2
    } else {
      R2_file <- reverse
    }
    # Normalize file path
    R2_file <- normalizePath(R2_file)

    # Determine output file
    if (is.null(fastqout_R2)) {
      message("No filename for R2 output file. No output file will be created.")
      outfile <- tempfile(pattern = "trimmed_R2", fileext = ".fq")
      temp_files <- c(temp_files, outfile)
    } else {
      message("Writing trimmed R2 sequences to file: ", fastqout_R2)
      outfile <- fastqout_R2
    }

    # Build argument string for command line to trim R2 reads
    args <- c("--fastq_filter", R2_file,
              "--fastq_stripright", stripright_R2,
              "--fastq_stripleft", stripleft_R2,
              "--fastq_minlen", minlen,
              "--threads", threads,
              "--fastqout", outfile)

    # Run vsearch
    vsearch_output <- system2(command = vsearch_executable,
                              args = args,
                              stdout = TRUE,
                              stderr = TRUE)

    # Read output into fastq object (tbl)
    trimmed_R2 <- microseq::readFastq(outfile)

    # Add R2 tibble as attribute to main tibble
    attr(trimmed_R1, "trimmed_R2") <- trimmed_R2

    # # Remove temp file if necessary
    # if (is.null(fastqout_R2)) {
    #   file.remove(outfile)
    # }
  }

  # Cleanup temporary files
  if (length(temp_files) > 0) {
    on.exit(
      file.remove(temp_files),
      add = TRUE)
  }

  return(trimmed_R1)
}

