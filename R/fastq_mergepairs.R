# TODO: Skrive dokumentasjon
# TODO: Skrive tester

#' Merge read-pairs
#'
#' @description Merging read-pairs with overlapping regions.
#'
#' @param fastq_input A FASTQ file path containing (forward) reads or a FASTQ object (tibble), see Details.
#' @param reverse A FASTQ file path containing reverse reads or a FASTQ object (tibble), see Details.
#' @param fastqout Name of the FASTQ output file for merged sequences. If \code{NULL} no FASTQ output file will be written to file. See Details.
#' @param log_file Name of the log file to capture messages from vsearch. If \code{NULL}, no log file is created.
#' @param threads Number of computational threads to be used by vsearch.
#'
#' @details The read-pairs in the input FASTQ-files (\code{fastq_input} and \code{reverse}) are merged if they have sufficient overlap, using vsearch.
#'
#' \code{fastq_input} and \code{reverse} can either be FASTQ files or FASTQ objects. If provided as tibbles, they must contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' If \code{fastqout} is specified, the merged reads are output to this file in FASTQ-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTQ-object, i.e. a tibble with columns \code{Header}, \code{Sequence} and \code{Quality}.
#'
#' If \code{log_file} is specified, the messages are output to this file. If unspecified (\code{NULL}) no log file is written.
#'
#' @return A tibble, \code{merged_fastq}, containing the merged FASTQ sequences, with columns \code{Header}, \code{Sequence} and \code{Quality}.
#'
#' The statistics from the merging, \code{statistics}, is an attribute of \code{merged_fastq}. This tibble contains merging statistics, including number of pairs, number of merged pairs, and length metrics.
#' The statistics can be accessed by running \code{attributes(merged_fastq)$statistics} or \code{attr(merged_fastq, "statistics")}.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_mergepairs <- function(fastq_input,
                                reverse,
                                fastqout = NULL,
                                log_file = NULL,
                                threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Handle input: file or tibble
  if (!is.character(fastq_input)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_fastq_file <- tempfile(pattern = "fastq_input_temp_", fileext = ".fq")
    temp_files <- c(temp_files, temp_fastq_file)
    microseq::writeFastq(fastq_input, temp_fastq_file)
    fastq_file <- temp_fastq_file

    # Capture original name for statistics table later
    fastq_input_name <- as.character(substitute(fastq_input))
  } else {
    fastq_file <- fastq_input

    # Capture original name for statistics table later
    fastq_input_name <- basename(fastq_input)
  }

  # Handle reverse: file or tibble
  if (!is.character(reverse)){
    # Ensure required columns exist
    required_cols_rev <- c("Header", "Sequence", "Quality")
    if (!all(required_cols_rev %in% colnames(reverse))) {
      stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_reverse_file <- tempfile(pattern = "reverse_temp_", fileext = ".fq")
    microseq::writeFastq(reverse, temp_reverse_file)
    reverse_file <- temp_reverse_file
    temp_files <- c(temp_files, temp_reverse_file)

    # Capture original name for statistics table later
    reverse_name <- as.character(substitute(reverse))
  } else {
    reverse_file <- reverse

    # Capture original name for statistics table later
    reverse_name <- basename(reverse)
  }

  # Check if input files exists
  if (!file.exists(fastq_file)) stop("Cannot find input FASTQ file: ", fastq_file)
  if (!is.null(reverse) && !file.exists(reverse_file)) stop("Cannot find reverse FASTQ file: ", reverse_file)

  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)
  reverse_file <- normalizePath(reverse_file)

  # Determine output file
  if (is.null(fastqout)) {
    message("No filename for output file. No output file will be created.")
    outfile <- tempfile(pattern = "merged", fileext = ".fq")
    temp_files <- c(temp_files, outfile)
  } else {
    message("Writing merged sequences to file:", fastqout)
    outfile <- fastqout
  }

  # Build argument string for command line
  args <- c("--fastq_mergepairs", fastq_file,
            "--reverse", reverse_file,
            "--threads", threads,
            "--fastqout", outfile)

  # Add log file if specified by user
  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Read output into fastq object (tbl)
  merged_fastq <- microseq::readFastq(outfile)

  if (!is.null(log_file)) {
    vsearch_output <- readLines(log_file)
  }

  # Output statistics in table
  statistics <- parse_merge_statistics(vsearch_output, fastq_input_name, reverse_name)

  # Add statistics as attribute to merging table
  attr(merged_fastq, "statistics") <- statistics

  # Cleanup temporary files
  if (length(temp_files) > 0) {
    on.exit(
      file.remove(temp_files),
      add = TRUE)
  }

  return(merged_fastq)
}

#' Parse statistics from output text in stdout from read merging to tibble
#'
#' @param output string of output from running vs_fastq_mergepairs
#' @param R1_file name of file with forward reads that was used in merging
#' @param R2_file name of file with reverse reads that was used in merging
#'
#' @return table with merging metrics
#' @noRd
parse_merge_statistics <- function(output, R1_file, R2_file) {

  # Extract values from output
  pairs <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Pairs$"), "\\d+"))
  merged <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Merged"), "\\d+"))
  too_many_diff <- as.numeric(stringr::str_extract(stringr::str_subset(output, "too many differences"), "\\d+"))
  alignment_low <- as.numeric(stringr::str_extract(stringr::str_subset(output, "alignment score too low"), "\\d+"))
  mean_frag_length <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Mean fragment length"), "\\d+\\.\\d+"))
  stddev_frag_length <- as.numeric(stringr::str_extract(stringr::str_subset(output, "Standard deviation of fragment length"), "\\d+\\.\\d+"))

  # Create table
  result_table <- data.frame(
    Tot_num_pairs = pairs,
    Merged = merged,
    Too_Many_Differences = too_many_diff,
    Low_Alignment_Score_or_score_drop_too_high = alignment_low,
    Mean_Fragment_Length = mean_frag_length,
    StdDev_Fragment_Length = stddev_frag_length,
    R1 = R1_file,
    R2 = R2_file
  )

  return(result_table)
}
