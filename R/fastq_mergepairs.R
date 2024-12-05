# TODO: Skrive dokumentasjon
# TODO: Skrive tester

#' Merge read-pairs
#'
#' @description Merging read-pairs with overlapping regions.
#'
#' @param fastq_file a FASTQ-file with forward reads (R1).
#' @param reverse a FASTQ-file with reverse reads (R2).
#' @param fastqout name of the FASTQ-file with the output or NULL, see Details.
#' @param log_file name of the log file with messages from running vsearch or NULL, see Details.
#' @param threads number of computational threads to be used by vsearch.
#'
#' @details The read-pairs in the input FASTQ-files (\code{fastq_file} and \code{reverse})
#' are merged if they have sufficient overlap, using vsearch.
#'
#' If \code{fastqout} is specified, the merged reads are output to this file in FASTQ-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTQ-object, i.e. a tibble with
#' columns \code{Header}, \code{Sequence} and \code{Quality}.
#'
#' If \code{log_file} is specified, the messages are output to this file.
#' If unspecified (\code{NULL}) no log file is written.
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
vs_fastq_mergepairs <- function(fastq_file,
                                reverse,
                                fastqout = NULL,
                                log_file = NULL,
                                threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Check is input files exist at given paths
  if (!file.exists(fastq_file)) stop("Cannot find input file: ", fastq_file)
  if (!file.exists(reverse)) stop("Cannot find reverse file: ", reverse)

  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)
  reverse <- normalizePath(reverse)

  # Create empty vector for collecting temporary files
  temp_files <- c()

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
            "--reverse", reverse,
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
  statistics <- parse_merge_statistics(vsearch_output, fastq_file, reverse)

  # Add statistics as attribute to merging table
  attr(merged_fastq, "statistics") <- statistics

  # # Remove temp file if necessary
  # if (is.null(fastqout)) {
  #   file.remove(outfile)
  # }

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

  if (!is.character(R1_file)) {
    R1_source <- deparse(substitute(R1_file))
  } else {
    R1_source <- basename(R1_file)
  }

  if (!is.character(R2_file)) {
    R2_source <- deparse(substitute(R2_source))
  } else {
    R2_source <- basename(R2_source)
  }

  # Create table
  result_table <- data.frame(
    Tot_num_pairs = pairs,
    Merged = merged,
    Too_Many_Differences = too_many_diff,
    Low_Alignment_Score_or_score_drop_too_high = alignment_low,
    Mean_Fragment_Length = mean_frag_length,
    StdDev_Fragment_Length = stddev_frag_length,
    R1_source = R1_source,
    R2_source = R2_source
  )

  return(result_table)
}
