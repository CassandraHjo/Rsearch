#' Merge paired-end sequence reads
#'
#' @description Merges paired-end sequence reads with overlapping regions into one sequence.
#'
#' @param fastq_input A FASTQ file path or a FASTQ object containing (forward) reads. See Details.
#' @param reverse A FASTQ file path or a FASTQ object containing (reverse) reads See Details.
#' @param output_format Desired output format of file or tibble: \code{"fasta"} or \code{"fastq"} (default).
#' @param fastaout Name of the FASTA output file with the merged reads. If \code{NULL} (default) no output will be written to file. See Details.
#' @param fastqout Name of the FASTQ output file with the merged reads. If \code{NULL} (default) no output will be written to file. See Details.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}. See Details.
#' @param log_file Name of the log file to capture messages from \code{vsearch}. If \code{NULL}, no log file is created. Defaults to \code{NULL}.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details The read-pairs in the input FASTQ-files (\code{fastq_input} and
#' \code{reverse}) are merged if they have sufficient overlap, using \code{vsearch}.
#'
#' \code{fastq_input} and \code{reverse} can either be FASTQ files or FASTQ
#' objects. FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#' Forward and reverse reads must appear in the same order and total number in
#' both files.
#'
#' If \code{fastaout} or \code{fastqout} is specified, the merged reads are
#' output to this file in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}) the results are returned as a FASTA or FASTQ
#' object, and no output is written to file. \code{output_format} has to match
#' the desired output files/objects.
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on
#' lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will
#' be discarded. By default, \code{minlen} is set to 0, which means that no
#' sequences are removed.
#' However, using the default value may allow empty sequences to remain in
#' the results.
#'
#' If \code{log_file} is specified, the messages and merging statistics are output to this file. If unspecified (\code{NULL}) no log file is written. If \code{fastqout} is specified, then \code{log_file} needs to be specified in order to get the merging statistics from \code{vsearch}.
#'
#' @return Tibble or \code{NULL}.
#'
#' If output files are not specified, a tibble containing the merged reads in the
#' format specified by \code{output_format} is returned. If an output file is
#' specified, results are written to file and nothing is returned.
#'
#' When a FASTA/FASTQ object is returned, the statistics from the merging,
#' \code{statistics}, is an attribute of the merging tibble (\code{merged_seqs}).
#' This tibble contains merging statistics, including number of pairs, number
#' of merged pairs, and length metrics. The statistics can be accessed by
#' running \code{attributes(merged_seqs)$statistics} or
#' \code{attr(merged_seqs, "statistics")}.
#'
#' @examples
#' \dontrun{
#' # Read example FASTQ files
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"), "R2_sample1_small.fq")
#'
#' # Define other arguments
#' output_format <- "fastq"
#'
#' # Execute merging, with tibble as output
#' merge_seqs <- vs_fastq_mergepairs(fastq_input = fastq_input,
#'                                   reverse = reverse,
#'                                   output_format = output_format)
#'
#' # Extract merging statistics
#' statistics <- attr(merge_seqs, "statistics")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_mergepairs <- function(fastq_input,
                                reverse,
                                output_format = "fastq",
                                fastaout = NULL,
                                fastqout = NULL,
                                fasta_width = 0,
                                minlen = 0,
                                log_file = NULL,
                                threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
  }

  # If output_format is "fasta", fastqout can not be defined
  if (output_format == "fasta") {
    if (!is.null(fastqout)) {
      stop("When output_format is defined as 'fasta', 'fastqout' cannot be used. Use 'fastaout' instead.")
    }
  }

  # If output_format is "fastq", fastaout can not be defined
  if (output_format == "fastq") {
    if (!is.null(fastaout)) {
      stop("When output_format is defined as 'fastq', 'fastaout' cannot be used. Use 'fastqout' instead.")
    }
  }

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0) {
      file.remove(temp_files)
    }
  }, add = TRUE)

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
  fastq_file <- normalizePath(fastq_file) |>
    shQuote()
  reverse_file <- normalizePath(reverse_file) |>
    shQuote()

  # Determine output file
  if (output_format == "fasta") {

    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "merged_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (output_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "merged_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }


  # Build argument string for command line
  args <- c("--fastq_mergepairs", fastq_file,
            "--reverse", reverse_file,
            "--threads", threads,
            "--fastq_minlen", minlen)

  # Add output files based on output_format
  if (output_format == "fastq") {
    args <- c(args, "--fastqout", outfile_fastq)
  } else if (output_format == "fasta") {
    args <- c(args, "--fastaout", outfile_fasta, "--fasta_width", fasta_width)
  }

  # Add log file if specified
  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Handle output if output file is NULL
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    # Extract statistics
    statistics <- parse_merge_statistics(vsearch_output, fastq_input_name, reverse_name)

    # Create results tibble
    if (output_format == "fastq") {
      merged_seqs <- microseq::readFastq(outfile_fastq)
    } else if (output_format == "fasta") {
      merged_seqs <- microseq::readFasta(outfile_fasta)
    }

    # Add statistics as attribute to merging table
    attr(merged_seqs, "statistics") <- statistics
  }

  # Return results
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {
    return(merged_seqs)
  } else {
    return(invisible(NULL))
  }
}

#' Parse merging statistics from string to tibble
#'
#' @description This function transforms the output from \code{vsearch} when
#' running \code{vs_fastq_mergepairs()} into a tibble.
#' The most important statistics are included in the tibble such as number of
#' read pairs, merged reads, reasons that reads were not merged, and mean and
#' standard deviation of read lengths.
#'
#' @param output A string of output from merging reads with \code{vsearch}.
#' @param fastq_input The name of the file/object with forward (R1) reads that was used in the merging.
#' @param reverse The name of the file/object with reverse (R2) reads that was used in the merging.
#'
#' @return A tibble with merging metrics, including number of read pairs,
#' merged reads, reasons that reads were not merged, and mean and standard
#' deviation of read lengths.
#'
#' @noRd
#'
parse_merge_statistics <- function(output, fastq_input, reverse) {

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
    R1 = fastq_input,
    R2 = reverse
  )

  return(result_table)
}
