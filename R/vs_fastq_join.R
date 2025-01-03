#' Join paired-end sequence reads
#'
#' @description Joins paired-end sequence reads into one sequence with a gap between them.
#'
#' @param fastq_input A FASTQ file path or a FASTQ object containing (forward) reads. See Details.
#' @param reverse A FASTQ file path or a FASTQ object containing (reverse) reads See Details.
#' @param output_format Desired output format of file or tibble: \code{"fasta"} or \code{"fastq"} (default).
#' @param fastaout Name of the FASTA output file with the joined reads. If \code{NULL} (default) no output will be written to file. See Details.
#' @param fastqout Name of the FASTQ output file with the joined reads. If \code{NULL} (default) no output will be written to file. See Details.
#' @param join_padgap The padding sequence to use in the gap between the sequences. Defaults to \code{"NNNNNNNN"}.
#' @param join_padgapq The quality of the padding sequence. Defaults to \code{"IIIIIIII"}, corresponding to a base quality score of 40 (a very high quality score with error probability \code{0.0001}).
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param log_file Name of the log file to capture messages from \code{vsearch}. If \code{NULL}, no log file is created. Defaults to \code{NULL}.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details The read-pairs in the input FASTQ-files (\code{fastq_input} and \code{reverse}) are joined into one sequence by adding a gap between them with a padding sequence, using \code{vsearch}.
#' The resulting sequences consist of the forward read, the padding sequence and the reverse complement of the reverse read.
#'
#' \code{fastq_input} and \code{reverse} can either be FASTQ files or FASTQ objects. FASTQ objects are tibbles that contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#' Forward and reverse reads must appear in the same order and total number in both files.
#'
#' If \code{fastaout} or \code{fastqout} is specified, the joined reads are output to this file in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}) the results are returned as a FASTA or FASTQ object, and no output is written to file. \code{output_format} has to match the desired output files/objects.
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will be discarded. By default, \code{minlen} is set to 0, which means that no sequences are removed.
#' However, using the default value may allow empty sequences to remain in the results.
#'
#' If \code{log_file} is specified, the messages and joining statistics are output to this file. If unspecified (\code{NULL}) no log file is written.
#'
#' @return Tibble or \code{NULL}.
#'
#' If output files are not specified, a tibble containing the joined reads in the format specified by \code{output_format} is returned. If an output file is specified, results are written to file and nothing is returned.
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
#' # Execute joining, with tibble as output
#' join_seqs <- vs_fastq_join(fastq_input = fastq_input,
#'                            reverse = reverse,
#'                            output_format = output_format)
#'
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_join <- function(fastq_input,
                          reverse,
                          output_format = "fastq",
                          fastaout = NULL,
                          fastqout = NULL,
                          join_padgap = "NNNNNNNN",
                          join_padgapq = "IIIIIIII",
                          fasta_width = 0,
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

  } else {
    fastq_file <- fastq_input
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

  } else {
    reverse_file <- reverse
  }

  # Check if input files exists
  if (!file.exists(fastq_file)) stop("Cannot find input FASTQ file: ", fastq_file)
  if (!is.null(reverse) && !file.exists(reverse_file)) stop("Cannot find reverse FASTQ file: ", reverse_file)

  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)
  reverse_file <- normalizePath(reverse_file)

  # Determine output file
  if (output_format == "fasta") {

    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "joined_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (output_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "joined_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }

  # Build argument string for command line
  args <- c("--fastq_join", fastq_file,
            "--reverse", reverse_file,
            "--join_padgap", join_padgap,
            "--join_padgapq", join_padgapq,
            "--threads", threads)

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

    # Create results tibble
    if (output_format == "fastq") {
      joined_seqs <- microseq::readFastq(outfile_fastq)
    } else if (output_format == "fasta") {
      joined_seqs <- microseq::readFasta(outfile_fasta)
    }
  }

  # Return results
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {
    return(joined_seqs)
  } else {
    return(invisible(NULL))
  }
}
