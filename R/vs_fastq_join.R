#' Join paired-end sequence reads
#'
#' @description \code{vs_fastq_join} joins paired-end sequence reads into a
#' single sequence with a specified gap between them using \code{VSEARCH}.
#'
#' @param fastq_input (Required). A FASTQ file path, a FASTQ tibble object
#' (forward reads), or a paired-end tibble of class \code{"pe_df"}. See
#' \emph{Details}.
#' @param reverse (Optional). A FASTQ file path or a FASTQ tibble object
#' (reverse reads). Optional if \code{fastq_input} is a \code{"pe_df"} object.
#' See \emph{Details}.
#' @param output_format (Optional). Desired output format of the file or tibble:
#' \code{"fasta"} or \code{"fastq"} (default).
#' @param fastaout (Optional). Name of the FASTA output file with the joined
#' reads. If \code{NULL} (default), no output is written to a file. See
#' \emph{Details}.
#' @param fastqout (Optional). Name of the FASTQ output file with the joined
#' reads. If \code{NULL} (default), no output is written to a file. See
#' \emph{Details}.
#' @param join_padgap (Optional). Padding sequence to use in the gap between the
#' sequences. Defaults to \code{"NNNNNNNN"}.
#' @param join_padgapq (Optional). Quality of the padding sequence. Defaults to
#' \code{"IIIIIIII"}, corresponding to a base quality score of 40 (a very high
#' quality score with error probability \code{0.0001}).
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Only applies if the output file is in FASTA format. Defaults to
#' \code{0}, which eliminates wrapping.
#' @param log_file (Optional). Name of the log file to capture messages from
#' \code{VSEARCH}. If \code{NULL} (default), no log file is created.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Read pairs from the input FASTQ files (\code{fastq_input} and \code{reverse})
#' are joined into a single sequence by adding a gap with a specified padding
#' sequence. The resulting sequences consist of the forward read, the padding
#' sequence, and the reverse complement of the reverse read.
#'
#' \code{fastq_input} and \code{reverse} can either be file paths to FASTQ files
#' or FASTQ objects. FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}, see
#' \code{\link{readFastq}}. Forward and reverse reads must appear in the same
#' order and have the same total number of reads in both files.
#'
#' If \code{fastq_input} is an object of class \code{"pe_df"}, the reverse reads
#' are automatically extracted from its \code{"reverse"} attribute unless
#' explicitly provided via the \code{reverse} argument. This simplifies function
#' calls when using paired-end tibbles created by functions such as
#' \code{\link{fastx_synchronize}} or \code{\link{vs_fastx_trim_filt}}.
#'
#' If \code{fastaout} or \code{fastqout} is specified, the joined reads are
#' written to the respective file in either FASTA or FASTQ format.
#'
#' If both \code{fastaout} or \code{fastqout} are \code{NULL}, the results are
#' returned as a FASTA or FASTQ object, and no file is written.
#' \code{output_format} must match the desired output files/objects.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} is
#' discarded. By default, \code{minlen} is set to 0, which means that no
#' sequences are removed. However, using the default value may allow empty
#' sequences to remain in the results.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{fastaout} or \code{fastqout} is specified, the joined sequences are
#' written to the specified output file, and no tibble is returned.
#'
#' If \code{fastaout} or \code{fastqout} is \code{NULL}, a tibble containing the
#' joined reads in the format specified by \code{output_format} is returned.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "small_R1.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R2.fq")
#' output_format <- "fastq"
#'
#' # Execute joining and return a FASTQ tibble
#' join_seqs <- vs_fastq_join(fastq_input = fastq_input,
#'                            reverse = reverse,
#'                            output_format = output_format)
#'
#' # Execute joining and write joined sequences to file
#' vs_fastq_join(fastq_input = fastq_input,
#'               reverse = reverse,
#'               fastqout = "joined_sequences.fq",
#'               output_format = output_format)
#'
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_fastq_join vs_fasta_join vs_fastx_join fastq_join join
#'
#' @export
#'
vs_fastq_join <- function(fastq_input,
                          reverse = NULL,
                          output_format = "fastq",
                          fastaout = NULL,
                          fastqout = NULL,
                          join_padgap = "NNNNNNNN",
                          join_padgapq = "IIIIIIII",
                          fasta_width = 0,
                          log_file = NULL,
                          threads = 1,
                          vsearch_options = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from 'fasta' or 'fastq'.")
  }

  if (output_format == "fasta" && !is.null(fastqout)) {
    stop("When output_format is 'fasta', 'fastqout' cannot be used. Use 'fastaout' instead.")
  }

  if (output_format == "fastq" && !is.null(fastaout)) {
    stop("When output_format is 'fastq', 'fastaout' cannot be used. Use 'fastqout' instead.")
  }

  # Extract reverse from pe_df class if not explicitly provided
  if (is_pe_df(fastq_input) && is.null(reverse)) {
    reverse <- attr(fastq_input, "reverse")
    if (is.null(reverse)) {
      stop("fastq_input has class 'pe_df' but no 'reverse' attribute found.")
    }
  }

  # Early file existence checks
  if (is.character(fastq_input) && !file.exists(fastq_input)) {
    stop("Cannot find input FASTQ file: ", fastq_input)
  }

  if (is.character(reverse) && !file.exists(reverse)) {
    stop("Cannot find reverse FASTQ file: ", reverse)
  }

  # Collect temporary files for cleanup
  temp_files <- character()
  on.exit({
    if (length(temp_files) > 0 && all(file.exists(temp_files))) {
      file.remove(temp_files)
    }
  }, add = TRUE)

  # Process forward reads
  if (!is.character(fastq_input)) {
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_fastq_file <- tempfile("fastq_input_", fileext = ".fq")
    microseq::writeFastq(fastq_input, temp_fastq_file)
    temp_files <- c(temp_files, temp_fastq_file)
    fastq_file <- temp_fastq_file
  } else {
    fastq_file <- normalizePath(fastq_input)
  }

  # Process reverse reads
  if (is.null(reverse)) {
    stop("No reverse reads provided. Please supply reverse or use a 'pe_df' object.")
  }

  if (!is.character(reverse)) {
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(reverse))) {
      stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_reverse_file <- tempfile("reverse_input_", fileext = ".fq")
    microseq::writeFastq(reverse, temp_reverse_file)
    temp_files <- c(temp_files, temp_reverse_file)
    reverse_file <- temp_reverse_file
  } else {
    reverse_file <- normalizePath(reverse)
  }

  if (!file.exists(fastq_file)) stop("Cannot find input FASTQ file: ", fastq_file)
  if (!file.exists(reverse_file)) stop("Cannot find reverse FASTQ file: ", reverse_file)

  # Define output file paths
  if (output_format == "fastq") {
    outfile <- if (is.null(fastqout)) tempfile("joined_", fileext = ".fq") else fastqout
  } else {
    outfile <- if (is.null(fastaout)) tempfile("joined_", fileext = ".fa") else fastaout
  }
  if (is.null(fastqout) && is.null(fastaout)) temp_files <- c(temp_files, outfile)

  # Construct VSEARCH arguments
  args <- c("--fastq_join", shQuote(fastq_file),
            "--reverse", shQuote(reverse_file),
            "--join_padgap", join_padgap,
            "--join_padgapq", join_padgapq,
            "--threads", threads)

  if (output_format == "fastq") {
    args <- c(args, "--fastqout", outfile)
  } else {
    args <- c(args, "--fastaout", outfile, "--fasta_width", fasta_width)
  }

  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  check_vsearch_status(vsearch_output, args)

  # Read and return results if not written to disk
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {
    joined_seqs <- if (output_format == "fastq") {
      microseq::readFastq(outfile)
    } else {
      microseq::readFasta(outfile)
    }
    return(joined_seqs)
  } else {
    return(invisible(NULL))
  }
}

