#' Merge paired-end sequence reads
#'
#' @description \code{vs_fastq_mergepairs} merges paired-end sequence reads with
#' overlapping regions into one sequence using \code{VSEARCH}.
#'
#' @param fastq_input (Required). A FASTQ file path, a FASTQ tibble (forward
#' reads), or a paired-end tibble of class \code{"pe_df"}. See \emph{Details}.
#' @param reverse (Optional). A FASTQ file path or a FASTQ tibble (reverse
#' reads). Optional if \code{fastq_input} is a \code{"pe_df"} object. See
#' \emph{Details}.
#' @param output_format (Optional). Desired output format of file or tibble:
#' \code{"fasta"} (default) or \code{"fastq"}.
#' @param fastaout (Optional). Name of the FASTA output file with the merged
#' reads. If \code{NULL} (default), no output is written to file. See
#' \emph{Details}.
#' @param fastqout (Optional). Name of the FASTQ output file with the merged
#' reads. If \code{NULL} (default) no output is written to file. See
#' \emph{Details}.
#' @param minovlen (Optional). Minimum overlap between the merged reads. Must be
#' at least 5. Defaults to \code{10}.
#' @param minlen (Optional). Minimum number of bases a sequence must have to be
#' retained. Defaults to \code{0}. See \emph{Details}.
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Only applies if the output file is in FASTA format. Defaults to
#' \code{0}, which eliminates wrapping.
#' @param sample (Optional). Add the given sample identifier string to sequence
#' headers. For instance, if the given string is "ABC", the text ";sample=ABC"
#' will be added to the header. If \code{NULL} (default), no identifier is added.
#' @param log_file (Optional). Name of the log file to capture messages from
#' \code{VSEARCH}. If \code{NULL} (default), no log file is created.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details
#' Read pairs from the input FASTQ files (\code{fastq_input} and \code{reverse})
#' are merged into a single sequence by overlapping regions. The resulting
#' sequences consist of the merged forward and reverse reads with the specified
#' minimum overlap.
#'
#' \code{fastq_input} and \code{reverse} can either be file paths to FASTQ files
#' or FASTQ objects. FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}, see
#' \code{\link[microseq]{readFastq}}. Forward and reverse reads must appear in
#' the same order and have the same total number of reads in both files.
#'
#' If \code{fastq_input} is an object of class \code{"pe_df"}, the reverse reads
#' are automatically extracted from its \code{"reverse"} attribute unless
#' explicitly provided via the \code{reverse} argument. This allows streamlined
#' input handling for paired-end tibbles created by
#' \code{\link{fastx_synchronize}} or \code{\link{vs_fastx_trim_filt}}.
#'
#' If \code{fastaout} or \code{fastqout} is specified, the merged reads are
#' written to the respective file in either FASTA or FASTQ format.
#'
#' If both \code{fastaout} or \code{fastqout} are \code{NULL}, the results are
#' returned as a FASTA or FASTQ object, and no file is written.
#'
#' \code{output_format} has to match the desired output files/objects.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will
#' be discarded. Default \code{minlen} is 0, meaning no sequences are removed.
#' However, using the default value may allow empty sequences to remain in
#' the results.
#'
#' If \code{log_file} is \code{NULL} and \code{fastqout} or \code{fastaout} is
#' specified, merging statistics from \code{VSEARCH} will not be captured.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{fastaout} or \code{fastqout} is specified , the merged sequences are
#' written to the specified output file, and no tibble is returned.
#'
#' If \code{fastaout} or \code{fastqout} is \code{NULL}, a tibble containing the
#' merged reads in the format specified by \code{output_format} is returned.
#'
#' The \code{"statistics"} attribute of the returned tibble (when
#' \code{fastaout} or \code{fastqout} is \code{NULL}) is a tibble with the
#' following columns:
#' \itemize{
#'   \item \code{Tot_num_pairs}: Total number of read pairs before merging.
#'   \item \code{Merged}: Number of read pairs that merged.
#'   \item \code{Mean_Read_Length_before_merging}: Mean read length before
#'   merging (R1 and R2).
#'   \item \code{Mean_Read_Length_after_merging}: Mean read length after
#'   merging.
#'   \item \code{StdDev_Read_Length}: Standard deviation of read length
#'   after merging.
#'   \item \code{R1}: Name of the file/object with forward (R1) reads used in
#'   the merging.
#'   \item \code{R2}: Name of the file/object with reverse (R2) reads used in
#'   the merging.
#' }
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
#' # Merge sequences and return a FASTQ tibble
#' merge_seqs <- vs_fastq_mergepairs(fastq_input = fastq_input,
#'                                   reverse = reverse,
#'                                   output_format = output_format)
#'
#' # Extract merging statistics
#' statistics <- attr(merge_seqs, "statistics")
#'
#' # Merge sequences and write sequences to a FASTQ file
#' vs_fastq_mergepairs(fastq_input = fastq_input,
#'                     reverse = reverse,
#'                     output_format = output_format,
#'                     fastqout = "merged_sequences.fq")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_fastq_mergepairs vs_fastx_mergepairs vs_fasta_mergepairs
#' vs_mergepairs fastq_mergepairs mergepairs
#'
#' @export
#'
vs_fastq_mergepairs <- function(fastq_input,
                                reverse = NULL,
                                output_format = "fasta",
                                fastaout = NULL,
                                fastqout = NULL,
                                minovlen = 10,
                                minlen = 0,
                                fasta_width = 0,
                                sample = NULL,
                                log_file = NULL,
                                threads = 1,
                                vsearch_options = NULL,
                                tmpdir = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
  }

  if (output_format == "fasta" && !is.null(fastqout)) {
    stop("When output_format is 'fasta', 'fastqout' cannot be used. Use 'fastaout' instead.")
  }

  if (output_format == "fastq" && !is.null(fastaout)) {
    stop("When output_format is 'fastq', 'fastaout' cannot be used. Use 'fastqout' instead.")
  }

  # Extract reverse if fastq_input is a pe_df object
  if (is_pe_df(fastq_input) && is.null(reverse)) {
    reverse <- attr(fastq_input, "reverse")
    if (is.null(reverse)) {
      stop("fastq_input has class 'pe_df' but no 'reverse' attribute found.")
    }
  }

  # Early checks for file existence
  if (is.character(fastq_input) && !file.exists(fastq_input)) {
    stop("Cannot find input FASTQ file: ", fastq_input)
  }
  if (is.character(reverse) && !file.exists(reverse)) {
    stop("Cannot find reverse FASTQ file: ", reverse)
  }

  temp_files <- c()
  on.exit({
    if (length(temp_files) > 0) file.remove(temp_files[file.exists(temp_files)])
  }, add = TRUE)

  # Forward read handling
  if (!is.character(fastq_input)) {
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_fastq_file <- tempfile("fastq_input_",
                                tmpdir = tmpdir,
                                fileext = ".fq")
    microseq::writeFastq(fastq_input, temp_fastq_file)
    temp_files <- c(temp_files, temp_fastq_file)
    fastq_file <- temp_fastq_file
    fastq_input_name <- deparse(substitute(fastq_input))[1]
  } else {
    fastq_file <- normalizePath(fastq_input)
    fastq_input_name <- basename(fastq_input)
  }

  # Reverse read handling
  if (is.null(reverse)) {
    stop("No reverse reads provided. Please supply reverse or use a 'pe_df' object.")
  }

  if (!is.character(reverse)) {
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(reverse))) {
      stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_reverse_file <- tempfile("reverse_input_",
                                  tmpdir = tmpdir,
                                  fileext = ".fq")
    microseq::writeFastq(reverse, temp_reverse_file)
    temp_files <- c(temp_files, temp_reverse_file)
    reverse_file <- temp_reverse_file
    reverse_name <- deparse(substitute(reverse))[1]
  } else {
    reverse_file <- normalizePath(reverse)
    reverse_name <- basename(reverse)
  }

  # Output file setup
  if (output_format == "fasta") {
    outfile <- if (is.null(fastaout)) tempfile("merged_", tmpdir = tmpdir, fileext = ".fa") else fastaout
  } else {
    outfile <- if (is.null(fastqout)) tempfile("merged_", tmpdir = tmpdir, fileext = ".fq") else fastqout
  }
  if (is.null(fastaout) && is.null(fastqout)) temp_files <- c(temp_files, outfile)

  # Build VSEARCH args
  args <- c("--fastq_mergepairs", shQuote(fastq_file),
            "--reverse", shQuote(reverse_file),
            "--fastq_minovlen", minovlen,
            "--threads", threads,
            "--fastq_minlen", minlen)

  if (output_format == "fastq") {
    args <- c(args, "--fastqout", outfile)
  } else {
    args <- c(args, "--fastaout", outfile, "--fasta_width", fasta_width)
  }

  if (!is.null(sample)) {
    args <- c(args, "--sample", sample)
  }

  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  check_vsearch_status(vsearch_output, args)

  # Read results if not written to file
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    if (file.size(outfile) == 0) {
      stop(paste("Output file is empty. No reads were merged.",
                 "Arguments used in merging: ", paste(args, collapse = " ")))
    }

    merged_seqs <- if (output_format == "fastq") {
      microseq::readFastq(outfile)
    } else {
      microseq::readFasta(outfile)
    }

    stats.tbl <- calculate_merge_statistics(fastq_file,
                                            reverse_file,
                                            merged_seqs,
                                            fastq_input_name,
                                            reverse_name)

    attr(merged_seqs, "statistics") <- stats.tbl
    return(merged_seqs)
  }

  return(invisible(NULL))
}

#' Calculate merging statistics
#'
#' @description Calculates important merging statistics after running
#' \code{vs_fastq_mergepairs()},like number of read pairs, merged reads, and
#' mean and standard deviation of read lengths.
#'
#' @param fastq_file File path to FASTQ file containing the forward reads (R1)
#' used as input for the merging.
#' @param reverse_file File path to FASTQ file containing the reverse reads (R2)
#' used as input for the merging.
#' @param merged_seqs Tibble containing the merged sequences resulting from the
#' merging process.
#' @param fastq_file_name Name of the file/object with forward (R1) reads used
#' in the merging.
#' @param reverse_file_name The name of the file/object with reverse (R2) reads
#' used in the merging.
#'
#' @return A tibble with merging statistics, including:
#' \itemize{
#'   \item \code{Tot_num_pairs}: Total number of read pairs before merging.
#'   \item \code{Merged}: Number of read pairs that merged.
#'   \item \code{Mean_Read_Length_before_merging}: Mean read length before
#'   merging (R1 and R2).
#'   \item \code{Mean_Read_Length_after_merging}: Mean read length after
#'   merging.
#'   \item \code{StdDev_Read_Length}: Standard deviation of read length
#'   after merging.
#'   \item \code{R1}: Name of the file/object with forward (R1) reads used in
#'   the merging.
#'   \item \code{R2}: Name of the file/object with reverse (R2) reads used in
#'   the merging.
#' }
#'
#' @return A tibble with merging statistics.
#'
#' @noRd
#'
calculate_merge_statistics <- function(fastq_file,
                                       reverse_file,
                                       merged_seqs,
                                       fastq_file_name,
                                       reverse_file_name) {

  # Read forward and reverse FASTQ files
  r1 <- microseq::readFastq(fastq_file)
  r2 <- microseq::readFastq(reverse_file)

  # Calculate statistics
  pairs <- nrow(r1)
  merged <- nrow(merged_seqs)
  mean_length_before <- round(((mean(nchar(r1$Sequence)) + mean(nchar(r2$Sequence))) / 2), 2)
  mean_length_merged_reads <- round(mean(nchar(merged_seqs$Sequence)), 2)
  sd_read_length_merged_reads <- round(stats::sd(nchar(merged_seqs$Sequence)), 2)


  # Create statistics table
  result_table <- tibble::tibble(
    Tot_num_pairs = pairs,
    Merged = merged,
    Mean_Read_Length_before_merging = mean_length_before,
    Mean_Read_Length_after_merging = mean_length_merged_reads,
    StdDev_Read_Length = sd_read_length_merged_reads,
    R1 = fastq_file_name,
    R2 = reverse_file_name
  )

  return(result_table)
}
