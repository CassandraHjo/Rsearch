#' Trim and/or filter sequences in FASTA/FASTQ format
#'
#' @description \code{vs_fastx_trim_filt} trims and/or filters FASTA/FASTQ
#' sequences using \code{VSEARCH}. This function processes both forward and
#' reverse reads (if provided) and allows for various filtering criteria based
#' on sequence quality, length, abundance, and more.
#'
#' @param fastx_input (Required). A FASTA/FASTQ file path or FASTA/FASTQ object
#' containing (forward) reads. See \emph{Details}.
#' @param reverse (Optional). A FASTA/FASTQ file path or object containing
#' reverse reads. If \code{fastx_input} is a \code{"pe_df"} object and
#' \code{reverse} is not provided, the reverse reads will be extracted from its
#' \code{"reverse"} attribute.
#' @param output_format (Optional). Desired output format of file or tibble:
#' \code{"fasta"} or \code{"fastq"} (default). If \code{fastx_input} is a FASTA
#' file path or a FASTA object, \code{output_format} cannot be \code{"fastq"}.
#' @param fastaout (Optional). Name of the FASTA output file for the sequences
#' given in \code{fastx_input}. If \code{NULL} (default), no FASTA sequences are
#' written to file. See \emph{Details}.
#' @param fastqout (Optional). Name of the FASTQ output file for the sequences
#' given in \code{fastx_input}. If \code{NULL} (default), no FASTQ sequences are
#' written to file. See \emph{Details}.
#' @param fastaout_rev (Optional). Name of the FASTA output file for the reverse
#' sequences. If \code{NULL} (default), no FASTA sequences are written to file.
#' See \emph{Details}.
#' @param fastqout_rev (Optional). Name of the FASTQ output file for the reverse
#' sequences. If \code{NULL} (default), no FASTQ sequences are written to file.
#' See \emph{Details}.
#' @param trunclen (Optional). Truncate sequences to the specified length.
#' Shorter sequences are discarded. If \code{NULL} (default), the trimming is
#' not applied.
#' @param truncqual (Optional). Truncate sequences starting from the first base
#' with a quality score of the specified value or lower. Defaults to \code{1}.
#' @param truncee (Optional). Truncate sequences so that their total expected
#' error does not exceed the specified value. If \code{NULL} (default), the
#' trimming is not applied.
#' @param truncee_rate (Optional). Truncate sequences so that their average
#' expected error per base is not higher than the specified value. The
#' truncation will happen at first occurrence. The average expected error per
#' base is calculated as the total expected number of errors divided by the
#' length of the sequence after truncation. If \code{NULL} (default), the
#' trimming is not applied.
#' @param stripright (Optional). Number of bases stripped from the right end of
#' the reads. Defaults to \code{0}.
#' @param stripleft (Optional). Number of bases stripped from the left end of
#' the reads. Defaults to \code{0}.
#' @param maxee_rate (Optional). Threshold for average expected error. Numeric
#' value ranging form \code{0.0} to \code{1.0}. Defaults to \code{0.01}. See
#' \emph{Details}.
#' @param minlen (Optional). Minimum number of bases a sequence must have to be
#' retained. Defaults to \code{0}. See \emph{Details}.
#' @param maxlen (Optional). Maximum number of bases a sequences can have to be
#' retained. If \code{NULL} (default), the filter is not applied.
#' @param maxns (Optional). Maximum number of N's for a given sequence.
#' Sequences with more N's than the specified number are discarded. Defaults to
#' \code{0}.
#' @param minsize (Optional). Minimum abundance for a given sequence. Sequences
#' with lower abundance are discarded. If \code{NULL} (default), the filter is
#' not applied.
#' @param maxsize (Optional). Maximum abundance for a given sequence. Sequences
#' with higher abundance are discarded. If \code{NULL} (default), the filter is
#' not applied.
#' @param minqual (Optional). Minimum base quality for a read to be retained. A
#' read is discarded if it contains bases with a quality score below the given
#' value. Defaults to \code{0}, meaning no reads are discarded.
#' @param relabel (Optional). Relabel sequences using the given prefix and a
#' ticker to construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 (Optional). If \code{TRUE} (default), relabel sequences
#' using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Defaults to \code{0}, which eliminates wrapping.
#' @param sample (Optional). Add the given sample identifier string to sequence
#' headers. For instance, if the given string is "ABC", the text ";sample=ABC"
#' will be added to the header. If \code{NULL} (default), no identifier is added.
#' @param stats (Optional). If \code{TRUE} (default), a tibble with statistics
#' about the filtering is added as an attribute of the returned tibble. If
#' \code{FALSE}, no statistics are added.
#' @param log_file (Optional). Name of the log file to capture messages from
#' \code{VSEARCH}. If \code{NULL} (default), no log file is created.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Reads from the input files/objects (\code{fastx_input} and \code{reverse})
#' are trimmed and/or filtered based on the specified criteria using
#' \code{VSEARCH}.
#'
#' \code{fastx_input} and \code{reverse} can either be file paths to FASTA/FASTQ
#' files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
#' columns \code{Header} and \code{Sequence}, see \code{\link{readFasta}}. FASTQ
#' objects are tibbles that contain the columns \code{Header}, \code{Sequence},
#' and \code{Quality}, see \code{\link{readFastq}}.
#'
#' If \code{fastx_input} is an object of class \code{"pe_df"}, the reverse reads
#' are automatically extracted from its \code{"reverse"} attribute unless
#' explicitly provided via the \code{reverse} argument.
#'
#' If \code{reverse} is provided, it is processed alongside \code{fastx_input}
#' using the same trimming/filtering criteria.
#'
#' Note that if you want to trim/filter the forward and reverse reads
#' differently, you must pass them separately to this function, get two result
#' files/objects, and then use \code{\link{fastx_synchronize}} to synchronize
#' the read pairs again.
#'
#' If \code{fastaout} and \code{fastaout_rev} or \code{fastqout} and
#' \code{fastqout_rev} are specified, trimmed and/or filtered sequences are
#' written to these files in the specified format.
#'
#' If output files are \code{NULL}, results are returned as a tibbles. When
#' returning tibbles, the reverse sequences (if provided) are attached as an
#' attribute named \code{"reverse"}.
#'
#' When reverse reads are returned as an attribute, the primary tibble is also
#' assigned the S3 class \code{"pe_df"} to indicate that it represents
#' paired-end data. This class tag can be used by downstream tools to recognize
#' paired-end tibbles.
#'
#' Note that certain options are not compatible with both file formats. For
#' instance, options that trim or filter sequences based on quality scores are
#' unavailable when the input is of type \code{"fasta"}. Visit the
#' \code{VSEARCH}
#' \href{https://github.com/torognes/vsearch?tab=readme-ov-file#getting-help}{documentation}
#' for more details.
#'
#' Sequences with an average expected error greater than the specified
#' \code{maxee_rate} are discarded. For a given sequence, the average expected
#' error is the sum of error probabilities for all the positions in the sequence,
#' divided by the length of the sequence.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will
#' be discarded. By default, \code{minlen} is set to 0, which means that no
#' sequences are removed. However, using the default value may allow empty
#' sequences to remain in the results.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If output files are specified, the results are written directly to the
#' specified output files, and no tibble is returned.
#'
#' If output files (\code{fastaout}/\code{fastqout} and
#' \code{fastaout_rev}/\code{fastqout_rev}) are \code{NULL}, a tibble containing
#' the trimmed and/or filtered reads from \code{fastx_input} in the format
#' specified by \code{output_format} is returned.
#'
#' If \code{reverse} is provided, a tibble containing the trimmed and/or
#' filtered reverse sequences is attached as an attribute, named
#' \code{"reverse"} to the returned table.
#'
#' When the reverse reads are present, the returned tibble is assigned the
#' class \code{"pe_df"}, identifying it as paired-end data.
#'
#' The \code{"statistics"} attribute of the returned tibble (when
#' output files are \code{NULL}) is a tibble with the
#' following columns:
#' \itemize{
#'   \item \code{Kept_Sequences}: Number of retained sequences.
#'   \item \code{Discarded_Sequences}: Number of discarded sequences.
#'   \item \code{fastx_source}: Name of the file/object with forward (R1) reads.
#'   \item \code{reverse_source}: (If \code{reverse} is specified) Name of the
#'   file/object with reverse (R2) reads.
#' }
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "small_R1.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small_R1.fq")
#' output_format <- "fastq"
#' maxee_rate <- 0.01
#' minlen <- 0
#'
#' # Trim/filter sequences and return a FASTQ tibble
#' filt_seqs <- vs_fastx_trim_filt(fastx_input = fastx_input,
#'                                 reverse = reverse,
#'                                 output_format = output_format,
#'                                 maxee_rate = maxee_rate,
#'                                 minlen = minlen)
#'
#' # Extract tibbles
#' R1_filt <- filt_seqs
#' R2_filt <- attr(filt_seqs, "reverse")
#'
#' # Extract filtering statistics
#' statistics <- attr(filt_seqs, "statistics")
#'
#' # Trim/filter sequences and write results to FASTQ files
#' vs_fastx_trim_filt(fastx_input = fastx_input,
#'                    reverse = reverse,
#'                    fastqout = "filt_R1.fq",
#'                    fastqout_rev = "filt_R2.fq",
#'                    output_format = output_format,
#'                    maxee_rate = maxee_rate,
#'                    minlen = minlen)
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_fastx_trim_filt vs_fastq_trim_filt vs_fasta_trim_filt
#' fastx_trim_filt trim_filt
#'
#' @export
#'
vs_fastx_trim_filt <- function(fastx_input,
                               reverse = NULL,
                               output_format = "fastq",
                               fastaout = NULL,
                               fastqout = NULL,
                               fastaout_rev = NULL,
                               fastqout_rev = NULL,
                               trunclen = NULL,
                               truncqual = 1,
                               truncee = NULL,
                               truncee_rate = NULL,
                               stripright = 0,
                               stripleft = 0,
                               maxee_rate = 0.01,
                               minlen = 0,
                               maxlen = NULL,
                               maxns = 0,
                               minsize = NULL,
                               maxsize = NULL,
                               minqual = 0,
                               relabel = NULL,
                               relabel_sha1 = FALSE,
                               fasta_width = 0,
                               sample = NULL,
                               stats = TRUE,
                               log_file = NULL,
                               threads = 1,
                               vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
  }

  # If output_format is "fastq", fastaout and fastaout_rev can not be defined
  if (output_format == "fastq") {
    if (!is.null(fastaout) || !is.null(fastaout_rev)) {
      stop("When output_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
    }
  }

  # If output_format is "fasta", fastqout and fastqout_rev can not be defined
  if (output_format == "fasta") {
    if (!is.null(fastqout) || !is.null(fastqout_rev)) {
      stop("When output_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
    }
  }

  # Handle case when fastx_input is a pe_df and reverse is NULL
  if (is_pe_df(fastx_input) && is.null(reverse)) {
    reverse <- attr(fastx_input, "reverse")
    if (is.null(reverse)) {
      stop("fastx_input has class 'pe_df' but no 'reverse' attribute found.")
    }
  }

  # If reverse is specified, ensure paired output parameters are both NULL or both character strings
  if (!is.null(reverse)) {
    if (output_format == "fasta") {
      # Check that both fastaout and fastaout_rev are NULL or both are character strings
      if ((is.null(fastaout) && !is.null(fastaout_rev)) ||
          (!is.null(fastaout) && is.null(fastaout_rev))) {
        stop("When 'reverse' is specified and output_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
      }
    }

    if (output_format == "fastq") {
      # Check that both fastqout and fastqout_rev are NULL or both are character strings
      if ((is.null(fastqout) && !is.null(fastqout_rev)) ||
          (!is.null(fastqout) && is.null(fastqout_rev))) {
        stop("When 'reverse' is specified and output_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
      }
    }
  }

  # Create empty vector for collecting temporary files
  temp_files <- character()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0 && is.character(temp_files)) {
      existing_files <- temp_files[file.exists(temp_files)]
      if (length(existing_files) > 0) {
        file.remove(existing_files)
      }
    }
  }, add = TRUE)

  # Handle input for primary sequences
  fasta_input_detected <- FALSE
  if (!is.character(fastx_input)){
    if ("Quality" %in% colnames(fastx_input)){

      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_file_primary <- tempfile(pattern = "primary_input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_primary)
      microseq::writeFastq(fastx_input, temp_file_primary)

      fastx_file <- temp_file_primary

      # Capture original name for statistics table later
      fastx_input_name <- as.character(substitute(fastx_input))
    } else {

      if (output_format == "fastq") {
        stop("Invalid output_format when input tibble is of type 'fasta'")
      }

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_primary <- tempfile(pattern = "primary_input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file_primary)
      microseq::writeFasta(fastx_input, temp_file_primary)

      fastx_file <- temp_file_primary

      # Capture original name for statistics table later
      fastx_input_name <- as.character(substitute(fastx_input))

      fasta_input_detected <- TRUE
    }
  } else {
    if (!file.exists(fastx_input)) stop("Cannot find input file: ", fastx_input)

    fastx_file <- fastx_input

    # Capture original name for statistics table later
    fastx_input_name <- basename(fastx_input)
  }

  # Disable FASTQ-only options for FASTA input
  if (!is.character(fastx_input)) {
    if (!"Quality" %in% colnames(fastx_input)) {
      fasta_input_detected <- TRUE
    }
  } else {
    # For character path input, check if input appears to be a FASTA file by reading the first line
    first_line <- tryCatch(readLines(fastx_input, n = 1), error = function(e) "")
    if (grepl("^>", first_line)) {
      fasta_input_detected <- TRUE
    }
  }

  if (isTRUE(fasta_input_detected)) {
    if (!is.null(maxee_rate)) {
      warning("maxee_rate is ignored for FASTA input and will be set to NULL")
      maxee_rate <- NULL
    }
    if (!is.null(truncqual)) {
      warning("truncqual is ignored for FASTA input and will be set to NULL")
      truncqual <- NULL
    }
    if (!is.null(truncee)) {
      warning("truncee is ignored for FASTA input and will be set to NULL")
      truncee <- NULL
    }
    if (!is.null(truncee_rate)) {
      warning("truncee_rate is ignored for FASTA input and will be set to NULL")
      truncee_rate <- NULL
    }
  }

  # Handle input for reverse sequences
  if (!is.null(reverse)){
    if (!is.character(reverse)){
      if ("Quality" %in% colnames(reverse)){

        # Validate tibble
        required_cols <- c("Header", "Sequence", "Quality")
        if (!all(required_cols %in% colnames(reverse))) {
          stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
        }

        temp_file_reverse <- tempfile(pattern = "reverse_input", fileext = ".fq")
        temp_files <- c(temp_files, temp_file_reverse)
        microseq::writeFastq(reverse, temp_file_reverse)

        reverse_file <- temp_file_reverse

        # Capture original name for statistics table later
        reverse_name <- as.character(substitute(reverse))

      } else {

        if (output_format == "fastq") {
          stop("Invalid output_format when input tibble is of type 'fasta'")
        }

        # Validate tibble
        required_cols <- c("Header", "Sequence")
        if (!all(required_cols %in% colnames(reverse))) {
          stop("Reverse FASTA object must contain columns: Header and Sequence")
        }

        temp_file_reverse <- tempfile(pattern = "reverse_input", fileext = ".fa")
        temp_files <- c(temp_files, temp_file_reverse)
        microseq::writeFasta(reverse, temp_file_reverse)

        reverse_file <- temp_file_reverse

        # Capture original name for statistics table later
        reverse_name <- as.character(substitute(reverse))
      }
    } else {
      if (!file.exists(reverse)) stop("Cannot find reverse file: ", reverse)

      reverse_file <- reverse

      # Capture original name for statistics table later
      reverse_name <- basename(reverse)
    }
  }

  # Handle output for primary sequences
  if (output_format == "fasta") {
    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "filtered_primary_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (output_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "filtered_primary_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }

  # Handle output for reverse sequences
  if (!is.null(reverse)) {
    if (output_format == "fasta") {
      if (is.null(fastaout_rev)) {
        outfile_fasta_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fa")
        temp_files <- c(temp_files, outfile_fasta_rev)
      } else {
        outfile_fasta_rev <- fastaout_rev
      }
    }

    if (output_format == "fastq") {
      if (is.null(fastqout_rev)) {
        outfile_fastq_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fq")
        temp_files <- c(temp_files, outfile_fastq_rev)
      } else {
        outfile_fastq_rev <- fastqout_rev
      }
    }
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  if (!is.null(reverse)) {
    reverse_file <- normalizePath(reverse_file)
  }

  # Build argument string for command line
  args <- c("--fastx_filter", shQuote(fastx_file),
            "--fastq_minlen", minlen,
            "--threads", threads)

  # Add reverse to arguments if provided
  if (!is.null(reverse)) {
    args <- c(args, "--reverse", shQuote(reverse_file))
  }

  # Add trimming and filtering arguments if provided
  if (!is.null(maxee_rate)) {
    args <- c(args, "--fastq_maxee_rate", maxee_rate)
  }

  if (!is.null(maxlen)) {
    args <- c(args, "--fastq_maxlen", maxlen)
  }

  if (!is.null(maxns)) {
    args <- c(args, "--fastq_maxns", maxns)
  }

  if (!is.null(maxsize)) {
    args <- c(args, "--maxsize", maxsize)
  }

  if (!is.null(minsize)) {
    args <- c(args, "--minsize", minsize)
  }

  if (!is.null(trunclen)) {
    args <- c(args, "--fastq_trunclen", trunclen)
  }

  if (!is.null(truncqual)) {
    args <- c(args, "--fastq_truncqual", truncqual)
  }

  if (!is.null(truncee)) {
    args <- c(args, "--fastq_truncee", truncee)
  }

  if (!is.null(truncee_rate)) {
    args <- c(args, "--fastq_truncee_rate", truncee_rate)
  }

  if (stripright > 0) {
    args <- c(args, "--fastq_stripright", stripright)
  }

  if (stripleft > 0) {
    args <- c(args, "--fastq_stripleft", stripleft)
  }

  if (minqual > 0) {
    args <- c(args, "--fastq_minqual", minqual)
  }

  if (relabel_sha1) {
    args <- c(args, "--relabel_sha1", "")
  }

  if (!is.null(relabel)) {
    args <- c(args, "--relabel", relabel)
  }

  # Add output files based on output_format
  if (output_format == "fastq") {
    args <- c(args, "--fastqout", outfile_fastq)
    if (!is.null(reverse)) {
      args <- c(args, "--fastqout_rev", outfile_fastq_rev)
    }
  } else if(output_format == "fasta") {
    args <- c(args, "--fastaout", outfile_fasta, "--fasta_width", fasta_width)
    if (!is.null(reverse)) {
      args <- c(args, "--fastaout_rev", outfile_fasta_rev)
    }
  }

  # Add sample identifier if specified
  if (!is.null(sample)) {
    args <- c(args, "--sample", sample)
  }

  # Add log file if specified
  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Handle output if output files are NULL
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    # Process primary sequences
    if (output_format == "fasta") {
      filt_seqs <- microseq::readFasta(outfile_fasta)
    } else if (output_format == "fastq") {
      filt_seqs <- microseq::readFastq(outfile_fastq)
    }

    # Process reverse sequences if provided
    if (!is.null(reverse)) {
      if (output_format == "fasta") {
        filt_reverse <- microseq::readFasta(outfile_fasta_rev)
      } else if (output_format == "fastq") {
        filt_reverse <- microseq::readFastq(outfile_fastq_rev)
      }
    }

    # Extract statistics
    if (isTRUE(stats)) {

      statistics <- calculate_trim_filt_statistics(
        input_file = fastx_file,
        output_file = if (output_format == "fastq") outfile_fastq else outfile_fasta,
        format = output_format,
        fastx = fastx_input_name,
        reverse = if (!is.null(reverse)) reverse_name else NULL,
        reverse_file = if (!is.null(reverse)) reverse_file else NULL,
        reverse_output_file = if (!is.null(reverse)) {
          if (output_format == "fastq") outfile_fastq_rev else outfile_fasta_rev
        } else NULL
      )

      attr(filt_seqs, "statistics") <- statistics
    }

    # Add reverse table as attribute and class tag to the primary table
    if (!is.null(reverse)) {
      attr(filt_seqs, "reverse") <- filt_reverse

      # Add class label
      class(filt_seqs) <- c("pe_df", class(filt_seqs))
    }
  }

  # Return results
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    return(filt_seqs)
  } else {
    return(invisible(NULL))
  }
}

#' Calculate trimming and filtering statistics
#'
#' @param input_file Path to original FASTA/FASTQ file.
#' @param output_file Path to filtered FASTA/FASTQ file.
#' @param format Either "fasta" or "fastq".
#' @param reverse_file Optional reverse input file.
#' @param reverse_output_file Optional reverse output file.
#' @param fastx Source name of input file/object (used for reporting).
#' @param reverse Source name of reverse file/object (used for reporting).
#'
#' @return A tibble with calculated statistics, including:
#' \itemize{
#'   \item \code{Kept_Sequences}: Number of retained sequences.
#'   \item \code{Discarded_Sequences}: Number of discarded sequences.
#'   \item \code{fastx_source}: Name of the file/object with forward (R1) reads.
#'   \item \code{reverse_source}: (If \code{reverse} is specified) Name of the
#'   file/object with reverse (R2) reads.
#' }
#'
#' @noRd
calculate_trim_filt_statistics <- function(input_file,
                                           output_file,
                                           format,
                                           fastx,
                                           reverse = NULL,
                                           reverse_file = NULL,
                                           reverse_output_file = NULL) {

  # Read input and output files
  read_fun <- if (format == "fastq") microseq::readFastq else microseq::readFasta

  input_seqs <- read_fun(input_file)
  output_seqs <- read_fun(output_file)

  n_input <- nrow(input_seqs)
  if (n_input == 0) {
    warning("No sequences in input file: ", input_file)
  }
  n_output <- nrow(output_seqs)

  # Prepare result
  result_table <- tibble::tibble(
    Kept_Sequences = n_output,
    Discarded_Sequences = n_input - n_output,
    fastx_source = fastx
  )

  if (!is.null(reverse) && nzchar(reverse)) {
    result_table$reverse_source <- reverse
  }

  return(result_table)
}

