#' Search for exact full-length matches
#'
#' @description \code{vs_search_exact} searches for exact full-length matches to
#' the query sequences in the database of target sequences using \code{VSEARCH}.
#' Only 100% exact matches are reported, making this command much faster than
#' \code{\link{vs_usearch_global}}.
#'
#' @param fastx_input A FASTA/FASTQ file path or FASTA/FASTQ tibble object
#' containing the query sequences. See \emph{Details}.
#' @param db A FASTA/FASTQ file path or FASTA/FASTQ tibble object containing the
#' target sequences.
#' @param userout A character string specifying the name of the output file for
#' the alignment results. If \code{NULL} (default), no output is written to a
#' file and the results are returned as a tibble with the columns specified in
#' \code{userfields}. See \emph{Details}.
#' @param userfields Fields to include in the output file. Defaults to
#' \code{"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"}. See
#' \emph{Details}.
#' @param strand Specifies which strand to consider when comparing sequences.
#' Can be either \code{"plus"} (default) or \code{"both"}.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Identifies exact full-length matches between query and target sequences
#' using \code{VSEARCH}. Only 100% identical matches are reported, ensuring high
#' specificity.
#'
#' \code{fastx_input} and \code{db} can either be file paths to a FASTA/FASTQ
#' files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
#' columns \code{Header} and \code{Sequence}.FASTQ objects are tibbles that
#' contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' \code{userfields} specifies the fields to include in the output file. Fields
#' must be given as a character string separated by \code{"+"}. The default
#' value of \code{userfields} equals
#' \code{"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"}, which
#' gives a blast-like tab-separated format of twelve fields. See the
#' 'Userfields' section in the \code{VSEARCH} manual for more information.
#'
#' If \code{userout} is specified the alignment results are written to the
#' specified file, and no tibble is returned.
#'
#' If \code{userout} is \code{NULL} a tibble containing the alignment results
#' with the fields specified by \code{userfields} is returned.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{userout} is specified the alignment results are written to the
#' specified file, and no tibble is returned.
#'
#' If \code{userout} is \code{NULL} a tibble containing the alignment results
#' with the fields specified by \code{userfields} is returned.
#'
#' @seealso \code{\link{vs_usearch_global}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#' db <- microseq::readFastq(fastx_input)[1:80, ]
#' userout <- "userout.txt"
#'
#' # Search for exact full-length matches with default parameters, with file as output
#' vs_search_exact(fastx_input = fastx_input,
#'                 db = db,
#'                 userout = userout)
#'
#' # Read results, and give column names
#' outfile_search <- read.delim(userout,
#'                              sep = "\t",
#'                              header = FALSE)
#'
#' colnames(outfile_alignment) <- c("query", "target", "id", "alnlen",
#'                                  "mism", "opens", "qlo", "qhi",
#'                                  "tlo", "thi", "evalue", "bits")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_search_exact search_exact
#'
#' @export
#'
vs_search_exact <- function(fastx_input,
                            db,
                            userout = NULL,
                            userfields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
                            strand = "plus",
                            threads = 1,
                            vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
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

  # Handle input query sequences
  if (!is.character(fastx_input)){
    if ("Quality" %in% colnames(fastx_input)){

      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_file_fastx <- tempfile(pattern = "fastx_input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_fastx)
      microseq::writeFastq(fastx_input, temp_file_fastx)

      fastx_file <- temp_file_fastx

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_fastx <- tempfile(pattern = "fastx_input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file_fastx)
      microseq::writeFasta(fastx_input, temp_file_fastx)

      fastx_file <- temp_file_fastx

    }
  } else {
    if (!file.exists(fastx_input)) stop("Cannot find input file: ", fastx_input)

    fastx_file <- fastx_input
  }

  # Handle input target sequences
  if (!is.character(db)){
    if ("Quality" %in% colnames(db)){

      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(db))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFastq(db, temp_file_db)

      db_file <- temp_file_db

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(db))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFasta(db, temp_file_db)

      db_file <- temp_file_db

    }
  } else {
    if (!file.exists(db)) stop("Cannot find input file: ", db)

    db_file <- db
  }

  # Handle output
  if (is.null(userout)) {
    userout_file <- tempfile(pattern = "userout", fileext = ".txt")
    temp_files <- c(temp_files, userout)
  } else {
    userout_file <- userout
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  db_file <- normalizePath(db_file)

  # Build argument string for command line
  args <- c("--search_exact", shQuote(fastx_file),
            "--db", shQuote(db_file),
            "--userout", userout_file,
            "--userfields", userfields,
            "--threads", threads,
            "--strand", strand)

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(userout)) {

    # Read userout file
    userout_df <- utils::read.delim(userout_file,
                                    sep = "\t",
                                    header = FALSE)

    # Set column names
    columns <- unlist(strsplit(userfields, "\\+"))
    colnames(userout_df) <- columns
  }

  # Return results
  if (is.null(userout)) { # Return tibble
    return(userout_df)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
