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
#' @param blast6out Name of the output file for the search results in a
#' blast-like tab-separated format of twelve fields, with one line per
#' query-target matching.
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
#' Results are written to the file specified by \code{blast6out} in a blast-like
#' tab-separated format containing twelve fields:
#' \itemize{
#'   \item \code{Query ID}
#'   \item \code{Target ID}
#'   \item \code{Percent Identity}
#'   \item \code{Alignment Length}
#'   \item \code{Number of Mismatches}
#'   \item \code{Number of Gap Openings}
#'   \item \code{Query Start}
#'   \item \code{Query End}
#'   \item \code{Target Start}
#'   \item \code{Target End}
#'   \item \code{E-value}
#'   \item \code{Bit Score}
#' }
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return
#' \code{NULL}. Results are written directly to the file specified by
#' \code{blast6out}).
#'
#' @seealso \code{\link{vs_usearch_global}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "R1_sample1_small.fq")
#' db <- microseq::readFastq(fastx_input)[1:80, ]
#' blast6out <- "blast6out.txt"
#'
#' # Search for exact full-length matches with default parameters, with file as output
#' vs_search_exact(fastx_input = fastx_input,
#'                 db = db,
#'                 blast6out = blast6out)
#'
#' # Read results, and give column names
#' outfile_search <- read.delim(blast6out,
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
                            blast6out,
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

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  db_file <- normalizePath(db_file)

  # Build argument string for command line
  args <- c("--search_exact", shQuote(fastx_file),
            "--db", shQuote(db_file),
            "--blast6out", blast6out,
            # "--userout", userout,
            # "--userfields", "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
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

  # Return results
  return(invisible(NULL))

}
