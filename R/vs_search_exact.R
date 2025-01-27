#' Search for exact full-length matches
#'
#' @description Searches for exact full-length matches to the query sequences in
#' the database of target sequences.
#'
#' @param fastx_input A FASTA/FASTQ file path or object containing the query
#' sequences. See details.
#' @param db A FASTA/FASTQ file path or object containing the target sequences
#' in FASTQ/FASTA format.
#' @param blast6out Name of the output file for the search results in a blast-like
#' tab-separated format of twelve fields, with one line per query-target matching.
#' @param strand \code{"plus"} (default) or \code{"both"}.
#' When comparing sequences only check the plus strand or both strands.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See Details.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#'
#' @details Searches for exact full-length matches to the query sequences in
#' the database of target sequences, using \code{VSEARCH}. Only 100% exact
#' matches are reported. This command is much faster than
#' \code{\link{vs_usearch_global}}.
#'
#' \code{fastx_input} can either be FASTA/FASTQ files or objects. FASTA objects
#' are tibbles that contain the columns \code{Header} and \code{Sequence}.
#' FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#'
#' \code{vsearch_options} can be used to pass additional arguments to \code{VSEARCH},
#' that are not implemented in \code{Rsearch}. See the \code{VSEARCH} manual for
#' additional arguments, and how to use them.
#'
#' @return
#' \code{NULL} (Output is written to file specified by \code{blast6out}).
#'
#' @seealso \code{\link{vs_usearch_global}}
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
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
                            threads = 1,
                            strand = "plus",
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
