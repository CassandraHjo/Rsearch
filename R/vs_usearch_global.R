#' Global pairwise alignment
#'
#' @description
#' Compare target sequences to the query sequences in FASTA or FASTQ format,
#' using global pairwise alignment.
#'
#' @param fastx_input A FASTA/FASTQ file path or object containing the query
#' sequences. See details.
#' @param db A FASTA/FASTQ file path or object containing the target sequences
#' in FASTQ/FASTA format.
#' @param blast6out Name of the output file for the search results in a blast-like
#' tab-separated format of twelve fields, with one line per query-target matching.
#' @param id The pairwise identity threshold. Defaults to \code{0.8}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}.
#' Defaults to \code{1}.
#'
#' @details
#'
#' \code{fastx_input} can either be FASTA/FASTQ files or objects. FASTA objects
#' are tibbles that contain the columns \code{Header} and \code{Sequence}.
#' FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#' The pairwise identity \code{id} is defined as
#' the number of (matching columns) / (alignment length - terminal gaps)
#'
#'
#' @returns
#' \code{NULL} (Output is written to file specified by \code{blast6out}).
#'
#' @export
#'
vs_usearch_global <- function(fastx_input,
                              db,
                              blast6out,
                              # userout,
                              id = 0.8,
                              threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

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
      print(paste("db_file:", db_file))

    }
  } else {
    if (!file.exists(db)) stop("Cannot find input file: ", db)

    db_file <- db
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file) |>
    shQuote()

  db_file <- normalizePath(db_file) |>
    shQuote()

  # Build argument string for command line
  args <- c("--usearch_global", fastx_file,
            "--db", db_file,
            "--id", id,
            "--blast6out", blast6out,
            # "--userout", userout,
            # "--userfields", "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
            "--threads", threads)

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Return results
  return(invisible(NULL))

}
