#' Global pairwise alignment
#'
#' @description
#' \code{vs_usearch_global} performs global pairwise alignment of query
#' sequences against target sequences using \code{VSEARCH}.
#'
#' @param fastx_input A FASTA/FASTQ file path or FASTA/FASTQ object. See
#' \emph{Details}.
#' @param database A FASTA/FASTQ file path or FASTA/FASTQ tibble object containing the
#' target sequences.
#' @param userout A character string specifying the name of the output file for
#' the alignment results. If \code{NULL} (default), no output is written to a
#' file and the results are returned as a tibble with the columns specified in
#' \code{userfields}. See \emph{Details}.
#' @param otutabout A character string specifying the name of the output file in
#' an OTU table format. If \code{NULL} (default), no output is written to a file.
#' If \code{TRUE}, the output is returned as a tibble. See \emph{Details}.
#' @param userfields Fields to include in the output file. Defaults to
#' \code{"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"}. See
#' \emph{Details}.
#' @param gapopen Penalties for gap opening. Defaults to \code{"20I/2E"}. See
#' \emph{Details}.
#' @param gapext Penalties for gap extension. Defaults to \code{"2I/1E"}. See
#' \emph{Details}.
#' @param id Pairwise identity threshold. Defines the minimum identity required
#' for matches. Defaults to \code{0.7}.
#' @param strand Specifies which strand to consider when comparing sequences.
#' Can be either \code{"plus"} (default) or \code{"both"}.
#' @param maxaccepts Maximum number of matching target sequences to accept
#' before stopping the search for a given query. Defaults to \code{1}. Only
#' works when \code{strand} is set to \code{"plus"} (default).
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param vsearch_options Additional arguments to pass to \code{VSEARCH}.
#' Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Performs global pairwise alignment between query and target sequences using
#' \code{VSEARCH}, and reports matches based on the specified pairwise identity
#' threshold (\code{id}). Only alignments that meet or exceed the identity
#' threshold are included in the output.
#'
#' \code{fastx_input} and \code{database} can either be file paths to a FASTA/FASTQ
#' files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
#' columns \code{Header} and \code{Sequence}, see \code{\link{readFasta}}. FASTQ
#' objects are tibbles that contain the columns \code{Header}, \code{Sequence},
#' and \code{Quality}, see \code{\link{readFastq}}.
#'
#' \code{userfields} specifies the fields to include in the output file. Fields
#' must be given as a character string separated by \code{"+"}. The default
#' value of \code{userfields} equals
#' \code{"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"}, which
#' gives a blast-like tab-separated format of twelve fields. See the
#' 'Userfields' section in the \code{VSEARCH} manual for more information.
#'
#' \code{otutabout} gives the option to output the results in an OTU
#' table format with tab-separated columns. The first line will start with
#' the string \"#OTU ID\" and is followed by a tab-separated list of all sample
#' identifiers (\"sample=X\"). The following lines, one for each OTU, start with
#' the OTU identifier and are followed by a tab-separated list of abundances for
#' that OTU in each sample. If \code{otutabout} is a character string, the output
#' is written to the specified file. If \code{otutabout} is \code{TRUE}, the
#' function returns the OTU table as a tibble.
#'
#' Pairwise identity (\code{id}) is calculated as the number of matching columns
#' divided by the alignment length minus terminal gaps.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' Visit the \code{VSEARCH}
#' \href{https://github.com/torognes/vsearch?tab=readme-ov-file#getting-help}{documentation}
#' for information about defining \code{gapopen} and \code{gapext}.
#'
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{userout} is specified the alignment results are written to the
#' specified file, and no tibble is returned. If \code{userout} is \code{NULL} a
#' tibble containing the alignment results with the fields specified by
#' \code{userfields} is returned.
#'
#' If \code{otutabout} is \code{TRUE}, an OTU table is returned as a tibble.
#' If \code{otutabout} is a character string, the output is written to the file,
#' and no tibble is returned.
#'
#' If neither \code{userout} nor \code{otutabout} is specified, a tibble
#' containing the alignment results is returned.
#'
#' @examples
#' \dontrun{
#' # You would typically use something else as database
#' query_file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small.fasta")
#' db <- query_file
#'
#' # Run global pairwise alignment with default parameters and write results to file
#' vs_usearch_global(fastx_input = query_file,
#'                   database = db,
#'                   userout = "delete_me.txt")
#'
#' # Read results, and give column names
#' result.tbl <- read.table("delete_me.txt",
#'                          sep = "\t",
#'                          header = FALSE,
#'                          col.names = c("query", "target", "id", "alnlen",
#'                                        "mism", "opens", "qlo", "qhi",
#'                                        "tlo", "thi", "evalue", "bits"))
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_usearch_global usearch_global global_alignment
#'
#' @export
#'
vs_usearch_global <- function(fastx_input,
                              database,
                              userout = NULL,
                              otutabout = NULL,
                              userfields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
                              gapopen = "20I/2E",
                              gapext = "2I/1E",
                              id = 0.7,
                              strand = "plus",
                              maxaccepts = 1,
                              threads = 1,
                              vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  # Ensure only one output format is specified
  if (!is.null(userout) && !is.null(otutabout)) {
    stop("Only one of 'userout' or 'otutabout' can be specified.")
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
  if (!is.character(database)){
    if ("Quality" %in% colnames(database)){

      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(database))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFastq(database, temp_file_db)

      db_file <- temp_file_db

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(database))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFasta(database, temp_file_db)

      db_file <- temp_file_db

    }
  } else {
    if (!file.exists(database)) stop("Cannot find input file: ", database)

    db_file <- database
  }

  # Determine output file based on user input
  if (!is.null(userout)) {
    outfile <- userout
  } else if (!is.null(otutabout)) {
    outfile <- ifelse(is.character(otutabout), otutabout, tempfile(pattern = "otutable", fileext = ".tsv"))
  } else {
    outfile <- tempfile(pattern = "userout", fileext = ".txt")
  }

  if (is.null(userout) && (is.null(otutabout) || !is.character(otutabout))) {
    temp_files <- c(temp_files, outfile)
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  db_file <- normalizePath(db_file)

  # Build argument string for command line
  args <- c("--usearch_global", shQuote(fastx_file),
            "--db", shQuote(db_file),
            "--id", id,
            "--threads", threads,
            "--strand", strand,
            "--gapopen", gapopen,
            "--gapext", gapext)

  if (!is.null(userout)) {
    args <- c(args, "--userout", outfile, "--userfields", userfields)
  } else if (!is.null(otutabout)) {
    args <- c(args, "--otutabout", outfile)
  } else {
    args <- c(args, "--userout", outfile, "--userfields", userfields) # Default output
  }

  # Add maxaccepts if strand is "plus"
  if (strand == "plus") {
    args <- c(args, "--maxaccepts", maxaccepts)
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Determine return output
  if (!is.null(userout)) {
    return(invisible(NULL)) # No return if userout is specified
  } else if (!is.null(otutabout)) {
    if (is.character(otutabout)) {
      return(invisible(NULL)) # File output only
    } else {
      return(suppressMessages(readr::read_delim(outfile))) # Return as tibble
    }
  } else {
    userout_df <- suppressMessages(readr::read_delim(outfile, col_names = FALSE))
    columns <- unlist(strsplit(userfields, "\\+"))
    colnames(userout_df) <- columns
    return(userout_df)
  }
}
