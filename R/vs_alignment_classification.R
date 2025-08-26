#' Taxonomic classification with LCA
#'
#' @description \code{vs_alignment_classification} assigns taxonomy by global
#' alignment and Last Common Ancestor (LCA) consensus of database hits using
#' \code{VSEARCH}.
#'
#' @param fastx_input (Required). A FASTA/FASTQ file path or FASTA/FASTQ object.
#' See \emph{Details}.
#' @param database (Required). A FASTA/FASTQ file path or FASTA/FASTQ tibble
#' object containing the target sequences.
#' @param lcaout (Optional). A character string specifying the name of the
#' output file. If \code{NULL} (default), no output is
#' written to a file and the results are returned as a tibble with the columns
#' \code{query_id} and \code{taxonomy}.
#' @param lca_cutoff (Optional). Adjust the fraction of matching hits required
#' for the last common ancestor (LCA). Defaults to \code{1.0}, which requires
#' all hits to match at each taxonomic rank for that rank to be included. If a
#' lower cutoff value is used, e.g. 0.95, a small fraction of non-matching hits
#' are allowed while that rank will still be reported. The argument to this
#' option must be between \code{0.5} and \code{1.0}.
#' @param top_hits_only (Optional). If \code{TRUE}, only the top hits with an
#' equally high percentage of identity between the query and database sequence
#' sets are written to the output. Defaults to \code{FALSE}.
#' @param gapopen (Optional). Penalties for gap opening. Defaults to
#' \code{"20I/2E"}. See \emph{Details}.
#' @param gapext (Optional). Penalties for gap extension. Defaults to
#' \code{"2I/1E"}. See \emph{Details}.
#' @param id (Optional). Pairwise identity threshold. Defines the minimum
#' identity required for matches. Defaults to \code{0.7}.
#' @param strand (Optional). Specifies which strand to consider when comparing
#' sequences. Can be either \code{"plus"} (default) or \code{"both"}.
#' @param maxaccepts (Optional). Maximum number of matching target sequences to
#' accept before stopping the search for a given query. Defaults to \code{2}.
#' Must be larger than \code{1} for information to be useful.
#' @param maxrejects (Optional). Maximum number of non-matching target sequences
#' to consider before stopping the search for a given query. Defaults to 32. If
#' \code{maxaccepts} and \code{maxrejects} are both set to 0, the complete
#' database is searched.
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
#' Performs global sequence alignment against a reference database and assigns
#' taxonomy using the Last Common Ancestor (LCA) approach, reporting the deepest
#' taxonomic level consistently supported by the majority of hits.
#'
#' \code{fastx_input} and \code{database} can either be file paths to a
#' FASTA/FASTQ files or FASTA/FASTQ objects. FASTA objects are tibbles that
#' contain the columns \code{Header} and \code{Sequence}, see
#' \code{\link[microseq]{readFasta}}. FASTQ objects are tibbles that contain the
#' columns \code{Header}, \code{Sequence}, and \code{Quality}, see
#' \code{\link[microseq]{readFastq}}.
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
#' If \code{lcaout} is specified the results are written to the specified file.
#' If \code{lcaout} is \code{NULL} a data.frame is returned.
#'
#' The data.frame contains the classification results for each query sequence.
#' Both the \code{Header} and \code{Sequence} columns of \code{fasta_input} are
#' copied into this table, and in addition are also the columns for each rank.
#' The ranks depend on the database file used, but are typically domain, phylum,
#' class, order,family, genus and species.
#'
#' @examples
#' \dontrun{
#' # Example files
#' db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "sintax_db.fasta")
#' fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small.fasta")
#'
#' tax.tbl <- vs_alignment_classification(fastx_input = fasta.file,
#'                                        database = db.file)
#' View(tax.tbl)
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_alignment_classification alignment_classification lca
#' lca_classification
#'
#' @export
#'
vs_alignment_classification <- function(fastx_input,
                                        database,
                                        lcaout = NULL,
                                        lca_cutoff = 1.0,
                                        top_hits_only = FALSE,
                                        gapopen = "20I/2E",
                                        gapext = "2I/1E",
                                        id = 0.7,
                                        strand = "plus",
                                        maxaccepts = 2,
                                        maxrejects = 32,
                                        threads = 1,
                                        vsearch_options = NULL,
                                        tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  # Validate lca_cutoff
  if (lca_cutoff > 1.0 || lca_cutoff <= 0.5) {
    stop("Invalid value for 'lca_cutoff'. Must be between 0.5 and 1.0.")
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

      temp_file_fastx <- tempfile(pattern = "fastx_input",
                                  tmpdir = tmpdir,
                                  fileext = ".fq")
      temp_files <- c(temp_files, temp_file_fastx)
      microseq::writeFastq(fastx_input, temp_file_fastx)

      fastx_file <- temp_file_fastx

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_fastx <- tempfile(pattern = "fastx_input",
                                  tmpdir = tmpdir,
                                  fileext = ".fa")
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

      temp_file_db <- tempfile(pattern = "db_input",
                               tmpdir = tmpdir,
                               fileext = ".fq")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFastq(database, temp_file_db)

      db_file <- temp_file_db

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(database))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file_db <- tempfile(pattern = "db_input",
                               tmpdir = tmpdir,
                               fileext = ".fa")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFasta(database, temp_file_db)

      db_file <- temp_file_db

    }
  } else {
    if (!file.exists(database)) stop("Cannot find input file: ", database)

    db_file <- database
  }

  # Determine output file based on user input
  if (!is.null(lcaout)) {
    outfile <- lcaout
  } else {
    outfile <- tempfile(pattern = "lcaout",
                        tmpdir = tmpdir,
                        fileext = ".txt")

    temp_files <- c(temp_files, outfile)
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  db_file <- normalizePath(db_file)

  # Build argument string for command line
  args <- c("--usearch_global", shQuote(fastx_file),
            "--db", shQuote(db_file),
            "--lcaout", shQuote(outfile),
            "--lca_cutoff", lca_cutoff,
            "--id", id,
            "--threads", threads,
            "--strand", strand,
            "--gapopen", gapopen,
            "--gapext", gapext,
            "--maxaccepts", maxaccepts,
            "--maxrejects", maxrejects)

  # Add top_hits_only if specified
  if (top_hits_only) {
    args <- c(args, "--top_hits_only", "")
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

  # Determine return output
  if (!is.null(lcaout)) {
    return(invisible(NULL)) # No return if lcaout is specified
  } else {
    # lcaout_df <- suppressMessages(
    #   readr::read_delim(outfile,
    #                     delim = "\t",
    #                     col_names = c("query_id", "taxonomy")))

    # The output table
    if (is.character(fastx_input)){
      fastx_input <- microseq::readFastq(fastx_input) |>
        dplyr::mutate(Header = stringr::str_remove(Header, "^>"))
    }

    lcaout_df <- fastx_input |>
      dplyr::select(Header, Sequence)
    lcaout_df <- suppressMessages(
      readr::read_delim(outfile,
                        delim = "\t",
                        col_names = c("Header", "taxonomy"))) |>
      dplyr::mutate(domain = stringr::str_extract(taxonomy, "(?<=d:)[^,]+")) |>
      dplyr::mutate(phylum = stringr::str_extract(taxonomy, "(?<=p:)[^,]+")) |>
      dplyr::mutate(class = stringr::str_extract(taxonomy, "(?<=c:)[^,]+")) |>
      dplyr::mutate(order = stringr::str_extract(taxonomy, "(?<=o:)[^,]+")) |>
      dplyr::mutate(family = stringr::str_extract(taxonomy, "(?<=f:)[^,]+")) |>
      dplyr::mutate(genus = stringr::str_extract(taxonomy, "(?<=g:)[^,]+")) |>
      dplyr::mutate(species = stringr::str_extract(taxonomy, "(?<=s:)[^,]+")) |>
      dplyr::select(-taxonomy) |>
      dplyr::right_join(lcaout_df, by = "Header") |>
      dplyr::relocate(Sequence, .after = tidyr::last_col())
    return(lcaout_df)
  }
}
