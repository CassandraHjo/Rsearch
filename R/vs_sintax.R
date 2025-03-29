#' Taxonomic classification
#'
#' @description \code{vs_sintax} classifies sequences using the Sintax algorithm
#' implemented in \code{VSEARCH}.
#'
#' @param fasta_input A FASTA file path or a FASTA object with reads to
#' classify, see \emph{Details}.
#' @param database A FASTA file path or a FASTA object containing the
#' reference database in FASTA format. The sequences need to be annotated
#' with taxonomy, see \emph{Details}.
#' @param outfile Name of the output file. If \code{NULL}
#' (default), results are returned as a data.frame.
#' @param cutoff minimum level of bootstrap support (0.0-1.0) for the classifications.
#' Defaults to \code{0.0}.
#' @param strand Specifies which strand to consider when comparing sequences.
#' Can be either \code{"plus"} (default) or \code{"both"}.
#' @param randseed Seed for the random number generator used in the Sintax
#' algorithm. Defaults to \code{NULL}.
#' @param logfile Name of the log file to capture messages from \code{VSEARCH}.
#' If \code{NULL} (default), no log file is created.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details The sequences in the input file are classified according to the
#' Sintax algorithm, using \code{VSEARCH}, see https://www.biorxiv.org/content/10.1101/074161v1
#'
#' \code{fasta_input} can either be a file path to a FASTA file or a
#' FASTA object.
#'
#' \code{database} can either be a file path to a FASTA file or a
#' FASTA object. FASTA objects are tibbles that contain the
#' text columns \code{Header} and \code{Sequence}, see \code{\link{readFasta}}.
#' The \code{Header} texts of this file must follow the sintax-pattern, see
#' \code{\link{make_sintax_database}}.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @returns If \code{outfile} is \code{NULL} a data.frame is returned. If it contains
#' a file name (text) the data.frame is written to that file with tab-separated columns.
#'
#' The data.frame contains the classification results for each input sequence.
#' Both the \code{Header} and \code{Sequence} columns of \code{input_fasta} are
#' copied into this table, and in addition are also the columns for each rank. The ranks
#' depend on the database file used, but are typically domain, phylum, class, order,
#' family, genus and species. For each classification is also a bootstrap support score.
#' These are in separate columns with corresponding names, i.e. domain_score,
#' phylum_score, etc.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_sintax sintax classify
#'
#' @export
#'
vs_sintax <- function(fasta_input,
                      database,
                      outfile = NULL,
                      cutoff = NULL,
                      strand = "plus",
                      randseed = NULL,
                      logfile = NULL,
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

  # Check if FASTA input is file or tibble
  if (!is.character(fasta_input)){
    temp_file <- tempfile(pattern = "input", fileext = ".fa")
    temp_files <- c(temp_files, temp_file)
    microseq::writeFasta(fasta_input, temp_file)
    fasta_file <- temp_file

    # Capture original name for statistics table later
    fasta_input_name <- as.character(substitute(fasta_input))
  } else {
    fasta_file <- fasta_input

    # Capture original name for statistics table later
    fasta_input_name <- basename(fasta_input)
    fasta_input <- microseq::readFasta(fasta_file)
  }

  # Check if input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Handle input data base sequences
  if (!is.character(database)){
    # Validate tibble
    required_cols <- c("Header", "Sequence")
    if (!all(required_cols %in% colnames(database))) {
      stop("FASTA data base must contain columns: Header and Sequence")
    }

    temp_file_db <- tempfile(pattern = "db_input", fileext = ".fa")
    temp_files <- c(temp_files, temp_file_db)
    microseq::writeFasta(db, temp_file_db)

    db_file <- temp_file_db

  } else {
    if (!file.exists(database)) stop("Cannot find input file: ", database)

    db_file <- database
  }

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)
  db_file <- normalizePath(db_file)

  # The temporary outfile
  tmp_outfile <- tempfile(pattern = "tmp_output", fileext = ".txt")
  temp_files <- c(temp_files, tmp_outfile)

  # Build argument string for command line
  args <- c("--sintax", shQuote(fasta_file),
            "--db", shQuote(db_file),
            "--threads", threads,
            "--strand", strand,
            "--tabbedout", tmp_outfile)

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add log file if specified
  if (!is.null(logfile)){
    args <- c(args, "--log", logfile)
  }

  # Add cutoff if specified
  if (!is.null(cutoff)){
    args <- c(args, "--sintax_cutoff", cutoff)
  }

  # Add random seed for sintax if specified
  if (!is.null(randseed)){
    args <- c(args, "--randseed", randseed)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # The output table
  out.tbl <- fasta_input
  out.tbl <- utils::read.table(tmp_outfile, sep = "\t", col.names = c("Header", "taxonomy", "plus")) |>
    dplyr::select(-plus) |>
    dplyr::mutate(domain = stringr::str_extract(taxonomy, "d:.+?\\)")) |>
    dplyr::mutate(domain_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(domain, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(domain = stringr::str_remove_all(domain, "d:|\\(.+")) |>
    dplyr::mutate(phylum = stringr::str_extract(taxonomy, "p:.+?\\)")) |>
    dplyr::mutate(phylum_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(phylum, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(phylum = stringr::str_remove_all(phylum, "p:|\\(.+")) |>
    dplyr::mutate(class = stringr::str_extract(taxonomy, "c:.+?\\)")) |>
    dplyr::mutate(class_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(class, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(class = stringr::str_remove_all(class, "c:|\\(.+")) |>
    dplyr::mutate(order = stringr::str_extract(taxonomy, "o:.+?\\)")) |>
    dplyr::mutate(order_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(order, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(order = stringr::str_remove_all(order, "o:|\\(.+")) |>
    dplyr::mutate(family = stringr::str_extract(taxonomy, "f:.+?\\)")) |>
    dplyr::mutate(family_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(family, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(family = stringr::str_remove_all(family, "f:|\\(.+")) |>
    dplyr::mutate(genus = stringr::str_extract(taxonomy, "g:.+?\\)")) |>
    dplyr::mutate(genus_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(genus, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(genus = stringr::str_remove_all(genus, "g:|\\(.+")) |>
    dplyr::mutate(species = stringr::str_extract(taxonomy, "s:.+?\\)")) |>
    dplyr::mutate(species_score = suppressWarnings(as.numeric(stringr::str_remove_all(stringr::str_extract(species, "\\(.+\\)"), "\\(|\\)")))) |>
    dplyr::mutate(species = stringr::str_remove_all(species, "s:|\\(.+")) |>
    dplyr::select(-taxonomy) |>
    dplyr::right_join(out.tbl, by = "Header") |>
    dplyr::relocate(Sequence, .after = tidyr::last_col())

  # Write out.tbl to file
  if (!is.null(outfile)) {

    # Read output into R tibble
    utils::write.table(out.tbl,
                       file = outfile,
                       sep = "\t")
  }

  # Return results
  if (is.null(outfile)) { # Return tibble
    return(out.tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}


#' Make Sintax database
#'
#' @description Creates a properly formatted fasta file for the use as a Sintax database.
#'
#' @param taxonomy_table A data.frame with sequences and proper information for making
#' a Sintax database, see \emph{Details}.
#' @param outfile Name of database file to create (a fasta file).
#'
#' @details The Sintax algorithm is used
#' by \code{VSEARCH} to assign taxonomic information to 16S sequences. It requires a
#' database, which is nothing but a fasta file of 16S sequences with properly formatted
#' \code{Header}-lines.
#'
#' The \code{taxonomy_table} provided as input here must have the columns:
#'
#' \itemize{
#'  \item \code{Header} - short unique text for each sequence
#'  \item \code{Sequence} - the sequences
#'  \item Columns \code{domain}, \code{phylum}, \code{class}, \code{order},
#'   \code{family}, \code{genus}, \code{species}. Text columns with taxon names.
#' }
#'
#' In some taxonomies the domain rank is named kingdom, but here we use the
#' word domain. You may very well have empty (NA) entries in the taxonomy columns
#' of the table.
#'
#' @returns No return in R, but a fasta file (\code{outfile}) with properly
#' formatted \code{Header} lines is created.
#'
#' @references \url{https://www.biorxiv.org/content/10.1101/074161v1}
#'
#' @aliases sintax
#'
#' @export
#'
make_sintax_db <- function(taxonomy_table, outfile){
  if(!exists("Header", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named Header, with a unique text for each sequence")
  }
  if(!exists("Sequence", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named Sequence, with the sequences")
  }
  if(!exists("domain", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named domain")
  }
  if(!exists("phylum", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named phylum")
  }
  if(!exists("class", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named class")
  }
  if(!exists("order", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named order")
  }
  if(!exists("family", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named family")
  }
  if(!exists("genus", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named genus")
  }
  if(!exists("species", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named species")
  }

  sintax.tbl <- taxonomy_table |>
    dplyr::mutate(Header = stringr::str_c(Header, ";tax=",
                                          "d:", domain, ",",
                                          "p:", phylum, ",",
                                          "c:", class, ",",
                                          "o:", order, ",",
                                          "f:", family, ",",
                                          "g:", genus, ",",
                                          "s:", species, ";"))

  writeFasta(sintax.tbl, out.file = outfile)
  return(invisible(NULL))
}
