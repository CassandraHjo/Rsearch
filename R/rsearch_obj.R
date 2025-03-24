#' Create an Rsearch object
#'
#' @description \code{rsearch_obj} prepares read count, taxonomy, and metadata
#' tables into a format suitable for generating a \code{phyloseq} object.
#'
#' @param readcount_data A file path or a data frame containing OTU count data.
#' The first column must contain OTU identifiers, and subsequent columns must
#' contain read count for each sample. This corresponds to output from the
#' \code{otutabout} argument in \code{cluster_size}.
#' @param tax_data A file path or a data frame containing the taxonomy
#' annotations. The first column must contain OTU identifiers. The second column
#' the OTU names and the second column must contain taxonomic lineage in sintax
#' format. This corresponds to output from the \code{tabbedout} argument in
#' \code{vs_sintax}.
#' @param metadata A file path or a data frame containing metadata for each
#' sample. One column must contain sample identifiers that match the column names
#' in the read count data.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{metadata} that contains sample identifiers. Defaults to
#' \code{"SampleID"}.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item \code{readcount.mat}: A matrix of OTU abundances with OTUs as rows
#'   and samples as columns.
#'   \item \code{tax.mat}: A matrix of taxonomic classifications with OTUs as
#'   rows and ranks (domain to genus) as columns.
#'   \item \code{metadata.df}: A data frame containing metadata for the samples.
#' }
#'
#' OTUs and samples not found in both the abundance and taxonomy/metadata tables
#' are discarded.
#'
#' @examples
#' \dontrun{
#' # Define inputs
#' readcount_data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "SOME_DATA.tsv")
#' tax_data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "SOME_DATA.tsv")
#' metadata <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "SOME_DATA.tsv")
#'
#' # Create Rsearch object
#' obj <- rsearch_obj(readcount_data = readcount_data,
#'                    tax_data = tax_data,
#'                    metadata = metadata,
#'                    sample_id_col = "SampleID")
#'
#' }
#'
#' @export
#'
rsearch_obj <- function(readcount_data,
                        tax_data,
                        metadata,
                        sample_id_col = "SampleID"){

  # Prepare read count data

  if (is.character(readcount_data)) {

    # Read from file
    readcount_data <- suppressMessages(readr::read_delim(readcount_data,
                                                         delim = "\t"))
  }

  otu.names <- dplyr::pull(readcount_data, 1) # Extract OTU names
  readcount.mat <- as.matrix(readcount_data[, -1]) # Extract abundance data
  rownames(readcount.mat) <- otu.names # Set OTU names as rownames

  # Prepare taxonomy data

  if (is.character(tax_data)) {

    # Read from file
    tax_data <- suppressMessages(readr::read_delim(tax_data, delim = "\t",
                                                   col_names = c("Header",
                                                                 "taxonomy",
                                                                 "strand")))
  } else {

    # Assume data frame
    colnames(tax_data) <- c("Header", "taxonomy", "strand")
  }

  tax_data <- tax_data |>
    dplyr::select(1:2) |>
    tidyr::extract(
      col = taxonomy,
      into = c("domain", "phylum", "class", "order", "family", "genus"),
      regex = "d:([^,]+)?(?:,p:([^,]+))?(?:,c:([^,]+))?(?:,o:([^,]+))?(?:,f:([^,]+))?(?:,g:([^,]+))?",
      remove = FALSE
    ) |>
    dplyr::mutate(across(-Header, ~str_remove(., '^[dpcofg]:'))) |>
    dplyr::mutate(across(-Header, ~str_remove(., '\\(.+\\)'))) |>
    dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

  tax.mat <- tax_data |>
    select(-c(Header, taxonomy)) |>
    as.matrix()

  rownames(tax.mat) <- dplyr::pull(tax_data, Header) # Set OTU names as rownames

  # Prepare metadata

  if (is.character(metadata)) {

    # Read from file
    metadata <- suppressMessages(readr::read_delim(metadata, delim = "\t"))
  }

  metadata.df <- metadata |>
    select(-all_of(sample_id_col)) |> # Remove sample ID column
    data.frame() # Convert to data frame

  rownames(metadata.df) <- dplyr::pull(metadata, sample_id_col) # Set sample IDs as rownames

  # Match samples between abundance and metadata

  common_samples <- intersect(colnames(readcount.mat), rownames(metadata.df))
  readcount.mat <- readcount.mat[, common_samples, drop = FALSE]
  metadata.df <- metadata.df[common_samples, , drop = FALSE]

  # Match OTUs between abundance and taxonomy

  common_otus <- intersect(rownames(readcount.mat), rownames(tax.mat))
  readcount.mat <- readcount.mat[common_otus, , drop = FALSE]
  tax.mat <- tax.mat[common_otus, , drop = FALSE]

  # Return list
  return(list(readcount.mat = readcount.mat,
              tax.mat = tax.mat,
              metadata.df = metadata.df))
}


