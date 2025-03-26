#' Create an Rsearch object
#'
#' @description \code{rsearch_obj} creates a list containing the basic
#' data structures from the processing of a set of samples.
#'
#' @param readcount_data A file path or a data frame containing OTU count data,
#' typically the output from \code{\link{vs_cluster_size}} or similar.
#' The first column must contain OTU identifiers, and subsequent columns must
#' contain read count for each sample.
#' @param sequence_data A file path or a data frame containing the taxonomy
#' for each OTU, typically the output from \code{\link{vs_sintax}} or similar.
#' The first column must contain OTU identifiers and the column must contain
#' the taxonomic lineage in sintax format.
#' @param sample_data A file path or a data frame containing data about each
#' sample. One column must contain sample identifiers that match the column names
#' in the read count data.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{sample_data} that contains sample identifiers. Defaults to
#' \code{"sample_id"}.
#'
#' @details The processing of a set of samples results in some basic data structures.
#' This function collects these into a single list object. We use a basic \code{list}`
#' in the R language instead of constructing a special class, making it fully open
#' and easy to work with.
#'
#' First, for any collection of samples we typically have a table of data about
#' these samples. The input \code{sample_data} provides this, either as an already
#' existing data.frame/tibble in R or as the name of a tab-separated text file.
#' It is assumed the samples are in the rows, and one of the columns *must* contain
#' a text as a unique identifier for each sample. This column name is specified in
#' \code{sample_id_col}.
#'
#' From clustering (\code{\link{vs_cluster_size}} or denoising (\code{\link{vs_unoise}})
#' of the reads we get a table of centroid sequences representing each OTU. The
#' input \code{sequence_data} contains these, where the *first* column contains
#' the OTU identifiers. The \code{Sequence} column must be among the remaining columns. In addition,
#' we typically also do a taxonomic classification of each sequence (see \code{\link{vs_sintax}})
#' and columns describing the taxonomy for each sequence may also typically be in here.
#'
#' Finally, by assigning all reads in each sample to the OTUs (\code{\link{vs_usearch_global}})
#' we get a read count table (input \code{readcount_data}) and this must be a
#' data.frame/tibble with one row for each sample and one column for each OTU.
#' The first column must be the OTU identifiers corresponding to those in the first
#' column of \code{sequence_data}. The remaining columns must have names corresponding to the
#' sample identifiers in the \code{sample_data}. OTUs and samples not found across
#' all data structures are discarded.
#'
#' To convert this tables into a \code{\link{phyloseq}} object, use
#' \code{\link{rsearch2phyloseq}}.
#'
#' @return A straightforward named list with three elements:
#' \itemize{
#'   \item \code{readcount.mat}: A numeric matrix of OTU abundances with OTUs as rows
#'   and samples as columns.
#'   \item \code{sequence.df}: A data.frame with one row for each OTU sequence and
#'   \item \code{sampledata.df}: A data frame containing data about the samples.
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Define inputs
#' readcount.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                             "SOME_DATA.tsv")
#' tax.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                       "SOME_DATA.tsv")
#' sample.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "SOME_DATA.tsv")
#'
#' # Create Rsearch object
#' obj <- rsearch_obj(readcount_data = readcount.data,
#'                    tax_data = tax.data,
#'                    sample_data = sample.data,
#'                    sample_id_col = "SampleID")
#'
#' }
#'
#' @export
#'
rsearch_obj <- function(readcount_data,
                        sequence_data,
                        sample_data,
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
              sampledata.df = metadata.df))
}


