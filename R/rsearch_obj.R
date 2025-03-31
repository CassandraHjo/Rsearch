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
#' The first column must contain OTU identifiers and the columns must contain
#' the taxonomic lineage in sintax format.
#' @param sample_data A file path or a data frame containing data about each
#' sample. One column must contain sample identifiers that match the column
#' names in the read count data.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{sample_data} that contains sample identifiers. Defaults to
#' \code{"sample_id"}.
#'
#' @details The processing of a set of samples results in some basic data
#' structures. This function collects these into a single list object. We use a
#' basic \code{list} in the R language instead of constructing a special class,
#' making it fully open and easy to work with.
#'
#' First, for any collection of samples we typically have a table of data about
#' these samples. The input \code{sample_data} provides this, either as an
#' already existing data.frame/tibble in R or as the name of a tab-separated
#' text file. It is assumed the samples are in the rows, and one of the columns
#' \strong{must} contain a text as a unique identifier for each sample. This
#' column name is specified in \code{sample_id_col}.
#'
#' From clustering (\code{\link{vs_cluster_size}}) or denoising
#' (\code{\link{vs_cluster_unoise}}) of the reads we get a table of centroid
#' sequences representing each OTU. The input \code{sequence_data} contains
#' these, where the \emph{first} column contains the OTU identifiers. The
#' \code{Sequence} column must be among the remaining columns. In addition, we
#' typically also do a taxonomic classification of each sequence
#' (see \code{\link{vs_sintax}}) and columns describing the taxonomy for each
#' sequence may also typically be in here.
#'
#' Finally, by assigning all reads in each sample to the OTUs
#' (\code{\link{vs_usearch_global}}) we get a read count table
#' (input \code{readcount_data}) and this must be a data.frame/tibble with one
#' row for each sample and one column for each OTU. The first column must be the
#' OTU identifiers corresponding to those in the first column of
#' \code{sequence_data}. The remaining columns must have names corresponding to
#' the sample identifiers in the \code{sample_data}. OTUs and samples not found
#' across all data structures are discarded.
#'
#' To convert this tables into a \code{\link{phyloseq}} object, use
#' \code{\link{rsearch2phyloseq}}.
#'
#' @return A straightforward named list with three elements:
#' \itemize{
#'   \item \code{readcount.mat}: A numeric matrix of OTU abundances with OTUs as
#'   rows and samples as columns.
#'   \item \code{sequence.df}: A data.frame with one row for each OTU sequence
#'   and
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
                        sample_id_col = "sample_id"){

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

  if (is.character(sequence_data)) {

    # Read from file
    sequence_data <- suppressMessages(utils::read.delim(sequence_data,
                                                        header = TRUE,
                                                        sep = "\t"))
  }

  sequence_data <- sequence_data |>
    dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

  # Prepare metadata

  if (is.character(sample_data)) {

    # Read from file
    sample_data <- suppressMessages(utils::read.delim(sample_data,
                                                      sep = "\t",
                                                      header = TRUE))
  }

  # Match samples between read count data and metadata

  common_samples <- intersect(colnames(readcount.mat),
                              dplyr::pull(sample_data, sample_id_col))
  readcount.mat <- readcount.mat[, common_samples, drop = FALSE]
  sampledata.df <- sample_data |>
    dplyr::filter(.data[[sample_id_col]] %in% common_samples)

  # Match OTUs between abundance and taxonomy

  common_otus <- intersect(rownames(readcount.mat), sequence_data$Header)
  readcount.mat <- readcount.mat[common_otus, , drop = FALSE]
  sequence.df <- sequence_data |>
    dplyr::filter(Header %in% common_otus)

  # Return list
  return(list(readcount.mat = readcount.mat,
              sequence.df = sequence.df,
              sampledata.df = sampledata.df))
}

#' Convert Rsearch to phyloseq object
#'
#' @description Creating a phyloseq object from a Rsearch object.
#'
#' @param rsearch.obj A Rsearch object, see \code{\link{rsearch_obj}}.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{sampledata.df} that contains sample identifiers. Defaults to
#' \code{"sample_id"}.
#'
#' @details This function converts an Rsearch object, which is a simple
#' \code{list}, to a \code{\link{phyloseq}} object from the \code{phyloseq} R
#' package.
#'
#' @return A \code{\link{phyloseq}} object.
#'
#' @export
#'
rsearch2phyloseq <- function(rsearch.obj, sample_id_col = "sample_id"){
  otu.table <- rsearch.obj$readcount.mat

  sample.data <- rsearch.obj$sampledata.df
  rownames(sample.data) <- rsearch.obj$sampledata.df[,sample_id_col]

  taxonomy.tbl <- dplyr::select(rsearch.obj$sequence.df, -c(Header, Sequence))

  if(ncol(taxonomy.tbl) > 0){
    tax.mat <- as.matrix(taxonomy.tbl)
    rownames(tax.mat) <- rsearch.obj$sequence.df$Header
    ps.obj <- phyloseq::phyloseq(phyloseq::otu_table(otu.table,
                                                     taxa_are_rows = T),
                                 phyloseq::sample_data(sample.data),
                                 phyloseq::tax_table(tax.mat))
  } else {
    ps.obj <- phyloseq::phyloseq(phyloseq::otu_table(otu.table,
                                                     taxa_are_rows = T),
                                 phyloseq::sample_data(sample.data))
  }
  return(ps.obj)
}

#' @title Convert phyloseq to Rsearch object
#'
#' @description Creating a simple list from a phyloseq object.
#'
#' @param phyloseq.obj A phyloseq object, see \code{\link{phyloseq}}.
#'
#' @details This function converts a phyloseq object to a simple
#' \code{\link{list}} giving the entries the names as in
#' \code{\link{rsearch_obj}}.
#'
#' This may be convenient for some wrangling on the data, and then perhaps
#' converting it back to a phyloseq object again with
#' \code{\link{rsearch2phyloseq}}.
#'
#' @return A \code{list} with entries as in a Rsearch object, except that the
#' \code{sequence.tbl} do not contain sequences, only taxonomy.
#'
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @export
#'
phyloseq2rsearch <- function(phyloseq.obj){
  lst <- list(
    sampledata.df = as.data.frame(as.matrix(phyloseq::sample_data(phyloseq.obj))),
    readcount.mat = as.matrix(as.data.frame(phyloseq::otu_table(phyloseq.obj))),
    sequence.df = as.data.frame(phyloseq::tax_table(phyloseq.obj)))
  return(lst)
}
