#' Create an Rsearch object
#'
#' @description \code{rsearch_obj} Standardizes and creates a list containing
#' three elements with data structures that can be used as input to build a
#' phyloseq object in the phyloseq package.
#'
#' @param readcount_data A file path or a data frame (or tibble) containing OTU
#' count data, typically the output from \code{\link{vs_cluster_size}} or
#' similar. This must have one row per OTU and one column per sample. The first
#' column must contain OTU identifiers corresponding to those in the first
#' column of \code{sequence_data}, and the remaining columns must have names
#' matching the sample identifiers in \code{sample_data}. OTUs and samples not
#' found across all data structures are discarded.
#' @param sequence_data A file path or a data frame (or tibble) containing
#' centroid sequences representing each OTU, typically obtained from
#' clustering (\code{\link{vs_cluster_size}}) or denoising
#' (\code{\link{vs_cluster_unoise}}). The first column must contain OTU
#' identifiers. One of the remaining columns must be named \code{Sequence},
#' containing the actual DNA sequences. Additional columns may include taxonomic
#' classification data, e.g. from \code{\link{vs_sintax}}.
#' @param sample_data A file path or a data frame (or tibble) containing
#' metadata about each sample. Samples are assumed to be in rows, and one of
#' the columns \strong{must} contain a unique identifier for each sample that
#' matches the column names in \code{readcount_data}.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{sample_data} that contains the unique sample identifiers.
#' This column will be used to match sample metadata to read count data.
#' Defaults to \code{"sample_id"}.
#'
#' @details This function accepts three datasets—read count data, sequence data,
#' and sample metadata—standardizes them, and returns a streamlined input suitable
#' for constructing a phyloseq object using the phyloseq package. The function
#' combines these into a single list object. The implementation uses a standard
#' \code{list} in R rather than a specialized class providing an open and easily
#' accessible structure.
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
#' @examples
#' \dontrun{
#' # Define inputs
#' readcount.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                             "SOME_DATA.tsv")
#' sequence.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                       "SOME_DATA.tsv")
#' sample.data <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "SOME_DATA.tsv")
#'
#' # Create Rsearch object
#' obj <- rsearch_obj(readcount_data = readcount.data,
#'                    sequence_data = sequence.data,
#'                    sample_data = sample.data,
#'                    sample_id_col = "SampleID")
#'
#' # Convert Rsearch object to phyloseq object
#' phy_obj <- rsearch2phyloseq(obj, sample_id_col = "SampleID")
#'
#' # Convert phyloseq object to Rsearch object
#' rsearch_obj <- phyloseq2rsearch(phy_obj)
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
#' @examples
#' # HER_TRENGS_KODE
#'
#' @seealso
#' \code{\link{rsearch_obj}}
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

#' Convert phyloseq to Rsearch object
#'
#' @description Creating a simple list from a phyloseq object.
#'
#' @param phyloseq.obj A phyloseq object, see \code{\link{phyloseq}}.
#'
#' @details This function converts a phyloseq object to a simple
#' \code{\link{list}} with three elements as dataframes (or tibbles). The
#' entries are named according to the structure used in
#' \code{\link{rsearch_obj}}
#'
#' @return A \code{list} with entries as in a Rsearch object, except that the
#' \code{sequence.tbl} do not contain sequences, only taxonomy.
#'
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @examples
#' # HER_TRENGS_KODE
#'
#'
#' @seealso
#' \code{\link{rsearch_obj}}
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
