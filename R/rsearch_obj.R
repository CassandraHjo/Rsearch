#' Create Rsearch object
#'
#' @description \code{rsearch_obj} standardizes and organizes data into an
#' Rsearch object.
#'
#' creates a list containing
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
#' (\code{\link{vs_cluster_unoise}}). The first column must be called
#' \code{Header} and contain OTU identifiers. One of the remaining columns must
#' be named \code{Sequence}, containing the actual DNA sequences. Additional
#' columns may include taxonomic classification data, e.g. from
#' \code{\link{vs_sintax}}.
#' @param sample_data A file path or a data frame (or tibble) containing
#' metadata about each sample. Samples are assumed to be in rows, and one of
#' the columns \strong{must} contain a unique identifier for each sample that
#' matches the column names in \code{readcount_data}.
#' @param sample_id_col A character string specifying the name of the column in
#' \code{sample_data} that contains the unique sample identifiers.
#' This column will be used to match sample metadata to read count data.
#' Defaults to \code{"sample_id"}.
#'
#' @details This function standardizes and organizes data into an
#' Rsearch object: a structured three key data components used or generated
#' during the Rsearch workflow: read count data, sequence data, and sample data.
#'
#' The function accepts three datasets—read count data, sequence data,
#' and sample metadata, and returns a streamlined input
#' suitable for constructing a phyloseq object using the
#' \code{\link{rsearch2phyloseq}} function. The implementation uses a
#' standard \code{list} in R rather than a specialized class providing an open
#' and easily accessible structure.
#'
#' To convert this object into a \code{\link{phyloseq}} object, use
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
#' readcount.dta <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                            "readcount_data.tsv")
#' sequence.dta <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                           "sequence_data.tsv")
#' sample.dta <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                         "sample_data.tsv")
#'
#' # Create Rsearch object
#' obj <- rsearch_obj(readcount_data = readcount.dta,
#'                    sequence_data = sequence.dta,
#'                    sample_data = sample.dta,
#'                    sample_id_col = "sample_id")
#'
#' # Convert Rsearch object to phyloseq object
#' phy_obj <- rsearch2phyloseq(obj, sample_id_col = "sample_id")
#'
#' # Convert phyloseq object to Rsearch object
#' rsearch_obj <- phyloseq2rsearch(phy_obj)
#'
#' }
#'
#' @seealso \link{rsearch2phyloseq} \link{phyloseq2rsearch}
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
    sequence_data <- suppressMessages(readr::read_delim(sequence_data,
                                                        delim = "\t"))
  }

  sequence_data <- sequence_data |>
    dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))

  # Prepare metadata

  if (is.character(sample_data)) {

    # Read from file
    sample_data <- suppressMessages(readr::read_delim(sample_data,
                                                      delim = "\t"))
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

#' Convert Rsearch object to phyloseq object
#'
#' @description \code{rsearch2phyloseq} converts an Rsearch object to a phyloseq
#' object.
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
#' @references
#' \url{https://joey711.github.io/phyloseq/}
#'
#' @examples
#' \dontrun{
#' # Convert Rsearch object to phyloseq object
#' phy_obj <- rsearch2phyloseq(obj, sample_id_col = "sample_id")
#' }
#'
#' @seealso
#' \code{\link{rsearch_obj}}
#'
#' @export
#'
rsearch2phyloseq <- function(rsearch.obj, sample_id_col = "sample_id"){

  otu.table <- rsearch.obj$readcount.mat
  sample.dta <- as.data.frame(rsearch.obj$sampledata.df)
  rownames(sample.dta) <- sample.dta[[sample_id_col]]

  taxonomy.tbl <- dplyr::select(rsearch.obj$sequence.df, -c(Header, Sequence))

  if(ncol(taxonomy.tbl) > 0){
    tax.mat <- as.matrix(taxonomy.tbl)
    rownames(tax.mat) <- rsearch.obj$sequence.df$Header
    ps.obj <- phyloseq::phyloseq(phyloseq::otu_table(otu.table,
                                                     taxa_are_rows = T),
                                 phyloseq::sample_data(sample.dta),
                                 phyloseq::tax_table(tax.mat))
  } else {
    ps.obj <- phyloseq::phyloseq(phyloseq::otu_table(otu.table,
                                                     taxa_are_rows = T),
                                 phyloseq::sample_data(sample.dta))
  }
  return(ps.obj)
}

#' Convert phyloseq object to Rsearch object
#'
#' @description Creating an Rsearch object (list) from a phyloseq object.
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
#' @references
#' \url{https://joey711.github.io/phyloseq/}
#'
#' @examples
#' \dontrun{
#' # Convert phyloseq object to Rsearch object
#' rsearch_obj <- phyloseq2rsearch(phy_obj)
#'
#' # Extract read count data
#' rsearch_obj$readcount.mat
#'
#' # Extract sample data
#' rsearch_obj$sampledata.df
#'
#' # Extract sequence data
#' rsearch_obj$sequence.df
#' }
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
