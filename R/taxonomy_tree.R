#' Make a taxonomy tree
#'
#' @description Creates a phylo object based on taxonomy
#'
#' @param taxonomy_table (Required). A data.frame with sequences and taxonomy
#' information, see \emph{Details}.
#' @param confidence (Optional). A threshold value used to replace taxa with
#' confidence scores below this to \code{NA}.
#'
#' @details In some data analyses involving OTU data a phylogenetic tree
#' describing the relatedness of the OTUs is required. To construct such trees
#' you typically need to make a multiple alignment of the sequences behind each
#' OTU, which is a huge job.
#'
#' An alternative is then to simply use the taxonomy, and create a
#' 'taxonomy-tree' instead of a phylogenetic tree. This function creates such a
#' tree from a taxonomy table of the same format as output by
#' \code{\link{vs_sintax}}.
#'
#' Distances between two OTUs reflect how high up in the taxonomy they have a
#' common taxon, i.e if they are of the same species the distance is 0, if
#' different species but same genus the distance is 1 etc. Note that \code{NA}s
#' in the taxonomy are not matched, increasing the distances, i.e if two OTUs
#' have \code{NA} as species and genus, but share family, the distance is 2.
#'
#' The \code{confidence} sets a threshold for replacing low-confidence taxa to
#' \code{NA}. For this to work the \code{taxonomy_table} must have columns with
#' such confidence scores i.e. columns domain_score, phylum_score,
#' ...species_score. If the species_score is below \code{confidence} the
#' corresponding species name is set to \code{NA}, and similar for all ranks.
#' The default is to ignore this confidence (\code{confidence = NULL}).
#'
#' From these distances a Neighbor Joining tree is built using
#' \code{\link[ape]{nj}}.
#'
#' @returns A phylo object, see \code{\link[ape]{nj}}.
#'
#' @references \url{https://www.biorxiv.org/content/10.1101/074161v1}
#'
#' @examples
#' \dontrun{
#' # Assign taxonomy with sintax
#' db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "sintax_db.fasta")
#' fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small.fasta")
#' tax.tbl <- vs_sintax(fasta_input = fasta.file, database = db.file)
#'
#' # Making tree
#' tax.tree <- taxonomy_tree(tax.tbl)
#' }
#'
#' @aliases taxonomy_tree
#'
#' @export
#'
taxonomy_tree <- function(taxonomy_table,
                          confidence = NULL){

  D.mat <- taxonomy_distance(taxonomy_table, confidence)
  tree <- ape::nj(D.mat)
  return(tree)
}

#' Creates a distance object based on taxonomy information.
#'
#' @description Creates a distance matrix based on taxonomy information
#'
#' @param taxonomy_table (Required). A data.frame with sequences and taxonomy
#' information, see \emph{Details}.
#' @param confidence (Optional). A threshold value used to replace taxa with
#' confidence scores below this to \code{NA}.
#'
#' @details In some data analyses involving OTU data, it is often useful to
#' quantify the relatedness of OTUs based on their taxonomy. This function
#' creates a distance matrix from a taxonomy table of the same format as output
#' by \code{\link{vs_sintax}}.
#'
#' Distances between two OTUs reflect how high up in the taxonomy they have a
#' common taxon, i.e if they are of the same species the distance is 0, if
#' different species but same genus the distance is 1 etc. Note that \code{NA}s
#' in the taxonomy are not matched, increasing the distances, i.e if two OTUs
#' have \code{NA} as species and genus, but share family, the distance is 2.
#'
#' The \code{confidence} sets a threshold for replacing low-confidence taxa to
#' \code{NA}. For this to work the \code{taxonomy_table} must have columns with
#' such confidence scores i.e. columns domain_score, phylum_score,
#' ...species_score. If the species_score is below \code{confidence} the
#' corresponding species name is set to \code{NA}, and similar for all ranks.
#' The default is to ignore this confidence (\code{confidence = NULL}).
#'
#' @returns A \code{dist} object containing taxonomic distances between OTUs.
#'
#' @references \url{https://www.biorxiv.org/content/10.1101/074161v1}
#'
#' @examples
#' \dontrun{
#' # Assign taxonomy with sintax
#' db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "sintax_db.fasta")
#' fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                      "small.fasta")
#' tax.tbl <- vs_sintax(fasta_input = fasta.file, database = db.file)
#'
#' # Calculate distance matrix
#' tax.dist <- taxonomy_distance(tax.tbl)
#'
#' # You can now directly use 'tax.dist' with functions like hclust or ape::nj
#' tax.tree <- ape::nj(tax.dist)
#' }
#'
#' @aliases taxonomy_distance
#'
#' @export
#'
taxonomy_distance <- function(taxonomy_table,
                              confidence = NULL){

  if(!exists("Header", where = taxonomy_table)){
    stop("The taxonomy_table must have a column named Header, with a unique text for each OTU")
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
  if(!is.null(confidence)){
    if(!exists("domain_score", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named domain")
    } else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(domain = ifelse(domain_score < confidence, NA, domain))
    }
    if(!exists("phylum", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named phylum")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(phylum = ifelse(phylum_score < confidence, NA, phylum))
    }
    if(!exists("class", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named class")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(class = ifelse(class_score < confidence, NA, class))
    }
    if(!exists("order", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named order")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(order = ifelse(order_score < confidence, NA, order))
    }
    if(!exists("family", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named family")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(family = ifelse(family_score < confidence, NA, family))
    }
    if(!exists("genus", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named genus")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(genus = ifelse(genus_score < confidence, NA, genus))
    }
    if(!exists("species", where = taxonomy_table)){
      stop("The taxonomy_table must have a column named species")
    }else {
      taxonomy_table <- taxonomy_table |>
        dplyr::mutate(species = ifelse(species_score < confidence, NA, species))
    }
  }

  tax.mat <- dplyr::select(taxonomy_table,
                           Header,
                           domain,
                           phylum,
                           class,
                           order,
                           family,
                           genus,
                           species) |>
    as.matrix()

  D.mat <- matrix(7, nrow = nrow(tax.mat), ncol = nrow(tax.mat))

  for(i in 2:8){
    idx <- which(outer(tax.mat[,i], tax.mat[,i], FUN = "=="))
    D.mat[idx] <- 8 - i
  }

  diag(D.mat) <- 0
  D.mat <- (D.mat + t(D.mat)) / 2

  rownames(D.mat) <- colnames(D.mat) <- tax.mat[,1]

  return(stats::as.dist(D.mat))
}

