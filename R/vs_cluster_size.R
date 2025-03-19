#' Cluster FASTA sequences
#'
#' @description \code{vs_cluster_size} clusters FASTA sequences from a given
#' file or object using \code{VSEARCH}´s \code{cluster_size} method. The
#' function automatically sorts sequences by decreasing abundance before
#' clustering.
#'
#' @param fasta_input A FASTA file path or a FASTA object containing reads to
#' cluster. See \emph{Details}.
#' @param centroids A character string specifying the name of the FASTA output
#' file for the cluster centroid sequences. If \code{NULL} (default), no output
#' is written to a file and the centroid sequences are returned as a FASTA
#' object. See \emph{Details}.
#' @param otutabout A character string specifying the name of the output file in
#' an OTU table format. If \code{NULL} (default), no output is written to a file.
#' If \code{TRUE}, the output is returned as a tibble. See \emph{Details}.
#' @param size_column If \code{TRUE}, a column with the size of each centroid is
#' added to the centroid output tibble.
#' @param id Pairwise identity threshold for sequence to be added to a
#' cluster. Defaults to \code{0.97}. See \emph{Details}.
#' @param strand Specifies which strand to consider when comparing sequences.
#' Can be either \code{"plus"} (default) or \code{"both"}.
#' @param sizein If \code{TRUE} (default), abundance annotations present in
#' sequence headers are taken into account.
#' @param sizeout If \code{TRUE} (default), abundance annotations are added to
#' FASTA headers.
#' @param relabel Relabel sequences using the given prefix and a ticker to
#' construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 If \code{TRUE} (default), relabel sequences using the
#' SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param fasta_width Number of characters per line in the output FASTA
#' file. Defaults to \code{0}, which eliminates wrapping.
#' @param sample Add the given sample identifier string to sequence headers. For
#' instance, if the given string is "ABC", the text ";sample=ABC" will be added
#' to the header. If \code{NULL} (default), no identifier is added.
#' @param log_file Name of the log file to capture messages from \code{VSEARCH}.
#' If \code{NULL} (default), no log file is created.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param vsearch_options Additional arguments to pass to \code{VSEARCH}.
#' Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Sequences are clustered based on the pairwise identity threshold specified by
#' \code{id}. Sequences are sorted by decreasing abundance before clustering.
#' The centroid of each cluster is the first sequence added to the cluster.
#'
#' \code{fasta_input} can either be a file path to a FASTA file or a FASTA
#' object. FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}.
#'
#' If neither \code{centroids} nor \code{otutabout} is specified (default), the
#' function returns the centroid sequences as a FASTA object.
#'
#' If \code{centroids} is specified, centroid sequences are written to the
#' specified file in FASTA format.
#'
#' \code{otutabout} gives the option to output the clustering results in an OTU
#' table format with tab-separated columns. The first line will start with
#' the string ’#OTU ID’ and is followed by a tab-separated list of all sample
#' identifiers ("sample=X"). The following lines, one for each OTU, starts with
#' the OTU identifier and is followed by a tab-separated list of abundances for
#' that OTU in each sample, in the order given on the first line. The OTU and
#' sample identifiers are extracted from the FASTA headers of the sequences.
#' If \code{otutabout} is a character string, the output is written to the
#' specified file. If \code{otutabout} is \code{TRUE}, the function returns
#' the OTU table as a tibble.
#'
#' \code{id} is a value between 0 and 1 that defines the minimum pairwise
#' identity required for a sequence to be added to a cluster. A sequence is not
#' added to a cluster if its pairwise identity with the centroid is bellow the
#' \code{id} threshold.
#' Pairwise identity is calculated as the number of matching columns divided by
#' the alignment length minus terminal gaps.
#'
#' If \code{log_file} is \code{NULL} and \code{centroids} is specified,
#' clustering statistics from \code{VSEARCH} will not be captured.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{centroids} is specified the centroid sequences are written to the
#' specified file, and no tibble is returned.
#'
#' If \code{otutabout} is \code{TRUE}, an OTU table is returned as a tibble.
#' If \code{otutabout} is a character string, the output is written to the file,
#' and no tibble is returned.
#'
#' If neither \code{centroids} nor \code{otutabout} is specified, a FASTA object
#' containing the centroid sequences is returned. The clustering statistics are
#' included as an attribute named \code{"statistics"}.
#'
#' The \code{"statistics"} attribute of the returned tibble (when
#' \code{centroids} is \code{NULL}) is a tibble with the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: Total number of nucleotides used as input for
#'   clustering.
#'   \item \code{min_length_input_seq}: Length of the shortest sequence used as
#'   input for clustering.
#'   \item \code{max_length_input_seq}: Length of the longest sequence used as
#'   input for clustering.
#'   \item \code{avg_length_input_seq}: Average length of the sequences used as
#'   input for clustering.
#'   \item \code{num_clusters}: Number of clusters generated.
#'   \item \code{min_size_cluster}: Size of the smallest cluster.
#'   \item \code{max_size_cluster}: Size of the largest cluster.
#'   \item \code{avg_size_cluster}: Average size of the clusters.
#'   \item \code{num_singletons}: Number of singletons after clustering.
#'   \item \code{input}: Name of the input file/object for the clustering.
#' }
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                                    "R1_sample1_small.fa")
#' centroids <- NULL
#'
#' # Cluster sequences and return a FASTA tibble
#' cluster_seqs <- vs_cluster_size(fasta_input = fasta_input,
#'                                 centroids = centroids)
#'
#' # Extract clustering statistics
#' statistics <- attr(cluster_seqs, "statistics")
#'
#' # Cluster sequences and write centroids to a file
#' vs_cluster_size(fasta_input = fasta_input,
#'                               centroids = "centroids_sequences.fa")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_cluster_size vs_cluster cluster_size cluster
#'
#' @export
#'
vs_cluster_size <- function(fasta_input,
                            centroids = NULL,
                            otutabout = NULL,
                            size_column = FALSE,
                            id = 0.97,
                            strand = "plus",
                            sizein = TRUE,
                            sizeout = TRUE,
                            relabel = NULL,
                            relabel_sha1 = FALSE,
                            fasta_width = 0,
                            sample = NULL,
                            log_file = NULL,
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
  if (!is.null(centroids) && !is.null(otutabout)) {
    stop("Only one of 'centroids' or 'otutabout' can be specified.")
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
  }

  # Check is input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)

  # Determine output file based on user input
  if (!is.null(centroids)) {
    outfile <- ifelse(is.character(centroids), centroids, tempfile(pattern = "centroids", fileext = ".fa"))
  } else if (!is.null(otutabout)) {
    outfile <- ifelse(is.character(otutabout), otutabout, tempfile(pattern = "otutable", fileext = ".tsv"))
  } else {
    outfile <- tempfile(pattern = "centroids", fileext = ".fa")
  }

  # Only add temporary files to temp_files
  if (is.null(centroids) && (is.null(otutabout) || !is.character(otutabout))) {
    temp_files <- c(temp_files, outfile)
  }

  # Build argument string for command line
  args <- c("--cluster_size", shQuote(fasta_file),
            "--id", id,
            "--threads", 1,
            "--strand", strand,
            "--fasta_width", fasta_width)

  if (!is.null(centroids)) {
    args <- c(args, "--centroids", outfile)
  } else if (!is.null(otutabout)) {
    args <- c(args, "--otutabout", outfile)
  } else {
    args <- c(args, "--centroids", outfile) # Default output
  }

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  # Add relabeling arguments if specified
  if (!is.null(relabel)){
    args <- c(args, "--relabel", relabel)
  }

  if (relabel_sha1){
    args <- c(args, "--relabel_sha1", "")
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add sample identifier if specified
  if (!is.null(sample)) {
    args <- c(args, "--sample", sample)
  }

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Determine return output
  if (!is.null(centroids)) {
    return(invisible(NULL)) # No return if centroids is specified
  } else if (!is.null(otutabout)) {
    if (is.character(otutabout)) {
      return(invisible(NULL)) # File output only
    } else {
      return(readr::read_delim(outfile)) # Return as tibble
    }
  } else {
    centroids_fasta <- microseq::readFasta(outfile)

    if (size_column) {
      centroids_fasta <- centroids_fasta |>
        dplyr::mutate(centroid_size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
        dplyr::mutate(centroid_size = as.numeric(centroid_size)) |>
        dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))
    }

    statistics <- calculate_cluster_statistics(centroids_fasta,
                                               fasta_file,
                                               fasta_input_name)
    attr(centroids_fasta, "statistics") <- statistics
    return(centroids_fasta)
  }
}

#' Calculate clustering statistics
#'
#' @description \code{calculate_cluster_statistics} calculates important
#' clustering statistics after running \code{vs_cluster_size}, including the
#' number of clusters, sequences, and nucleotides, as well as the lengths and
#' sizes of the sequences and clusters.
#'
#' @param centroids_fasta Output tibble from clustering with centroids.
#' Contains the columns: Header, Sequence, and centroid_size (if
#' \code{size_column} is specified in \code{vs_cluster_size}).
#' @param fasta_file File path to FASTA containing the input sequences to the
#' clustering.
#' @param fasta_input_name Name of the file/object with the input sequences
#' that was used in the clustering.
#'
#' @return A tibble with clustering statistics, including:
#' \itemize{
#'   \item \code{num_nucleotides}: Total number of nucleotides used as input for
#'   clustering.
#'   \item \code{min_length_input_seq}: Length of the shortest sequence used as
#'   input for clustering.
#'   \item \code{max_length_input_seq}: Length of the longest sequence used as
#'   input for clustering.
#'   \item \code{avg_length_input_seq}: Average length of the sequences used as
#'   input for clustering.
#'   \item \code{num_clusters}: Number of clusters generated.
#'   \item \code{min_size_cluster}: Size of the smallest cluster.
#'   \item \code{max_size_cluster}: Size of the largest cluster.
#'   \item \code{avg_size_cluster}: Average size of the clusters.
#'   \item \code{num_singletons}: Number of singletons after clustering.
#'   \item \code{input}: Name of the input file/object for the clustering.
#' }
#'
#' @return A tibble with clustering statistics.
#'
#' @noRd
#'
calculate_cluster_statistics <- function(centroids_fasta,
                                         fasta_file,
                                         fasta_input_name) {

  # Process clustering output
  if (!"centroid_size" %in% colnames(centroids_fasta)) {
    centroids_fasta <- centroids_fasta |>
      dplyr::mutate(centroid_size = stringr::str_extract(Header, "(?<=;size=)\\d+")) |>
      dplyr::mutate(centroid_size = as.numeric(centroid_size)) |>
      dplyr::mutate(Header = stringr::str_remove(Header, ";size=\\d+"))
  }
  # Make tibble from input sequences to the clustering
  input.df <- microseq::readFasta(fasta_file)

  # Calculate statistics
  num_nucleotides <- sum(nchar(input.df$Sequence))
  min_length_input_seq <- min(nchar(input.df$Sequence))
  max_length_input_seq <- max(nchar(input.df$Sequence))
  avg_length_input_seq <- mean(nchar(input.df$Sequence))
  num_clusters <- nrow(centroids_fasta)
  min_size_cluster <- min(centroids_fasta$centroid_size)
  max_size_cluster <- max(centroids_fasta$centroid_size)
  avg_size_cluster <- round(mean(centroids_fasta$centroid_size), 1)
  num_singletons <- sum(centroids_fasta$centroid_size == 1)

  # Create table
  result_table <- data.frame(
    num_nucleotides = num_nucleotides,
    min_length_input_seq = min_length_input_seq,
    max_length_input_seq = max_length_input_seq,
    avg_length_input_seq = avg_length_input_seq,
    num_clusters = num_clusters,
    min_size_cluster = min_size_cluster,
    max_size_cluster = max_size_cluster,
    avg_size_cluster = avg_size_cluster,
    num_singletons = num_singletons,
    input = fasta_input_name
  )

  return(result_table)
}

