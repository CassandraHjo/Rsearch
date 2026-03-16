# Cluster FASTA sequences

`vs_cluster_size` clusters FASTA sequences from a given file or object
using `VSEARCH`´s `cluster_size` method. The function automatically
sorts sequences by decreasing abundance before clustering.

## Usage

``` r
vs_cluster_size(
  fasta_input,
  centroids = NULL,
  otutabout = NULL,
  uc = NULL,
  clusters = NULL,
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
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fasta_input:

  (Required). A FASTA file path or a FASTA object containing reads to
  cluster. See *Details*.

- centroids:

  (Optional). A character string specifying the name of the FASTA output
  file for the cluster centroid sequences. If `NULL` (default), no
  output is written to a file and the centroid sequences are returned as
  a FASTA object. See *Details*.

- otutabout:

  (Optional). A character string specifying the name of the output file
  in an OTU table format. If `NULL` (default), no output is written to a
  file. If `TRUE`, the output is returned as a tibble. See *Details*.

- uc:

  (Optional). A character string specifying the name of the output file
  in a uclust-like format. If `NULL` (default), no output is written to
  a file. If `TRUE`, the output is returned as a tibble. See *Details*.

- clusters:

  (Optional). Output each cluster to a separate FASTA file using the
  prefix string provided and a ticker (0, 1, 2, etc.) to construct the
  path and the filenames. Defaults to `NULL`, meaning no such files are
  created.

- size_column:

  (Optional). If `TRUE`, a column with the size of each centroid is
  added to the centroid output tibble.

- id:

  (Optional). Pairwise identity threshold for sequence to be added to a
  cluster. Defaults to `0.97`. See *Details*.

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

- sizein:

  (Optional). If `TRUE` (default), abundance annotations present in
  sequence headers are taken into account.

- sizeout:

  (Optional). If `TRUE` (default), abundance annotations are added to
  FASTA headers.

- relabel:

  (Optional). Relabel sequences using the given prefix and a ticker to
  construct new headers. Defaults to `NULL`.

- relabel_sha1:

  (Optional). If `TRUE` (default), relabel sequences using the SHA1
  message digest algorithm. Defaults to `FALSE`.

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Defaults to `0`, which eliminates wrapping.

- sample:

  (Optional). Add the given sample identifier string to sequence
  headers. For instance, if the given string is "ABC", the text
  ";sample=ABC" will be added to the header. This option is only
  applicable when the output format is FASTA (`centroids`). If `NULL`
  (default), no identifier is added.

- log_file:

  (Optional). Name of the log file to capture messages from `VSEARCH`.
  If `NULL` (default), no log file is created.

- threads:

  (Optional). Number of computational threads to be used by `VSEARCH`.
  Defaults to `1`.

- vsearch_options:

  (Optional). Additional arguments to pass to `VSEARCH`. Defaults to
  `NULL`. See *Details*.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

## Value

A tibble or `NULL`.

If `centroids` is specified the centroid sequences are written to the
specified file, and no tibble is returned.

If `otutabout` is `TRUE`, an OTU table is returned as a tibble. If
`otutabout` is a character string, the output is written to the file,
and no tibble is returned.

If neither `centroids` nor `otutabout` is specified, a FASTA object with
the centroid sequences and additional column `otu_id` is returned. The
clustering statistics are included as an attribute named `"statistics"`.

The `"statistics"` attribute of the returned tibble (when `centroids` is
`NULL`) is a tibble with the following columns:

- `num_nucleotides`: Total number of nucleotides used as input for
  clustering.

- `min_length_input_seq`: Length of the shortest sequence used as input
  for clustering.

- `max_length_input_seq`: Length of the longest sequence used as input
  for clustering.

- `avg_length_input_seq`: Average length of the sequences used as input
  for clustering.

- `num_clusters`: Number of clusters generated.

- `min_size_cluster`: Size of the smallest cluster.

- `max_size_cluster`: Size of the largest cluster.

- `avg_size_cluster`: Average size of the clusters.

- `num_singletons`: Number of singletons after clustering.

- `input`: Name of the input file/object for the clustering.

## Details

Sequences are clustered based on the pairwise identity threshold
specified by `id`. Sequences are sorted by decreasing abundance before
clustering. The centroid of each cluster is the first sequence added to
the cluster.

`fasta_input` can either be a file path to a FASTA file or a FASTA
object. FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

If neither `centroids`, `otutabout`, nor `uc` is specified (default),
the function returns the centroid sequences as a FASTA object with an
additional column `otu_id`. This column contains the identifier
extracted from each sequence header.

If `centroids` is specified, centroid sequences are written to the
specified file in FASTA format.

`otutabout` gives the option to output the results in an OTU table
format with tab-separated columns. When writing to a file, the first
line starts with the string "#OTU ID", followed by a tab-separated list
of all sample identifiers (formatted as "sample=X"). Each subsequent
line, corresponding to an OTU, begins with the OTU identifier and is
followed by tab-separated abundances for that OTU in each sample. If
`otutabout` is a character string, the output is written to the
specified file. If `otutabout` is `TRUE`, the function returns the OTU
table as a tibble, where the first column is named `otu_id` instead of
"#OTU ID".

`uc` gives the option to output the results in a uclust-like format.
This is a tab-separated uclust-like format with 10 columns and 3
different types of entries (S, H or C). Each FASTA sequence i the input
can either be a cluster centroid (S) or a hit (H) assigned to a cluster.
Cluster records (C) summarize information (size, centroid label) for
each cluster. See the '–uc' section in the `VSEARCH` manual for more
information.

`id` is a value between 0 and 1 that defines the minimum pairwise
identity required for a sequence to be added to a cluster. A sequence is
not added to a cluster if its pairwise identity with the centroid is
below the `id` threshold. Pairwise identity is calculated as the number
of matching columns divided by the alignment length minus terminal gaps.

If `log_file` is `NULL` and `centroids` is specified, clustering
statistics from `VSEARCH` will not be captured.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
                                   "small.fasta")
centroids <- NULL

# Cluster sequences and return a FASTA tibble
cluster_seqs <- vs_cluster_size(fasta_input = fasta_input,
                                centroids = centroids)

# Extract clustering statistics
statistics <- attr(cluster_seqs, "statistics")

# Cluster sequences and write centroids to a file
vs_cluster_size(fasta_input = fasta_input,
                centroids = "centroids_sequences.fa")
} # }
```
