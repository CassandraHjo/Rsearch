# Denoising FASTA sequences

`vs_cluster_unoise` performs denoising of FASTA sequences from a given
file or object using `VSEARCH`´s `cluster_unoise` method.

## Usage

``` r
vs_cluster_unoise(
  fasta_input,
  otutabout = NULL,
  clusters = NULL,
  minsize = 8,
  unoise_alpha = 2,
  relabel = NULL,
  relabel_sha1 = FALSE,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fasta_input:

  (Required). A FASTA file path or a FASTA object containing reads to
  denoise. See *Details*.

- otutabout:

  (Optional). A character string specifying the name of the output file
  in an OTU table format. If `NULL` (default), the output is returned as
  a tibble in R. See *Details*.

- clusters:

  (Optional). Output each cluster to a separate FASTA file using the
  prefix string provided and a ticker (0, 1, 2, etc.) to construct the
  path and the filenames. Defaults to `NULL`, meaning no such files are
  created.

- minsize:

  (Optional). Minimum abundance of cluster centroids. Defaults to `8`.

- unoise_alpha:

  (Optional). Alpha value for the UNOISE algorithm. Defaults to `2`.

- relabel:

  (Optional). Relabel sequences using the given prefix and a ticker to
  construct new headers. Defaults to `NULL`.

- relabel_sha1:

  (Optional). If `TRUE` (default), relabel sequences using the SHA1
  message digest algorithm. Defaults to `FALSE`.

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

A read count table with one row for each cluster and one column for each
sample. If `otutabout` is a text it is assumed to be a file name, and
the results are written to this file. If no such text is supplied
(default), it is returned as a tibble.

The first two columns of this tibble lists the `Header` and `Sequence`
of the centroid sequences for each cluster.

The clustering statistics are included as an attribute named
`"statistics"` with the following columns:

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

Sequences are denoised according to the UNOISE version 3 algorithm by
Robert Edgar, but without the de novo chimera removal step. In this
algorithm, clustering of sequences depends both on their similarity and
their abundances. The abundance ratio (skew) is the abundance of a new
sequence divided by the abundance of the centroid sequence. This skew
must not be larger than beta if the sequences should be clustered
together. Beta is calculated as 2 raised to the power of minus 1 minus
alpha times the sequence distance. The sequence distance used is the
number of mismatches in the alignment, ignoring gaps. This means that
the abundance must be exponentially lower as the distance increases from
the centroid for a new sequence to be included in the cluster.

The argument `minsize` will affect the total number of clusters,
specifying the minimum copy number required for any centroid. A larger
value means (in general) fewer clusters.

`fasta_input` can either be a file path to a FASTA file or a FASTA
object. FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

The `Header` column **must** contain the size (copy number) for each
read. The size information must have the format ";size=X", where X is
the count for the given sequence. This is obtained by running all reads
through
[`vs_fastx_uniques`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.md)
with `sizeout = TRUE`.

You may use reads for a single sample or all reads from all samples as
input. In the latter case the `Header` must also contain sample
information on the format ";sample=xxx" where "xxx" is a unique sample
identifier text. Again, this is obtained by using
[`vs_fastx_uniques`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.md)
on the reads for each sample prior to this step. Use the
`sample = "xxx"` argument, where "xxx" is replaced with some unique text
for each sample.

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
# A small fasta file
fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "small.fasta")

# Denoise sequences and read counts
denoised.tbl <- vs_cluster_unoise(fasta_input = fasta_input)
head(denoised.tbl)

# Extract clustering statistics
statistics <- attr(denoised.tbl, "statistics")

# Cluster sequences and write results to a file
vs_cluster_unoise(fasta_input = fasta_input,
                  otutabout = "otutable.tsv")
} # }
```
