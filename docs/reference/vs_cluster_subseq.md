# Cluster FASTA sequences

`vs_cluster_subseq` clusters FASTA sequences from a given file or object
using `VSEARCH`´s `cluster_fast` method and 100 identity. The function
automatically sorts sequences by decreasing length before clustering.

## Usage

``` r
vs_cluster_subseq(
  fasta_input,
  centroids = NULL,
  strand = "plus",
  sizein = TRUE,
  fasta_width = 0,
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

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

- sizein:

  (Optional). If `TRUE` (default), abundance annotations present in
  sequence headers are taken into account.

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Defaults to `0`, which eliminates wrapping.

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

If `centroids` is not specified, a FASTA object is returned. This is a
`tibble` with columns `Header` and `Sequence`, and also the additional
column(s) `members` and, if `sizein = TRUE`, `size`.

## Details

After merging/dereplication some sequences may be sub-sequences of
longer sequences. This function will cluster such sequences at 100
(terminal gaps ignored), and keep the longest in each cluster as the
centroid.

`fasta_input` can either be a file path to a FASTA file or a FASTA
object. FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

If `sizein = TRUE` (default) the FASTA headers must contain text
matching the regular expression `"size=[0-9]+"` indicating the copy
number (=size) of each input sequence. This is then summed for each
cluster and added to the output. This text is typically added by
de-replication, see
[`vs_fastx_uniques`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.md).

The number of distinct sequences in each cluster is output as `members`.

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

# De-replicating
derep.tbl <- vs_fastx_uniques(fasta_input, output_format = "fasta")

# Clustering subsequences
cluster.tbl <- vs_cluster_subseq(fasta_input = derep.tbl)

# Cluster sequences and write centroids to a file
vs_cluster_subseq(fasta_input = derep.tbl,
                  centroids = "distinct.fa")
} # }
```
