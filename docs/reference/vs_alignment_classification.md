# Taxonomic classification with LCA

`vs_alignment_classification` assigns taxonomy by global alignment and
Last Common Ancestor (LCA) consensus of database hits using `VSEARCH`.

## Usage

``` r
vs_alignment_classification(
  fastx_input,
  database,
  lcaout = NULL,
  lca_cutoff = 1,
  top_hits_only = FALSE,
  gapopen = "20I/2E",
  gapext = "2I/1E",
  id = 0.7,
  strand = "plus",
  maxaccepts = 2,
  maxrejects = 32,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ object. See
  *Details*.

- database:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ tibble object
  containing the target sequences.

- lcaout:

  (Optional). A character string specifying the name of the output file.
  If `NULL` (default), no output is written to a file and the results
  are returned as a tibble with the columns `query_id` and `taxonomy`.

- lca_cutoff:

  (Optional). Adjust the fraction of matching hits required for the last
  common ancestor (LCA). Defaults to `1.0`, which requires all hits to
  match at each taxonomic rank for that rank to be included. If a lower
  cutoff value is used, e.g. 0.95, a small fraction of non-matching hits
  are allowed while that rank will still be reported. The argument to
  this option must be between `0.5` and `1.0`.

- top_hits_only:

  (Optional). If `TRUE`, only the top hits with an equally high
  percentage of identity between the query and database sequence sets
  are written to the output. Defaults to `FALSE`.

- gapopen:

  (Optional). Penalties for gap opening. Defaults to `"20I/2E"`. See
  *Details*.

- gapext:

  (Optional). Penalties for gap extension. Defaults to `"2I/1E"`. See
  *Details*.

- id:

  (Optional). Pairwise identity threshold. Defines the minimum identity
  required for matches. Defaults to `0.7`.

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

- maxaccepts:

  (Optional). Maximum number of matching target sequences to accept
  before stopping the search for a given query. Defaults to `2`. Must be
  larger than `1` for information to be useful.

- maxrejects:

  (Optional). Maximum number of non-matching target sequences to
  consider before stopping the search for a given query. Defaults to 32.
  If `maxaccepts` and `maxrejects` are both set to 0, the complete
  database is searched.

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

If `lcaout` is specified the results are written to the specified file.
If `lcaout` is `NULL` a data.frame is returned.

The data.frame contains the classification results for each query
sequence. Both the `Header` and `Sequence` columns of `fasta_input` are
copied into this table, and in addition are also the columns for each
rank. The ranks depend on the database file used, but are typically
domain, phylum, class, order,family, genus and species.

## Details

Performs global sequence alignment against a reference database and
assigns taxonomy using the Last Common Ancestor (LCA) approach,
reporting the deepest taxonomic level consistently supported by the
majority of hits.

`fastx_input` and `database` can either be file paths to a FASTA/FASTQ
files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
columns `Header` and `Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

Pairwise identity (`id`) is calculated as the number of matching columns
divided by the alignment length minus terminal gaps.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

Visit the `VSEARCH`
[documentation](https://github.com/torognes/vsearch?tab=readme-ov-file#getting-help)
for information about defining `gapopen` and `gapext`.

## References

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Example files
db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "sintax_db.fasta")
fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small.fasta")

tax.tbl <- vs_alignment_classification(fastx_input = fasta.file,
                                       database = db.file)
View(tax.tbl)
} # }
```
