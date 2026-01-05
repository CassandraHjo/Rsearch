# Taxonomic classification using the Sintax algorithm

`vs_sintax` classifies sequences using the Sintax algorithm implemented
in `VSEARCH`.

## Usage

``` r
vs_sintax(
  fasta_input,
  database,
  outfile = NULL,
  cutoff = 0,
  strand = "plus",
  sintax_random = TRUE,
  randseed = NULL,
  logfile = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fasta_input:

  (Required). A FASTA file path or a FASTA object with reads to
  classify, see *Details*.

- database:

  (Required). A FASTA file path or a FASTA object containing the
  reference database in FASTA format. The sequences need to be annotated
  with taxonomy, see *Details*.

- outfile:

  (Optional). Name of the output file. If `NULL` (default), results are
  returned as a data.frame.

- cutoff:

  (Optional). Minimum level of bootstrap support (0.0-1.0) for the
  classifications. Defaults to `0.0`.

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

- sintax_random:

  (Optional). If `TRUE` (default), the Sintax algorithm breaks ties
  between sequences with equally many kmer matches by a random draw.

- randseed:

  (Optional). Seed for the random number generator used in the Sintax
  algorithm. Defaults to `NULL`.

- logfile:

  (Optional). Name of the log file to capture messages from `VSEARCH`.
  If `NULL` (default), no log file is created.

- threads:

  (Optional). Number of computational threads to be used by `VSEARCH`.
  Defaults to `1`.

- vsearch_options:

  (Optional). A character string of additional arguments to pass to
  `VSEARCH`. Defaults to `NULL`. See *Details*.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

## Value

If `outfile` is `NULL` a data.frame is returned. If it contains a file
name (text) the data.frame is written to that file with tab-separated
columns.

The data.frame contains the classification results for each input
sequence. Both the `Header` and `Sequence` columns of `fasta_input` are
copied into this table, and in addition are also the columns for each
rank. The ranks depend on the database file used, but are typically
domain, phylum, class, order,family, genus and species. For each
classification is also a bootstrap support score. These are in separate
columns with corresponding names, i.e. domain_score, phylum_score, etc.

## Details

The sequences in the input file are classified according to the Sintax
algorithm, using `VSEARCH`, see
<https://www.biorxiv.org/content/10.1101/074161v1>.

`fasta_input` can either be a file path to a FASTA file or a FASTA
object. FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

`database` can either be a file path to a FASTA file or a FASTA object.
FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). The
`Header` texts of this file must follow the sintax-pattern, see
[`make_sintax_db`](https://cassandrahjo.github.io/Rsearch/reference/make_sintax_db.md).

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>
<https://www.biorxiv.org/content/10.1101/074161v1>

## Examples

``` r
if (FALSE) { # \dontrun{
# Example files
db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "sintax_db.fasta")
fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small.fasta")

tax.tbl <- vs_sintax(fasta_input = fasta.file, database = db.file)
View(tax.tbl)
} # }

```
