# Detect chimeras by comparing sequences to a reference database

`vs_uchime_ref` detects chimeras present in the FASTA sequences in using
`VSEARCH`'s `uchime_ref` algorithm.

## Usage

``` r
vs_uchime_ref(
  fasta_input,
  database,
  nonchimeras = NULL,
  chimeras = NULL,
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

  (Required). A FASTA file path or a FASTA object with reads. See
  *Details*.

- database:

  (Required). A FASTA file path or FASTA tibble object containing the
  reference sequences. These sequences are assumed to be chimera-free.

- nonchimeras:

  (Optional). Name of the FASTA output file for the non-chimeric
  sequences. If `NULL` (default), no output is written to file.

- chimeras:

  (Optional). Name of the FASTA output file for the chimeric sequences.
  If `NULL` (default), no output is written to file.

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
  ";sample=ABC" will be added to the header. If `NULL` (default), no
  identifier is added.

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

If `nonchimeras` and `chimeras` are specified, the resulting sequences
after chimera detection written directly to the specified files in FASTA
format, and no tibbles are returned.

If `nonchimeras` and `chimeras` are `NULL`, A FASTA object containing
non-chimeric sequences with an attribute `"chimeras"` containing a
tibble of chimeric sequences is returned. If no chimeras are found, the
`"chimeras"` attribute is an empty data frame.

Additionally, the returned tibble (when applicable) has an attribute
`"statistics"` containing a tibble with chimera detection statistics.

The statistics tibble has the following columns:

- `num_nucleotides`: Total number of nucleotides used as input for
  chimera detection.

- `num_sequences`: Total number of sequences used as input for chimera
  detection.

- `min_length_input_seq`: Length of the shortest sequence used as input
  for chimera detection.

- `max_length_input_seq`: Length of the longest sequence used as input
  for chimera detection.

- `avg_length_input_seq`: Average length of the sequences used as input
  for chimera detection.

- `num_non_chimeras`: Number of non-chimeric sequences.

- `num_chimeras`: Number of chimeric sequences.

- `input`: Name of the input file/object for the chimera detection.

## Details

Chimeras in the input FASTA sequences are detected using `VSEARCH`´s
`uchime_ref`.

`fasta_input` can either be a FASTA file or a FASTA object. FASTA
objects are tibbles that contain the columns `Header` and `Sequence`,
see [`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

`database` must be a FASTA file or a FASTA object with high-quality
non-chimeric sequences.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
query_file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                        "small.fasta")
db <- file.path(file.path(path.package("Rsearch"), "extdata"),
                "sintax_db.fasta")

# Detect chimeras with default parameters and return FASTA files
vs_uchime_ref(fasta_input = query_file,
              database = db,
              nonchimeras = "nonchimeras.fa",
              chimeras = "chimeras.fa")

# Detect chimeras with default parameters and return a FASTA tibble
nonchimeras.tbl <- vs_uchime_ref(fasta_input = query_file,
                                 database = db,
                                 nonchimeras = NULL,
                                 chimeras = NULL)

# Get chimeras tibble
chimeras.tbl <- attr(nonchimeras.tbl, "chimeras")

# Get statistics tibble
statistics.tbl <- attr(nonchimeras.tbl, "statistics")
} # }
```
