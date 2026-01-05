# Detect chimeras without external references (i.e. de novo)

`vs_uchime_denovo` detects chimeras present in the FASTA sequences in
using `VSEARCH`'s `uchime_denovo` algorithm. Automatically sorts
sequences by decreasing abundance to enhance chimera detection accuracy.

## Usage

``` r
vs_uchime_denovo(
  fasta_input,
  nonchimeras = NULL,
  chimeras = NULL,
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fasta_width = 0,
  sample = NULL,
  log_file = NULL,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fasta_input:

  (Required). A FASTA file path or a FASTA object with reads. If a
  tibble is provided, any columns in addition to `Header` and `Sequence`
  will be preserved in the output. See *Details*.

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

If `nonchimeras` and `chimeras` are `NULL`, a FASTA object containing
non-chimeric sequences is returned. This output tibble will include any
additional columns that were present in the `fasta_input` tibble. An
attribute named `"chimeras"` will contain a tibble of the chimeric
sequences, also with the additional columns preserved.

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
`uchime_denovo`. In de novo mode, input FASTA file/object must present
abundance annotations (i.e. a pattern \[;\]size=integer\[;\] in the
header). Input order matters for chimera detection, so it is recommended
to sort sequences by decreasing abundance.

`fasta_input` can either be a FASTA file or a FASTA object. FASTA
objects are tibbles that contain the columns `Header` and `Sequence`,
see [`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html).

When providing a tibble as `fasta_input`, you can include additional
columns with metadata (e.g., OTU IDs, sample origins). The function will
preserve these columns by joining them back to the results based on the
DNA sequence. This allows you to keep your metadata associated with your
sequences throughout the chimera detection process.

If `nonchimeras` and `chimeras` are specified, resulting non-chimeric
and chimeric sequences are written to these files in FASTA format.

If `nonchimeras` and `chimeras` are `NULL`, results are returned as a
FASTA-objects.

`nonchimeras` and `chimeras` must either both be specified or both be
`NULL`.

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
fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
                         "small_R1.fq")
nonchimeras <- "nonchimeras.fa"
chimeras <- "chimeras.fa"

# Detect chimeras with default parameters and return FASTA files
vs_uchime_denovo(fasta_input = fasta_input,
                 nonchimeras = nonchimeras,
                 chimeras = chimeras)

# Detect chimeras with default parameters and return a FASTA tibble
nonchimeras.tbl <- vs_uchime_denovo(fasta_input = fasta_input,
                                    nonchimeras = NULL,
                                    chimeras = NULL)

# Get chimeras tibble
chimeras.tbl <- attr(nonchimeras.tbl, "chimeras")

# Get statistics tibble
statistics.tbl <- attr(nonchimeras.tbl, "statistics")
} # }
```
