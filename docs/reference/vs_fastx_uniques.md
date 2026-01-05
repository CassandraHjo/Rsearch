# Dereplicate sequences

`vs_fastx_uniques` performs dereplication of sequences in a FASTA/FASTQ
file or object by merging identical sequences using `VSEARCH`.

## Usage

``` r
vs_fastx_uniques(
  fastx_input,
  output_format = "fasta",
  fastx_output = NULL,
  minuniquesize = 1,
  strand = "plus",
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fastq_qout_max = FALSE,
  fasta_width = 0,
  sample = NULL,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ object. See
  *Details*.

- output_format:

  (Optional). Desired output format of file or tibble: `"fasta"`
  (default) or `"fastq"`. If `fastx_input` is a FASTA file path or a
  FASTA object, `output_format` cannot be `"fastq"`.

- fastx_output:

  (Optional). Name of the output file for dereplicated reads from
  `fastx_input`. File can be in either FASTA or FASTQ format, depending
  on `output_format`. If `NULL` (default), no sequences are written to
  file. See *Details*.

- minuniquesize:

  (Optional). Minimum abundance value post-dereplication for a sequence
  not to be discarded. Defaults to `1`.

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

- fastq_qout_max:

  (Optional). If `TRUE`, the quality score will be the highest (best)
  quality score observed in each position. Defaults to `FALSE`.

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Defaults to `0`, which eliminates wrapping.

- sample:

  (Optional). Add the given sample identifier string to sequence
  headers. For nstance, if the given string is "ABC", the text
  ";sample=ABC" will be added to the header. If `NULL` (default), no
  identifier is added.

- vsearch_options:

  (Optional). A character string of additional arguments to pass to
  `VSEARCH`. Defaults to `NULL`. See *Details*.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

## Value

A tibble or `NULL`.

If `fastx_output` is specified, the dereplicated sequences are written
to the specified output file, and no tibble is returned.

If `fastx_output` `NULL`, a tibble containing the dereplicated reads in
the format specified by `output_format` is returned.

## Details

Sequences in the input file/object (`fastx_input`) are dereplicated by
merging identical sequences. Identical sequences are defined as
sequences with the same length and the same string of nucleotides (case
insensitive, T and U are considered the same).

`fastx_input` can either be a FASTA/FASTQ file or a FASTA/FASTQ object.
FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

By default, the quality scores in FASTQ output files will correspond to
the average error probability of the nucleotides in the each position.
If `fastq_qout_max = TRUE`, the quality score will be the highest (best)
quality score observed in each position.

If `fastx_output` is specified, the dereplicated sequences are output to
this file in format given by `output_format`. If `fastx_output` is
`NULL`, the dereplicated sequences are returned as a FASTA or FASTQ
object, depending on `output_format`.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
                         "small_R1.fq")
fastx_output <- NULL
output_format <- "fastq"

# Dereplicate sequences and return a FASTQ tibble
derep_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
                             fastx_output = fastx_output,
                             output_format = output_format)

# Dereplicate sequences and write derelicated sequences to a file
vs_fastx_uniques(fastx_input = fastx_input,
                 fastx_output = "dereplicated_sequences.fq",
                 output_format = output_format)
} # }
```
