# Subsample sequences

`vs_fastx_subsample` subsamples sequences in FASTA/FASTQ file or object
by randomly extracting sequences based on number or percentage using
`VSEARCH`.

## Usage

``` r
vs_fastx_subsample(
  fastx_input,
  output_format = "fastq",
  fastx_output = NULL,
  sample_pct = NULL,
  sample_size = NULL,
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  randseed = NULL,
  fasta_width = 0,
  sample = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ object. See
  *Details*.

- output_format:

  (Optional). Desired output format of file or tibble: `"fasta"` or
  `"fastq"` (default). If `fastx_input` is a FASTA file path or a FASTA
  object, `output_format` cannot be `"fastq"`.

- fastx_output:

  (Optional). Name of the output file for subsampled reads from
  `fastx_input`. File can be in either FASTA or FASTQ format, depending
  on `output_format`. If `NULL` (default), no sequences are written to
  file. See *Details*.

- sample_pct:

  (Optional). Percentage of the input sequences to be subsampled.
  Numeric value ranging from `0.0` to `100.0`. Defaults to `NULL`.

- sample_size:

  (Optional). The given number of sequences to extract. Must be a
  positive integer if specified. Defaults to `NULL`.

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

- randseed:

  (Optional). Random seed. Must be a positive integer. A given seed
  always produces the same output, which is useful for replicability.
  Defaults to `NULL`.

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Defaults to `0`, which eliminates wrapping.

- sample:

  (Optional). Add the given sample identifier string to sequence
  headers. For instance, if the given string is "ABC", the text
  ";sample=ABC" will be added to the header. If `NULL` (default), no
  identifier is added.

- threads:

  (Optional). Number of computational threads to be used by
  `VSEARCH`.Defaults to `1`.

- vsearch_options:

  Additional arguments to pass to `VSEARCH`. Defaults to `NULL`. See
  *Details*.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

## Value

A tibble or `NULL`.

If `fastx_output` is specified, the subsampled sequences are written to
the specified output file, and no tibble is returned.

If `fastx_output` `NULL`, a tibble containing the subsampled reads in
the format specified by `output_format` is returned.

## Details

Sequences in the input file/object (`fastx_input`) are subsampled by
randomly extracting a specified number or percentage of sequences.
Extraction is performed as random sampling with a uniform distribution
among the input sequences and without replacement.

`fastx_input` can either be a FASTA/FASTQ file or a FASTA/FASTQ object.
FASTA objects are tibbles that contain the columns `Header` and
`Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

Specify either `sample_size` or `sample_pct` to determine the number or
percentage of sequences to subsample. Only one of these parameters can
be specified at a time. If neither is specified, an error is thrown.

If `fastx_output` is specified, the sampled sequences are output to this
file in format given by `output_format`. If `fastx_output` is `NULL`,
the sample sequences are returned as a FASTA or FASTQ object, depending
on `output_format`.

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
sample_size <- 10

# Subsample sequences and return a FASTQ tibble
subsample_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   output_format = output_format,
                                   sample_size = sample_size)

# Subsample sequences and write subsampled sequences to a file
vs_fastx_subsample(fastx_input = fastx_input,
                   fastx_output = "subsample.fq",
                   output_format = output_format,
                   sample_size = sample_size)
} # }
```
