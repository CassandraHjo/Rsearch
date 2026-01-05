# Join paired-end sequence reads

`vs_fastq_join` joins paired-end sequence reads into a single sequence
with a specified gap between them using `VSEARCH`.

## Usage

``` r
vs_fastq_join(
  fastq_input,
  reverse = NULL,
  output_format = "fastq",
  fastaout = NULL,
  fastqout = NULL,
  join_padgap = "NNNNNNNN",
  join_padgapq = "IIIIIIII",
  fasta_width = 0,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastq_input:

  (Required). A FASTQ file path, a FASTQ tibble object (forward reads),
  or a paired-end tibble of class `"pe_df"`. See *Details*.

- reverse:

  (Optional). A FASTQ file path or a FASTQ tibble object (reverse
  reads). Optional if `fastq_input` is a `"pe_df"` object. See
  *Details*.

- output_format:

  (Optional). Desired output format of the file or tibble: `"fasta"` or
  `"fastq"` (default).

- fastaout:

  (Optional). Name of the FASTA output file with the joined reads. If
  `NULL` (default), no output is written to a file. See *Details*.

- fastqout:

  (Optional). Name of the FASTQ output file with the joined reads. If
  `NULL` (default), no output is written to a file. See *Details*.

- join_padgap:

  (Optional). Padding sequence to use in the gap between the sequences.
  Defaults to `"NNNNNNNN"`.

- join_padgapq:

  (Optional). Quality of the padding sequence. Defaults to `"IIIIIIII"`,
  corresponding to a base quality score of 40 (a very high quality score
  with error probability `0.0001`).

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Only applies if the output file is in FASTA format. Defaults to `0`,
  which eliminates wrapping.

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

If `fastaout` or `fastqout` is specified, the joined sequences are
written to the specified output file, and no tibble is returned.

If `fastaout` or `fastqout` is `NULL`, a tibble containing the joined
reads in the format specified by `output_format` is returned.

## Details

Read pairs from the input FASTQ files (`fastq_input` and `reverse`) are
joined into a single sequence by adding a gap with a specified padding
sequence. The resulting sequences consist of the forward read, the
padding sequence, and the reverse complement of the reverse read.

`fastq_input` and `reverse` can either be file paths to FASTQ files or
FASTQ objects. FASTQ objects are tibbles that contain the columns
`Header`, `Sequence`, and `Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html). Forward
and reverse reads must appear in the same order and have the same total
number of reads in both files.

If `fastq_input` is an object of class `"pe_df"`, the reverse reads are
automatically extracted from its `"reverse"` attribute unless explicitly
provided via the `reverse` argument. This simplifies function calls when
using paired-end tibbles created by functions such as
[`fastx_synchronize`](https://cassandrahjo.github.io/Rsearch/reference/fastx_synchronize.md)
or
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md).

If `fastaout` or `fastqout` is specified, the joined reads are written
to the respective file in either FASTA or FASTQ format.

If both `fastaout` or `fastqout` are `NULL`, the results are returned as
a FASTA or FASTQ object, and no file is written. `output_format` must
match the desired output files/objects.

Any input sequence with fewer bases than the value set in `minlen` is
discarded. By default, `minlen` is set to 0, which means that no
sequences are removed. However, using the default value may allow empty
sequences to remain in the results.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
fastq_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
                         "small_R1.fq")
reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R2.fq")
output_format <- "fastq"

# Execute joining and return a FASTQ tibble
join_seqs <- vs_fastq_join(fastq_input = fastq_input,
                           reverse = reverse,
                           output_format = output_format)

# Execute joining and write joined sequences to file
vs_fastq_join(fastq_input = fastq_input,
              reverse = reverse,
              fastqout = "joined_sequences.fq",
              output_format = output_format)

} # }
```
