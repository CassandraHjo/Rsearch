# Trim and/or filter sequences in FASTA/FASTQ format

`vs_fastx_trim_filt` trims and/or filters FASTA/FASTQ sequences using
`VSEARCH`. This function processes both forward and reverse reads (if
provided) and allows for various filtering criteria based on sequence
quality, length, abundance, and more.

## Usage

``` r
vs_fastx_trim_filt(
  fastx_input,
  reverse = NULL,
  output_format = "fastq",
  fastaout = NULL,
  fastqout = NULL,
  fastaout_rev = NULL,
  fastqout_rev = NULL,
  trunclen = NULL,
  truncqual = 1,
  truncee = NULL,
  truncee_rate = NULL,
  stripright = 0,
  stripleft = 0,
  maxee_rate = 0.01,
  minlen = 0,
  maxlen = NULL,
  maxns = 0,
  minsize = NULL,
  maxsize = NULL,
  minqual = 0,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fasta_width = 0,
  sample = NULL,
  stats = TRUE,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ object containing
  (forward) reads. See *Details*.

- reverse:

  (Optional). A FASTA/FASTQ file path or object containing reverse
  reads. If `fastx_input` is a `"pe_df"` object and `reverse` is not
  provided, the reverse reads will be extracted from its `"reverse"`
  attribute.

- output_format:

  (Optional). Desired output format of file or tibble: `"fasta"` or
  `"fastq"` (default). If `fastx_input` is a FASTA file path or a FASTA
  object, `output_format` cannot be `"fastq"`.

- fastaout:

  (Optional). Name of the FASTA output file for the sequences given in
  `fastx_input`. If `NULL` (default), no FASTA sequences are written to
  file. See *Details*.

- fastqout:

  (Optional). Name of the FASTQ output file for the sequences given in
  `fastx_input`. If `NULL` (default), no FASTQ sequences are written to
  file. See *Details*.

- fastaout_rev:

  (Optional). Name of the FASTA output file for the reverse sequences.
  If `NULL` (default), no FASTA sequences are written to file. See
  *Details*.

- fastqout_rev:

  (Optional). Name of the FASTQ output file for the reverse sequences.
  If `NULL` (default), no FASTQ sequences are written to file. See
  *Details*.

- trunclen:

  (Optional). Truncate sequences to the specified length. Shorter
  sequences are discarded. If `NULL` (default), the trimming is not
  applied.

- truncqual:

  (Optional). Truncate sequences starting from the first base with a
  quality score of the specified value or lower. Defaults to `1`.

- truncee:

  (Optional). Truncate sequences so that their total expected error does
  not exceed the specified value. If `NULL` (default), the trimming is
  not applied.

- truncee_rate:

  (Optional). Truncate sequences so that their average expected error
  per base is not higher than the specified value. The truncation will
  happen at first occurrence. The average expected error per base is
  calculated as the total expected number of errors divided by the
  length of the sequence after truncation. If `NULL` (default), the
  trimming is not applied.

- stripright:

  (Optional). Number of bases stripped from the right end of the reads.
  Defaults to `0`.

- stripleft:

  (Optional). Number of bases stripped from the left end of the reads.
  Defaults to `0`.

- maxee_rate:

  (Optional). Threshold for average expected error. Numeric value
  ranging form `0.0` to `1.0`. Defaults to `0.01`. See *Details*.

- minlen:

  (Optional). Minimum number of bases a sequence must have to be
  retained. Defaults to `0`. See *Details*.

- maxlen:

  (Optional). Maximum number of bases a sequences can have to be
  retained. If `NULL` (default), the filter is not applied.

- maxns:

  (Optional). Maximum number of N's for a given sequence. Sequences with
  more N's than the specified number are discarded. Defaults to `0`.

- minsize:

  (Optional). Minimum abundance for a given sequence. Sequences with
  lower abundance are discarded. If `NULL` (default), the filter is not
  applied.

- maxsize:

  (Optional). Maximum abundance for a given sequence. Sequences with
  higher abundance are discarded. If `NULL` (default), the filter is not
  applied.

- minqual:

  (Optional). Minimum base quality for a read to be retained. A read is
  discarded if it contains bases with a quality score below the given
  value. Defaults to `0`, meaning no reads are discarded.

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

- stats:

  (Optional). If `TRUE` (default), a tibble with statistics about the
  filtering is added as an attribute of the returned tibble. If `FALSE`,
  no statistics are added.

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

If output files are specified, the results are written directly to the
specified output files, and no tibble is returned.

If output files (`fastaout`/`fastqout` and
`fastaout_rev`/`fastqout_rev`) are `NULL`, a tibble containing the
trimmed and/or filtered reads from `fastx_input` in the format specified
by `output_format` is returned.

If `reverse` is provided, a tibble containing the trimmed and/or
filtered reverse sequences is attached as an attribute, named
`"reverse"` to the returned table.

When the reverse reads are present, the returned tibble is assigned the
class `"pe_df"`, identifying it as paired-end data.

The `"statistics"` attribute of the returned tibble (when output files
are `NULL`) is a tibble with the following columns:

- `Kept_Sequences`: Number of retained sequences.

- `Discarded_Sequences`: Number of discarded sequences.

- `fastx_source`: Name of the file/object with forward (R1) reads.

- `reverse_source`: (If `reverse` is specified) Name of the file/object
  with reverse (R2) reads.

## Details

Reads from the input files/objects (`fastx_input` and `reverse`) are
trimmed and/or filtered based on the specified criteria using `VSEARCH`.

`fastx_input` and `reverse` can either be file paths to FASTA/FASTQ
files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
columns `Header` and `Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

If `fastx_input` is an object of class `"pe_df"`, the reverse reads are
automatically extracted from its `"reverse"` attribute unless explicitly
provided via the `reverse` argument.

If `reverse` is provided, it is processed alongside `fastx_input` using
the same trimming/filtering criteria.

Note that if you want to trim/filter the forward and reverse reads
differently, you must pass them separately to this function, get two
result files/objects, and then use
[`fastx_synchronize`](https://cassandrahjo.github.io/Rsearch/reference/fastx_synchronize.md)
to synchronize the read pairs again.

If `fastaout` and `fastaout_rev` or `fastqout` and `fastqout_rev` are
specified, trimmed and/or filtered sequences are written to these files
in the specified format.

If output files are `NULL`, results are returned as a tibbles. When
returning tibbles, the reverse sequences (if provided) are attached as
an attribute named `"reverse"`.

When reverse reads are returned as an attribute, the primary tibble is
also assigned the S3 class `"pe_df"` to indicate that it represents
paired-end data. This class tag can be used by downstream tools to
recognize paired-end tibbles.

Note that certain options are not compatible with both file formats. For
instance, options that trim or filter sequences based on quality scores
are unavailable when the input is of type `"fasta"`. Visit the `VSEARCH`
[documentation](https://github.com/torognes/vsearch?tab=readme-ov-file#getting-help)
for more details.

Sequences with an average expected error greater than the specified
`maxee_rate` are discarded. For a given sequence, the average expected
error is the sum of error probabilities for all the positions in the
sequence, divided by the length of the sequence.

Any input sequence with fewer bases than the value set in `minlen` will
be discarded. By default, `minlen` is set to 0, which means that no
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
fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
                         "small_R1.fq")
reverse <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R1.fq")
output_format <- "fastq"
maxee_rate <- 0.01
minlen <- 0

# Trim/filter sequences and return a FASTQ tibble
filt_seqs <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                reverse = reverse,
                                output_format = output_format,
                                maxee_rate = maxee_rate,
                                minlen = minlen)

# Extract tibbles
R1_filt <- filt_seqs
R2_filt <- attr(filt_seqs, "reverse")

# Extract filtering statistics
statistics <- attr(filt_seqs, "statistics")

# Trim/filter sequences and write results to FASTQ files
vs_fastx_trim_filt(fastx_input = fastx_input,
                   reverse = reverse,
                   fastqout = "filt_R1.fq",
                   fastqout_rev = "filt_R2.fq",
                   output_format = output_format,
                   maxee_rate = maxee_rate,
                   minlen = minlen)
} # }
```
