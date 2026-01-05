# Merge paired-end sequence reads

`vs_fastq_mergepairs` merges paired-end sequence reads with overlapping
regions into one sequence using `VSEARCH`.

## Usage

``` r
vs_fastq_mergepairs(
  fastq_input,
  reverse = NULL,
  output_format = "fasta",
  fastaout = NULL,
  fastqout = NULL,
  minovlen = 10,
  minlen = 0,
  fasta_width = 0,
  sample = NULL,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastq_input:

  (Required). A FASTQ file path, a FASTQ tibble (forward reads), or a
  paired-end tibble of class `"pe_df"`. See *Details*.

- reverse:

  (Optional). A FASTQ file path or a FASTQ tibble (reverse reads).
  Optional if `fastq_input` is a `"pe_df"` object. See *Details*.

- output_format:

  (Optional). Desired output format of file or tibble: `"fasta"`
  (default) or `"fastq"`.

- fastaout:

  (Optional). Name of the FASTA output file with the merged reads. If
  `NULL` (default), no output is written to file. See *Details*.

- fastqout:

  (Optional). Name of the FASTQ output file with the merged reads. If
  `NULL` (default) no output is written to file. See *Details*.

- minovlen:

  (Optional). Minimum overlap between the merged reads. Must be at
  least 5. Defaults to `10`.

- minlen:

  (Optional). Minimum number of bases a sequence must have to be
  retained. Defaults to `0`. See *Details*.

- fasta_width:

  (Optional). Number of characters per line in the output FASTA file.
  Only applies if the output file is in FASTA format. Defaults to `0`,
  which eliminates wrapping.

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

If `fastaout` or `fastqout` is specified , the merged sequences are
written to the specified output file, and no tibble is returned.

If `fastaout` or `fastqout` is `NULL`, a tibble containing the merged
reads in the format specified by `output_format` is returned.

The `"statistics"` attribute of the returned tibble (when `fastaout` or
`fastqout` is `NULL`) is a tibble with the following columns:

- `Tot_num_pairs`: Total number of read pairs before merging.

- `Merged`: Number of read pairs that merged.

- `Mean_Read_Length_before_merging`: Mean read length before merging (R1
  and R2).

- `Mean_Read_Length_after_merging`: Mean read length after merging.

- `StdDev_Read_Length`: Standard deviation of read length after merging.

- `R1`: Name of the file/object with forward (R1) reads used in the
  merging.

- `R2`: Name of the file/object with reverse (R2) reads used in the
  merging.

## Details

Read pairs from the input FASTQ files (`fastq_input` and `reverse`) are
merged into a single sequence by overlapping regions. The resulting
sequences consist of the merged forward and reverse reads with the
specified minimum overlap.

`fastq_input` and `reverse` can either be file paths to FASTQ files or
FASTQ objects. FASTQ objects are tibbles that contain the columns
`Header`, `Sequence`, and `Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html). Forward
and reverse reads must appear in the same order and have the same total
number of reads in both files.

If `fastq_input` is an object of class `"pe_df"`, the reverse reads are
automatically extracted from its `"reverse"` attribute unless explicitly
provided via the `reverse` argument. This allows streamlined input
handling for paired-end tibbles created by
[`fastx_synchronize`](https://cassandrahjo.github.io/Rsearch/reference/fastx_synchronize.md)
or
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md).

If `fastaout` or `fastqout` is specified, the merged reads are written
to the respective file in either FASTA or FASTQ format.

If both `fastaout` or `fastqout` are `NULL`, the results are returned as
a FASTA or FASTQ object, and no file is written.

`output_format` has to match the desired output files/objects.

Any input sequence with fewer bases than the value set in `minlen` will
be discarded. Default `minlen` is 0, meaning no sequences are removed.
However, using the default value may allow empty sequences to remain in
the results.

If `log_file` is `NULL` and `fastqout` or `fastaout` is specified,
merging statistics from `VSEARCH` will not be captured.

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

# Merge sequences and return a FASTQ tibble
merge_seqs <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                  reverse = reverse,
                                  output_format = output_format)

# Extract merging statistics
statistics <- attr(merge_seqs, "statistics")

# Merge sequences and write sequences to a FASTQ file
vs_fastq_mergepairs(fastq_input = fastq_input,
                    reverse = reverse,
                    output_format = output_format,
                    fastqout = "merged_sequences.fq")
} # }
```
