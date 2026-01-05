# Synchronize FASTA and FASTQ files or objects

`fastx_synchronize` synchronizes sequences between two FASTA/FASTQ files
or objects by retaining only the common sequences present in both.

## Usage

``` r
fastx_synchronize(
  file1,
  file2 = NULL,
  file_format = "fastq",
  file1_out = NULL,
  file2_out = NULL
)
```

## Arguments

- file1:

  (Required). A FASTQ file path, a FASTQ tibble, or a paired-end tibble
  of class `"pe_df"`. See *Details*.

- file2:

  (Optional). A FASTQ file path or a FASTQ tibble. Optional if `file1`
  is a `"pe_df"` object. See *Details*.

- file_format:

  (Optional). Format of the input (`file1` and `file2`) and the desired
  output format: `"fasta"` or `"fastq"` (default). This determines the
  format for both outputs.

- file1_out:

  (Optional). Name of the output file for synchronized reads from
  `file1`. The file is in either FASTA or FASTQ format, depending on
  `file_format`. If `NULL` (default), no sequences are written to a
  file. See *Details*.

- file2_out:

  (Optional). Name of the output file for synchronized reads from
  `file2`. The file is in either FASTA or FASTQ format, depending on
  `file_format`. If `NULL` (default), no sequences are written to a
  file. See *Details*.

## Value

A tibble or `NULL`.

If both `file1_out` and `file2_out` are `NULL`, a tibble containing the
synchronized reads from `file1` is returned. The synchronized reads from
`file2` are accessible via the `"reverse"` attribute of the returned
tibble.

If both `file1_out` and `file2_out` are specified, the synchronized
sequences are written to the specified output files, and no tibble is
returned.

## Details

`file1` and `file2` can either be paths to FASTA/FASTQ files or tibble
objects containing the sequences. FASTA objects are tibbles that contain
the columns `Header` and `Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

If `file1` is an object of class `"pe_df"`, the second read tibble is
automatically extracted from its `"reverse"` attribute unless explicitly
provided via the `file2` argument. This allows streamlined input
handling for paired-end tibbles created by
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md).

Sequence IDs in the `Header` fields must be identical for each read pair
in both `file1` and `file2` for synchronization to work correctly.

If `file1_out` and `file2_out` are specified, the synchronized sequences
are written to these files in the format specified by `file_format`.

If `file1_out` and `file2_out` are `NULL`, the function returns a
FASTA/FASTQ object containing synchronized reads from `file1`. The
synchronized reads from `file2` are included as an attribute named
`"reverse"` in the returned tibble.

The returned tibble is assigned the S3 class `"pe_df"`, indicating that
it represents paired-end sequence data. Downstream functions can use
this class tag to distinguish paired-end tibbles from other tibbles.

Both `file1_out` and `file2_out` must either be `NULL` or both must be
character strings specifying the file paths.

## Examples

``` r
# Define arguments
file1 <- system.file("extdata/small_R1.fq", package = "Rsearch")
file2 <- system.file("extdata/small_R1.fq", package = "Rsearch")
file_format <- "fastq"
file1_out <- NULL
file2_out <- NULL

# Synchronize files and return as a tibble
sync_seqs <- fastx_synchronize(file1 = file1,
                               file2 = file2,
                               file_format = file_format,
                               file1_out = file1_out,
                               file2_out = file2_out)

# Extract tibbles with synchronized sequences
R1_sync <- sync_seqs
R2_sync <- attr(sync_seqs, "reverse")

# Synchronize files and write to output files

# Define output file names
out1 <- tempfile(fileext = ".fastq")
out2 <- tempfile(fileext = ".fastq")

fastx_synchronize(file1 = file1,
                  file2 = file2,
                  file_format = file_format,
                  file1_out = out1,
                  file2_out = out2)

```
