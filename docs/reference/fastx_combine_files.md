# Combine FASTA/FASTQ files in a directory into a single file or object

`fastx_combine_files` combines all FASTA or FASTQ files within a
specified directory into a single output file or a tibble object.

## Usage

``` r
fastx_combine_files(
  files_dir,
  output_file = NULL,
  file_ext = ".fq",
  file_format = "fastq",
  tmpdir = NULL
)
```

## Arguments

- files_dir:

  (Required). A character string specifying the path to the directory
  containing the files to be combined. Files must be uncompressed.

- output_file:

  (Optional). A character string specifying the name of the output file.
  If `NULL` (default), the combined data is returned as a FASTA/FASTQ
  object depending on `file_format` instead of being written to a file.

- file_ext:

  (Optional). File extension of the files to be combined. Defaults to
  `".fq"`.

- file_format:

  (Optional). Format of files to be combined and the desired output
  format: either `"fasta"` or `"fastq"` (default). See *Details*.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

## Value

A tibble or `NULL`.

If `output_file` is specified, the combined sequences are written to the
specified file.

If `output_file` is `NULL`, the combined sequences are returned as a
tibble in the format specified by `file_format`.

## Details

`files_dir` must contain uncompressed FASTA or FASTQ files matching the
specified `file_ext`.

All files with the specified `file_ext` in `files_dir` are concatenated
into a single output file or tibble.

A FASTA object is a tibble containing the columns `Header` and
`Sequence`. A FASTQ object is a tibble containing the columns `Header`,
`Sequence`, and `Quality`.

If `output_file` is specified, the combined sequences are written to
this file in the format specified by `file_format`.

If `output_file` is `NULL`, the combined sequences are returned as a
tibble in the format specified by `file_format`, and no file is written.

## Examples

``` r
# Define arguments
files_dir <- system.file("extdata", package = "Rsearch")
output_file <- NULL
file_ext <- ".fq"
file_format <- "fastq"

# Combine files and return tibble object
combined_files <- fastx_combine_files(files_dir = files_dir,
                                      output_file = output_file,
                                      file_ext = file_ext,
                                      file_format = file_format)

# Combine files and write to output file

# Define output file name
out <- tempfile(fileext = ".fastq")

fastx_combine_files(files_dir = files_dir,
                    output_file = out,
                    file_ext = file_ext,
                    file_format = file_format)
```
