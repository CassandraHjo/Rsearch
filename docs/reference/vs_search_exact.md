# Search for exact full-length matches

`vs_search_exact` searches for exact full-length matches of query
sequences in a database of target sequences using `VSEARCH`.

## Usage

``` r
vs_search_exact(
  fastx_input,
  database,
  userout = NULL,
  otutabout = NULL,
  userfields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
  strand = "plus",
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ tibble object
  containing the query sequences. See *Details*.

- database:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ tibble object
  containing the target sequences.

- userout:

  (Optional). A character string specifying the name of the output file
  for the alignment results. If `NULL` (default), no output is written
  to a file and the results are returned as a tibble with the columns
  specified in `userfields`. See *Details*.

- otutabout:

  (Optional). A character string specifying the name of the output file
  in an OTU table format. If `NULL` (default), no output is written to a
  file. If `TRUE`, the output is returned as a tibble. See *Details*.

- userfields:

  (Optional). Fields to include in the output file. Defaults to
  `"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"`. See
  *Details*.

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

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

A tibble or `NULL`.

If `userout` is `NULL` a tibble containing the alignment results with
the fields specified by `userfields` is returned. If `userout` is
specified the alignment results are written to the specified file, and
no tibble is returned.

If `otutabout` is `TRUE`, an OTU table is returned as a tibble. If
`otutabout` is a character string, the output is written to the file,
and no tibble is returned.

If neither `userout` nor `otutabout` is specified, a tibble containing
the alignment results is returned.

## Details

Identifies exact full-length matches between query and target sequences
using `VSEARCH`. Only 100 specificity and making this command much
faster than
[`vs_usearch_global`](https://cassandrahjo.github.io/Rsearch/reference/vs_usearch_global.md).

`fastx_input` and `database` can either be file paths to a FASTA/FASTQ
files or FASTA/FASTQ objects. FASTA objects are tibbles that contain the
columns `Header` and `Sequence`, see
[`readFasta`](https://rdrr.io/pkg/microseq/man/readFasta.html). FASTQ
objects are tibbles that contain the columns `Header`, `Sequence`, and
`Quality`, see
[`readFastq`](https://rdrr.io/pkg/microseq/man/readFastq.html).

`userfields` specifies the fields to include in the output file. Fields
must be given as a character string separated by `"+"`. The default
value of `userfields` equals
`"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits"`, which
gives a blast-like tab-separated format of twelve fields. See the
'Userfields' section in the `VSEARCH` manual for more information.

`otutabout` gives the option to output the results in an OTU table
format with tab-separated columns. When writing to a file, the first
line starts with the string "#OTU ID", followed by a tab-separated list
of all sample identifiers (formatted as "sample=X"). Each subsequent
line, corresponding to an OTU, begins with the OTU identifier and is
followed by tab-separated abundances for that OTU in each sample. If
`otutabout` is a character string, the output is written to the
specified file. If `otutabout` is `TRUE`, the function returns the OTU
table as a tibble, where the first column is named `otu_id` instead of
"#OTU ID".

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

## References

<https://github.com/torognes/vsearch>

## See also

[`vs_usearch_global`](https://cassandrahjo.github.io/Rsearch/reference/vs_usearch_global.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# You would typically use something else as database
query_file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                        "small.fasta")
db <- query_file

# Search for exact full-length matches with default parameters, with file as output
vs_search_exact(fastx_input = query_file,
                database = db,
                userout = "delete_me.txt")

# Read results, and give column names
result.tbl <- read.table("delete_me.txt",
                         sep = "\t",
                         header = FALSE,
                         col.names = c("query", "target", "id", "alnlen",
                                       "mism", "opens", "qlo", "qhi",
                                       "tlo", "thi", "evalue", "bits"))
} # }
```
