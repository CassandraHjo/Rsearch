# Global pairwise alignment

`vs_usearch_global` performs global pairwise alignment of query
sequences against target sequences using `VSEARCH`.

## Usage

``` r
vs_usearch_global(
  fastx_input,
  database,
  userout = NULL,
  otutabout = NULL,
  userfields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits",
  gapopen = "20I/2E",
  gapext = "2I/1E",
  id = 0.7,
  strand = "plus",
  maxaccepts = 1,
  maxrejects = 32,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (Required). A FASTA/FASTQ file path or FASTA/FASTQ object. See
  *Details*.

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

- gapopen:

  (Optional). Penalties for gap opening. Defaults to `"20I/2E"`. See
  *Details*.

- gapext:

  (Optional). Penalties for gap extension. Defaults to `"2I/1E"`. See
  *Details*.

- id:

  (Optional). Pairwise identity threshold. Defines the minimum identity
  required for matches. Defaults to `0.7`.

- strand:

  (Optional). Specifies which strand to consider when comparing
  sequences. Can be either `"plus"` (default) or `"both"`.

- maxaccepts:

  (Optional). Maximum number of matching target sequences to accept
  before stopping the search for a given query. Defaults to `1`.

- maxrejects:

  (Optional). Maximum number of non-matching target sequences to
  consider before stopping the search for a given query. Defaults to 32.
  If `maxaccepts` and `maxrejects` are both set to 0, the complete
  database is searched.

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

If `userout` is specified the alignment results are written to the
specified file, and no tibble is returned. If `userout` is `NULL` a
tibble containing the alignment results with the fields specified by
`userfields` is returned.

If `otutabout` is `TRUE`, an OTU table is returned as a tibble. If
`otutabout` is a character string, the output is written to the file,
and no tibble is returned.

## Details

Performs global pairwise alignment between query and target sequences
using `VSEARCH`, and reports matches based on the specified pairwise
identity threshold (`id`). Only alignments that meet or exceed the
identity threshold are included in the output.

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

Pairwise identity (`id`) is calculated as the number of matching columns
divided by the alignment length minus terminal gaps.

`vsearch_options` allows users to pass additional command-line arguments
to `VSEARCH` that are not directly supported by this function. Refer to
the `VSEARCH` manual for more details.

Visit the `VSEARCH`
[documentation](https://github.com/torognes/vsearch?tab=readme-ov-file#getting-help)
for information about defining `gapopen` and `gapext`.

## References

<https://github.com/torognes/vsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# You would typically use something else as database
query_file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small.fasta")
db <- query_file

# Run global pairwise alignment with default parameters and write results to file
vs_usearch_global(fastx_input = query_file,
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
