# Optimize read truncation with truncee_rate

`vs_optimize_truncee_rate` optimizes the truncation parameter
`truncee_rate` to achieve the best possible merging results. The
function iterates through a specified range of `truncee_rate` values to
identify the optimal value that maximizes the proportion of high-quality
merged read pairs.

## Usage

``` r
vs_optimize_truncee_rate(
  fastq_input,
  reverse = NULL,
  sample_size = 10000,
  minovlen = 10,
  truncee_rate_range = c(seq(0.002, 0.04, by = 0.002)),
  minlen = 1,
  min_size = 2,
  maxee_rate = 0.01,
  threads = 1,
  plot_title = TRUE,
  tmpdir = NULL,
  verbose = TRUE
)
```

## Arguments

- fastq_input:

  (Required). A FASTQ file path, FASTQ tibble (forward reads), or a
  paired-end tibble of class `"pe_df"`. See *Details*.

- reverse:

  (Optional). A FASTQ file path or FASTQ tibble (reverse reads).
  Optional if `fastq_input` is a `"pe_df"` object.

- sample_size:

  (Optional). Number of read pairs to randomly sample for the
  optimization process. Sampling speeds up the process on large
  datasets. Set to `NULL` to use all reads. Defaults to `10000`.

- minovlen:

  (Optional). Minimum overlap between the merged reads. Must be at
  least 5. Defaults to `10`.

- truncee_rate_range:

  (Optional). A numeric vector of `truncee_rate` values to test.
  Sequences are truncated so that their average expected error per base
  is lower than the specified value. Defaults to
  `(0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024, 0.026, 0.028, 0.030, 0.032, 0.034, 0.036, 0.038, 0.040)`.

- minlen:

  (Optional). Minimum number of bases a sequence must have to be
  retained. Defaults to `0`. See *Details*.

- min_size:

  (Optional). Minimum copy number (size) for a merged read to be
  included in the results. Defaults to `2`.

- maxee_rate:

  (Optional). Threshold for average expected error. Must range from
  `0.0` to `1.0`. Defaults to `0.01`. See *Details*.

- threads:

  (Optional). Number of computational threads to be used by `VSEARCH`.
  Defaults to `1`.

- plot_title:

  (Optional). If `TRUE` (default), a summary title will be displayed in
  the plot. Set to `FALSE` for no title.

- tmpdir:

  (Optional). Path to the directory where temporary files should be
  written when tables are used as input or output. Defaults to `NULL`,
  which resolves to the session-specific temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- verbose:

  (Optional). Logical. Controls verbosity of the function. If `TRUE`, a
  progress bar is displayed. Defaults to `TRUE`.

## Value

A data frame with the following columns:

- `truncee_rate_value`: Tested `truncee_rate` value.

- `merged_read_pairs`: Count of merged read-pairs with a copy number
  above `min_size` after dereplication.

- `R1_length`: Average length of R1-reads after trimming.

- `R2_length`: Average length of R2-reads after trimming.

The returned data frame has an attribute named `"plot"` containing a
[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object based on the returned data frame. The plot visualizes
`truncee_rate` values against `merged_read_pairs`, `R1_length`, and
`R2_length`, with the optimal `truncee_rate` value marked by a red
dashed line.

Additionally, the returned data frame has an attribute named
`"optimal_truncee_rate"` containing the optimal `truncee_rate` value.

## Details

The function uses
[`vs_fastq_mergepairs`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastq_mergepairs.md),
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md),
and
[`vs_fastx_uniques`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.md)
where the arguments to this functions are described in detail.

If `fastq_input` has class `"pe_df"`, the reverse reads will be
automatically extracted from the `"reverse"` attribute unless explicitly
provided in the `reverse` argument.

The best possible truncation option (`truncee_rate`) for merging is
measured by the number of merged read-pairs with a copy number above the
number specified by `min_size` after dereplication.

Changing `min_size` will affect the results. A low `min_size` will
include merged sequences with a lower copy number after dereplication,
and a higher `min_size` will filter out more reads and only count
high-frequency merged sequences.

## References

<https://github.com/torognes/vsearch>

## See also

[`vs_fastq_mergepairs`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastq_mergepairs.md),
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md),
[`vs_fastx_uniques`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R1.fq")
R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R2.fq")

# Run optimizing function
optimize.tbl <- vs_optimize_truncee_rate(fastq_input = R1.file,
                                         reverse = R2.file)

# Display plot
print(attr(optimize.tbl, "plot"))

} # }
```
