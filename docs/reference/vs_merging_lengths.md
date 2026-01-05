# Length statistics after merging

`vs_merging_lengths` computes length statistics for forward reads,
reverse reads, merged reads, and their overlaps before and after
merging.

## Usage

``` r
vs_merging_lengths(
  fastq_input,
  reverse = NULL,
  minovlen = 10,
  minlen = 0,
  threads = 1,
  plot_title = TRUE,
  tmpdir = NULL
)
```

## Arguments

- fastq_input:

  (Required). A FASTQ file path, a FASTQ tibble (forward reads), or a
  paired-end tibble of class `"pe_df"`. See *Details*.

- reverse:

  (Optional). A FASTQ file path or FASTQ tibble containing reverse
  reads. Optional if `fastq_input` is a `"pe_df"` object.

- minovlen:

  (Optional). Minimum overlap between the merged reads. Must be at
  least 5. Defaults to `10`.

- minlen:

  (Optional). Minimum number of bases a sequence must have to be
  retained. Defaults to `0`. See *Details*.

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

## Value

A tibble with the following columns:

- `length_1`: The length of the forward reads.

- `length_2`: The length of the reverse reads.

- `length_merged`: The length of the merged reads.

- `length_overlap`: The length of the overlap between the forward and
  reverse reads.

In case of missing values for the latter two columns, it means that the
corresponding reads were not merged.

The tibble includes additional attributes:

- `plot`:

  A
  [`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
  object visualizing the returned data frame.

- `statistics`:

  Additional statistics returned from
  [`vs_fastq_mergepairs`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastq_mergepairs.md).

## Details

The function uses
[`vs_fastq_mergepairs`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastq_mergepairs.md)
where the arguments to this function are described in detail.

If `fastq_input` is an object of class `"pe_df"`, the reverse reads are
automatically extracted from its `"reverse"` attribute unless explicitly
provided via the `reverse` argument. This allows streamlined input
handling for paired-end tibbles created by
[`fastx_synchronize`](https://cassandrahjo.github.io/Rsearch/reference/fastx_synchronize.md)
or
[`vs_fastx_trim_filt`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_trim_filt.md).

These length statistics are most typically used in order to tune the
filter and trimming of reads such that the merged reads are of high
quality.

## References

<https://github.com/torognes/vsearch>

## See also

[`vs_fastq_mergepairs`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastq_mergepairs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define arguments
R1.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R1.fq")
R2.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small_R2.fq")

# Run function
merging.tbl <- vs_merging_lengths(fastq_input = R1.file,
                                  reverse = R2.file)

# Display plot
merging_stats_plot <- attr(merging.tbl, "plot")
print(merging_stats_plot)

} # }
```
