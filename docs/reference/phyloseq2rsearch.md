# Convert phyloseq object to Rsearch object

Creating an Rsearch object (list) from a phyloseq object.

## Usage

``` r
phyloseq2rsearch(phyloseq.obj)
```

## Arguments

- phyloseq.obj:

  (Required). A phyloseq object, see
  [`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html).

## Value

A `list` with entries as in a Rsearch object, except that the
`sequence.tbl` do not contain sequences, only taxonomy.

## Details

This function converts a phyloseq object to a simple
[`list`](https://rdrr.io/r/base/list.html) with three elements as
dataframes (or tibbles). The entries are named according to the
structure used in
[`rsearch_obj`](https://cassandrahjo.github.io/Rsearch/reference/rsearch_obj.md)

## References

<https://joey711.github.io/phyloseq/>

## See also

[`rsearch_obj`](https://cassandrahjo.github.io/Rsearch/reference/rsearch_obj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert phyloseq object to Rsearch object
rsearch_obj <- phyloseq2rsearch(phy_obj)

# Extract read count data
rsearch_obj$readcount.mat

# Extract sample data
rsearch_obj$sampledata.df

# Extract sequence data
rsearch_obj$sequence.df
} # }
```
