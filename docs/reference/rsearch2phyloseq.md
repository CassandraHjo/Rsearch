# Convert Rsearch object to phyloseq object

`rsearch2phyloseq` converts an Rsearch object to a phyloseq object.

## Usage

``` r
rsearch2phyloseq(rsearch.obj, sample_id_col = "sample_id")
```

## Arguments

- rsearch.obj:

  (Required). An Rsearch object, see
  [`rsearch_obj`](https://cassandrahjo.github.io/Rsearch/reference/rsearch_obj.md).

- sample_id_col:

  (Optional). A character string specifying the name of the column in
  `sampledata.df` that contains sample identifiers. Defaults to
  `"sample_id"`.

## Value

A [`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html) object.

## Details

This function converts an Rsearch object, which is a simple `list`, to a
[`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html) object from
the `phyloseq` R package.

## References

<https://joey711.github.io/phyloseq/>

## See also

[`rsearch_obj`](https://cassandrahjo.github.io/Rsearch/reference/rsearch_obj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Rsearch object to phyloseq object
phy_obj <- rsearch2phyloseq(obj, sample_id_col = "sample_id")
} # }
```
