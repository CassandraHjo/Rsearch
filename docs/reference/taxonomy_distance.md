# Creates a distance object based on taxonomy information.

Creates a distance matrix based on taxonomy information

## Usage

``` r
taxonomy_distance(taxonomy_table, confidence = NULL)
```

## Arguments

- taxonomy_table:

  (Required). A data.frame with taxonomy information, see *Details*.

- confidence:

  (Optional). A threshold value used to replace taxa with confidence
  scores below this to `NA`.

## Value

A [`dist`](https://rdrr.io/r/stats/dist.html) object containing
taxonomic distances between OTUs.

## Details

In some data analyses involving OTU data, it is often useful to quantify
the relatedness of OTUs based on their taxonomy. This function creates a
distance matrix from a taxonomy table of the same format as output by
[`vs_sintax`](https://cassandrahjo.github.io/Rsearch/reference/vs_sintax.md).
This means `taxonomy_table` must have columns domain, phylum, class,
order, family, genus and species. It must also have a column Header with
a text that unique to each row. These are used as row/column names in
the returned [`dist`](https://rdrr.io/r/stats/dist.html) object.

Distances between two OTUs reflect how high up in the taxonomy they have
a common taxon, i.e if they are distinct OTUs but of the same species
the distance is 1, if they are different species but same genus the
distance is 2 etc. Note that `NA`s in the taxonomy are not matched,
increasing the distances, i.e if two OTUs have `NA` as species and
genus, but share family, the distance is 3 (implicitly assuming they are
different genera but same family).

The `confidence` sets a threshold for replacing low-confidence taxa to
`NA`. For this to work the `taxonomy_table` must have columns with such
confidence scores i.e. columns domain_score, phylum_score,
...species_score. If the species_score is below `confidence` the
corresponding species name is set to `NA`, and similar for all ranks.
The default is to ignore this confidence (`confidence = NULL`).

## References

<https://www.biorxiv.org/content/10.1101/074161v1>

## Examples

``` r
if (FALSE) { # \dontrun{
# Assign taxonomy with sintax
db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "sintax_db.fasta")
fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "small.fasta")
tax.tbl <- vs_sintax(fasta_input = fasta.file, database = db.file)

# Calculate distance matrix
tax.dist <- taxonomy_distance(tax.tbl)

# You can now directly use 'tax.dist' with functions like hclust or ape::nj
tax.tree <- ape::nj(tax.dist)
} # }
```
