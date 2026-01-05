# Make a taxonomy tree

Creates a phylo object based on taxonomy

## Usage

``` r
taxonomy_tree(taxonomy_table, confidence = NULL)
```

## Arguments

- taxonomy_table:

  (Required). A data.frame with taxonomy information, see *Details*.

- confidence:

  (Optional). A threshold value used to replace taxa with confidence
  scores below this to `NA`.

## Value

A phylo object, see [`nj`](https://rdrr.io/pkg/ape/man/nj.html).

## Details

In some data analyses involving OTU data a phylogenetic tree describing
the relatedness of the OTUs is required. To construct such trees you
typically need to make a multiple alignment of the sequences behind each
OTU, which is a huge job.

An alternative is then to simply use the taxonomy, and create a
'taxonomy-tree' instead of a phylogenetic tree. This function creates
such a tree from a taxonomy table of the same format as output by
[`vs_sintax`](https://cassandrahjo.github.io/Rsearch/reference/vs_sintax.md).

Distances between two OTUs reflect how high up in the taxonomy they have
a common taxon, i.e if they are of the same species the distance is 0,
if different species but same genus the distance is 1 etc. Note that
`NA`s in the taxonomy are not matched, increasing the distances, i.e if
two OTUs have `NA` as species and genus, but share family, the distance
is 2.

The `confidence` sets a threshold for replacing low-confidence taxa to
`NA`. For this to work the `taxonomy_table` must have columns with such
confidence scores i.e. columns domain_score, phylum_score,
...species_score. If the species_score is below `confidence` the
corresponding species name is set to `NA`, and similar for all ranks.
The default is to ignore this confidence (`confidence = NULL`).

From these distances a Neighbor Joining tree is built using
[`nj`](https://rdrr.io/pkg/ape/man/nj.html).

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

# Making tree
tax.tree <- taxonomy_tree(tax.tbl)
} # }
```
