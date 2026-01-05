# Make Sintax database

Creates a properly formatted FASTA file for the use as a Sintax
database.

## Usage

``` r
make_sintax_db(taxonomy_table, outfile)
```

## Arguments

- taxonomy_table:

  (Required). A data.frame with sequences and proper information for
  making a Sintax database, see *Details*.

- outfile:

  (Required). Name of database file to create (a FASTA file).

## Value

No return in R, but a FASTA file (`outfile`) with properly formatted
`Header` lines is created.

## Details

The Sintax algorithm is used by `VSEARCH` to assign taxonomic
information to 16S sequences. It requires a database, which is nothing
but a FASTA file of 16S sequences with properly formatted
`Header`-lines.

The `taxonomy_table` provided as input here must have the columns:

- `Header` - short unique text for each sequence

- `Sequence` - the sequences

- Columns `domain`, `phylum`, `class`, `order`, `family`, `genus`,
  `species`. Text columns with taxon names.

In some taxonomies the domain rank is named kingdom, but here we use the
word domain. You may very well have empty (NA) entries in the taxonomy
columns of the table.

## References

<https://www.biorxiv.org/content/10.1101/074161v1>

## Examples

``` r
if (FALSE) { # \dontrun{
# First, you need a table of the same format as output by vs_sintax:
db.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                     "sintax_db.fasta")
fasta.file <- file.path(file.path(path.package("Rsearch"), "extdata"),
                        "small.fasta")
tax.tbl <- vs_sintax(fasta_input = fasta.file, database = db.file)

# Inspect tax.tbl to see its columns. You replace the column content with
# your desired taxonomy.
# From such a tax.tbl you create the database file:
make_sintax_db(tax.tbl, outfile = "delete_ma.fasta")
} # }
```
