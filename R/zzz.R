.onLoad <- function(libname, pkgname){
  load(system.file("extdata/vsearch_executable.rds", package = "Rsearch"))
  options(Rsearch.vsearch_executable = vsearch_executable)
}

# Declare global variables to avoid R CMD check NOTE
utils::globalVariables(c(
  "vsearch_executable", "Header", "centroid_size", "tag", "domain", "phylum",
  "family", "genus", "species", "Lower", "Upper", "MedianQuality", "MeanQuality",
  "EE_rate", "Sequence", "Length", "Quality", "size", "num_reads",
  "domain_score", "phylum_score", "class_score", "order_score", "family_score",
  "genus_score", "species_score", "tax.tbl", "length_1", "length_2",
  "length_merged", "length_overlap", "merged_read_pairs", "R1_length",
  "R2_length", "truncee_rate_value", "value", "metric", "truncqual_value", "db",
  "plus", "taxonomy", ".data", "type", "centroid", "member", "members",
  "#OTU ID", "otu_id", "cluster_size"
))

