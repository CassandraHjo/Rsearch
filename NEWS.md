# Rsearch 1.1.0

## New features

* Added `taxonomy_distance()` for calculating pairwise taxonomic
  distances between OTUs from a taxonomy table.

* Added support for filtering taxonomy assignments by confidence scores
  through a new `confidence` argument in `taxonomy_tree()` and
  `taxonomy_distance()`.

* Extended `vs_cluster_unoise()` and `vs_cluster_size()` with new arguments for 
controlling output from clustering output: `clusters` and `uc`. The `uc` 
argument was also added to `vs_usearch_global()`.

* Added a verbose argument to functions that print progress messages or other 
information to the console.

* Updated `vs_optimize_truncee_rate()` and `vs_optimize_truncqual()` with an 
option to randomly sample 10,000 sequences when the input contains more than 
10,000 sequences, thereby reducing computation time.

## Documentation

* Updated function documentation and package website to describe the new 
taxonomy-distance functionality, clustering output options, and other newly 
introduced arguments.

## Tutorial

* Updated the tutorial to display the correct figure
