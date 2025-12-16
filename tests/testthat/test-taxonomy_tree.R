test_that("taxonomy_tree errors when required columns are missing", {
  df <- tibble::tibble(Header = "seq1", domain = "Bacteria")

  expect_error(taxonomy_tree(df), "must have a column named phylum")
})

test_that("taxonomy_tree returns a phylo object", {
  taxonomy_table <- readRDS(test_path("testdata", "tax_tbl.rds"))

  tree <- taxonomy_tree(taxonomy_table)

  expect_s3_class(tree, "phylo")
})

test_that("taxonomy_tree filters low-confidence taxa", {
  tbl <- readRDS(test_path("testdata", "tax_tbl.rds"))

  filtered <- taxonomy_tree(tbl, confidence = 0.95)

  expect_s3_class(filtered, "phylo")
})

test_that("taxonomy_tree errors when required columns are missing", {
  base_tbl <- tibble::tibble(
    Header = "seq1",
    domain = "Bacteria",
    phylum = "Firmicutes",
    class = "Bacilli",
    order = "Lactobacillales",
    family = "Lactobacillaceae",
    genus = "Lactobacillus",
    species = "casei"
  )

  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -Header)),
    "must have a column named Header"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -domain)),
    "must have a column named domain"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -phylum)),
    "must have a column named phylum"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -class)),
    "must have a column named class"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -order)),
    "must have a column named order"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -family)),
    "must have a column named family"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -genus)),
    "must have a column named genus"
  )
  expect_error(
    taxonomy_tree(dplyr::select(base_tbl, -species)),
    "must have a column named species"
  )
})
