test_that("vs_sintax fails with invalid strand", {
  fasta_file <- test_path("testdata", "small.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")

  expect_error(
    vs_sintax(fasta_input = fasta_file,
              database = db_file,
              strand = "invalid"),
    "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("vs_sintax errors if input fasta file does not exist", {
  fake_fasta <- tempfile(fileext = ".fa")
  db_file <- test_path("testdata", "sintax_db.fasta")

  expect_error(
    vs_sintax(fasta_input = fake_fasta, database = db_file),
    "Cannot find .*\\.fa")
})

test_that("vs_sintax errors if input database file does not exist", {
  fasta_input <- tibble::tibble(
    Header = "seq1",
    Sequence = "ATCGATCG"
  )
  fake_db <- tempfile(fileext = ".fa")

  expect_error(
    vs_sintax(fasta_input = fasta_input, database = fake_db),
    paste("Cannot find input file:", normalizePath(fake_db, mustWork = FALSE)))
})

test_that("vs_sintax returns expected columns in classification", {
  fasta_file <- test_path("testdata", "R1.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")

  result <- vs_sintax(fasta_input = fasta_file, database = db_file)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("Header", "domain", "phylum", "class", "order",
                    "family", "genus", "species") %in% names(result)))
})

test_that("vs_sintax writes output to file if outfile is specified", {
  fasta_file <- test_path("testdata", "R1.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")
  outfile <- withr::local_tempfile(fileext = ".txt")

  result <- vs_sintax(fasta_input = fasta_file,
                      database = db_file,
                      outfile = outfile)

  expect_null(result)
  expect_true(file.exists(outfile))
  out_tbl <- read.table(outfile, header = TRUE)
  expect_true("Header" %in% names(out_tbl))
})

test_that("vs_sintax creates logfile if specified", {
  fasta_file <- test_path("testdata", "R1.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")
  log_file <- withr::local_tempfile()

  vs_sintax(fasta_input = fasta_file, database = db_file, logfile = log_file)

  expect_true(file.exists(log_file))
})

test_that("vs_sintax accepts additional vsearch_options", {
  fasta_file <- test_path("testdata", "R1.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")

  result <- vs_sintax(fasta_input = fasta_file,
                      database = db_file,
                      vsearch_options = c(""))

  expect_s3_class(result, "data.frame")
})

test_that("vs_sintax errors when database tibble lacks required columns", {
  fasta_input <- tibble::tibble(
    Header = "seq1",
    Sequence = "ATCGATCG"
  )

  fake_db <- tibble::tibble(ID = "A", Seq = "ATCG")

  expect_error(
    vs_sintax(fasta_input = fasta_input, database = fake_db),
    "FASTA data base must contain columns: Header and Sequence"
  )
})

test_that("make_sintax_db produces a database", {

  input_df <- readRDS(test_path("testdata", "make_sintax_db.rds"))

  temp_out_db <- withr::local_tempfile()

  expect_invisible(make_sintax_db(taxonomy_table = input_df,
                                  outfile = temp_out_db))

  expect_true(file.exists(temp_out_db))
  expect_equal(
    microseq::readFasta(
      test_path("testdata", "output", "sintax_db_synthetic.fasta")),
    microseq::readFasta(temp_out_db))
})

test_that("make_sintax_db errors when required columns are missing", {
  input_df <- readRDS(test_path("testdata", "make_sintax_db.rds"))

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-Header), outfile),
    "must have a column named Header, with a unique text for each sequence"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-Sequence), outfile),
    "must have a column named Sequence, with the sequences"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-domain), outfile),
    "must have a column named domain"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-phylum), outfile),
    "must have a column named phylum"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-class), outfile),
    "must have a column named class"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-order), outfile),
    "must have a column named order"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-family), outfile),
    "must have a column named family"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-genus), outfile),
    "must have a column named genus"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-species), outfile),
    "must have a column named species"
  )
})

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

test_that("vs_sintax returns empty data.frame when no classification match", {
  fasta_input <- tibble::tibble(Header = "unknown", Sequence = "NNNNNNNNNN")
  db_file <- test_path("testdata", "sintax_db.fasta")

  result <- vs_sintax(fasta_input = fasta_input, database = db_file)

  expect_true(nrow(result) == 1)
  expect_true(all(is.na(result$domain)))
})

test_that("vs_sintax works with database as tibble", {
  fasta_input <- tibble::tibble(Header = "seq1", Sequence = "ACGTACGT")
  database <- tibble::tibble(
    Header = "ref1;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae,g:Lactobacillus,s:casei;",
    Sequence = "ACGTACGTACGT"
  )

  result <- vs_sintax(fasta_input = fasta_input, database = database)

  expect_s3_class(result, "data.frame")
  expect_true("domain" %in% names(result))
})

test_that("vs_sintax includes randseed argument", {
  fasta_file <- test_path("testdata", "R1.fasta")
  db_file <- test_path("testdata", "sintax_db.fasta")

  result <- vs_sintax(fasta_input = fasta_file, database = db_file, randseed = 123)

  expect_s3_class(result, "data.frame")
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


