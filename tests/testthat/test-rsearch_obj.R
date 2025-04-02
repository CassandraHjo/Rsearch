test_that("rsearch_obj constructs object correctly from data frames", {

  # Simulated inputs
  readcount_data <- tibble::tibble(
    OTU = c("OTU1", "OTU2", "OTU3"),
    Sample1 = c(10, 0, 5),
    Sample2 = c(0, 8, 2)
  )

  sequence_data <- tibble::tibble(
    Header = c("OTU1", "OTU2", "OTU3"),
    Sequence = c("ATCG", "ATCC", "ATGG"),
    Kingdom = c("k:Bacteria", "k:Bacteria", "k:Bacteria")
  )

  sample_data <- tibble::tibble(
    sample_id = c("Sample1", "Sample2"),
    Environment = c("soil", "water")
  )

  # Create object
  obj <- rsearch_obj(readcount_data, sequence_data, sample_data)

  expect_type(obj, "list")
  expect_named(obj, c("readcount.mat", "sequence.df", "sampledata.df"))
  expect_s3_class(obj$sequence.df, "data.frame")
  expect_s3_class(obj$sampledata.df, "data.frame")
  expect_true(is.matrix(obj$readcount.mat))

  # Verify sample/OTU alignment
  expect_equal(rownames(obj$readcount.mat), obj$sequence.df$Header)
  expect_equal(colnames(obj$readcount.mat), obj$sampledata.df$sample_id)
})

test_that("rsearch_obj handles subsetting to intersecting samples and OTUs", {

  readcount_data <- tibble::tibble(
    OTU = c("OTU1", "OTU2"),
    Sample1 = c(1, 2),
    Sample3 = c(3, 4) # Sample3 is not in sample_data
  )

  sequence_data <- tibble::tibble(
    Header = c("OTU1", "OTU3"), # OTU3 is not in readcount_data
    Sequence = c("AAAA", "CCCC")
  )

  sample_data <- tibble::tibble(
    sample_id = c("Sample1", "Sample2"), # Sample2 is not in readcount_data
    Treatment = c("A", "B")
  )

  obj <- rsearch_obj(readcount_data, sequence_data, sample_data)

  expect_equal(dim(obj$readcount.mat), c(1, 1))  # Only OTU1 and Sample1 should remain
  expect_equal(obj$sequence.df$Header, "OTU1")
  expect_equal(obj$sampledata.df$sample_id, "Sample1")
})

test_that("rsearch_obj accepts input files", {

  # Create temp input files
  readcount_file <- withr::local_tempfile(fileext = ".tsv")
  sequence_file <- withr::local_tempfile(fileext = ".tsv")
  sample_file <- withr::local_tempfile(fileext = ".tsv")

  readr::write_tsv(
    tibble::tibble(OTU = c("OTU1"), Sample1 = c(10)),
    readcount_file
  )

  readr::write_tsv(
    tibble::tibble(Header = "OTU1", Sequence = "ATCG"),
    sequence_file
  )

  readr::write_tsv(
    tibble::tibble(sample_id = "Sample1", pH = 7.1),
    sample_file
  )

  obj <- rsearch_obj(readcount_file, sequence_file, sample_file)

  expect_type(obj, "list")
  expect_equal(dim(obj$readcount.mat), c(1, 1))
  expect_equal(obj$sequence.df$Header, "OTU1")
  expect_equal(obj$sampledata.df$sample_id, "Sample1")
})

test_that("rsearch_obj strips size annotation from sequence headers", {

  readcount_data <- tibble::tibble(OTU = "OTU1", Sample1 = 42)
  sequence_data <- tibble::tibble(Header = "OTU1;size=999", Sequence = "AGTC")
  sample_data <- tibble::tibble(sample_id = "Sample1")

  obj <- rsearch_obj(readcount_data, sequence_data, sample_data)

  expect_equal(obj$sequence.df$Header, "OTU1")
})
