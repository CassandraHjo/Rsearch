test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_merging_lengths(fastq_input = R1,
                                   reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_merging_lengths(fastq_input = R1,
                                   reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("get merging lengths from merging two fastq files", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                           reverse = reverse)

  expect_equal(merging_lengths_df,
               readRDS(test_path("testdata", "output", "merging_lengths_sample1_fastq_files.rds")))

})

test_that("get merging lengths from merging two fastq tibbles", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                           reverse = reverse)

  expect_equal(merging_lengths_df,
               readRDS(test_path("testdata", "output", "merging_lengths_sample1_fastq_tibbles.rds")))

})
