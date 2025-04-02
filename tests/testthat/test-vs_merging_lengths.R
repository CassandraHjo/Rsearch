test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_merging_lengths(fastq_input = R1,
                                  reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) |>
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

  expected_df <- readRDS(test_path("testdata", "output", "merging_lengths_sample1_fastq_files.rds"))

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(merging_lengths_df, expected_df)
})

test_that("get merging lengths from merging two fastq tibbles", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                           reverse = reverse,
                                           plot_title = FALSE)

  expected_df <- readRDS(test_path("testdata", "output", "merging_lengths_sample1_fastq_tibbles.rds"))

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(merging_lengths_df, expected_df)

})

test_that("vs_merging_lengths creates histograms when R1 and R2 lengths vary", {
  # Load sample data
  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))

  # Artificially shorten half the sequences in both inputs to force multiple lengths
  n_half <- floor(nrow(fastq_input) / 2)

  fastq_input$Sequence[1:n_half] <- substr(fastq_input$Sequence[1:n_half], 1, 20)
  fastq_input$Quality[1:n_half] <- substr(fastq_input$Quality[1:n_half], 1, 20)

  reverse$Sequence[1:n_half] <- substr(reverse$Sequence[1:n_half], 1, 20)
  reverse$Quality[1:n_half] <- substr(reverse$Quality[1:n_half], 1, 20)

  # Run function
  suppressWarnings({
    merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                             reverse = reverse)
  })

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL

  # Check output is as expected
  expect_s3_class(merging_lengths_df, "tbl_df")
})

