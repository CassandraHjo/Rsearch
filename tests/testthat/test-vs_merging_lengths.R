test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))

  expect_error(vs_merging_lengths(fastq_input = R1,
                                  reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(-Header)

  expect_error(vs_merging_lengths(fastq_input = R1,
                                  reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("get merging lengths from merging two fastq files", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                           reverse = reverse)

  expected_df <- readRDS(test_path("testdata", "output", "merging_lengths_fq_files.rds"))

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(merging_lengths_df, expected_df)

  expect_error(vs_merging_lengths(fastq_input = fastq_input),
               "No reverse reads provided. Please supply reverse or use a 'pe_df' object.")
})

test_that("get merging lengths from merging two fastq tibbles", {

  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input,
                                           reverse = reverse,
                                           plot_title = FALSE)

  expected_df <- readRDS(test_path("testdata", "output", "merging_lengths_fq_tibbles.rds"))

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(merging_lengths_df, expected_df)

})

test_that("vs_merging_lengths creates histograms when R1 and R2 lengths vary", {
  # Load sample data
  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))

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

test_that("pe_df tibble can be merged", {

  fastq_input <- readRDS(test_path("testdata", "pe_df.rds"))

  merging_lengths_df <- vs_merging_lengths(fastq_input = fastq_input)

  expected_df <- readRDS(test_path("testdata", "output", "merging_lengths_pe_df.rds"))

  expect_s3_class(attr(merging_lengths_df, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(merging_lengths_df, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(merging_lengths_df, expected_df)

  attr(fastq_input, "reverse") <- NULL
  expect_error(vs_merging_lengths(fastq_input = fastq_input),
               "fastq_input has class 'pe_df' but no 'reverse' attribute found.")
})

