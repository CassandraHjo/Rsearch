test_that("error when wrong file_format", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastx"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                   reverse = R2,
                                   file_format = file_format),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error when file_format is 'fasta', and fastqout and fastqout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fasta"
  fastqout <- "some_file.fq"
  fastqout_rev <- "some_other_file.fq"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                   reverse = R2,
                                   file_format = file_format,
                                   fastqout = fastqout,
                                   fastqout_rev = fastqout_rev),
               "When file_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
})

test_that("error when file_format is 'fastq', and fastaout and fastaout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  fastaout <- "some_file.fa"
  fastaout_rev <- "some_other_file.fa"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                   reverse = R2,
                                   file_format = file_format,
                                   fastaout = fastaout,
                                   fastaout_rev = fastaout_rev),
               "When file_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fasta"
  fastaout <- "some_file.fa"
  fastaout_rev <- NULL

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                   reverse = R2,
                                   file_format = file_format,
                                   fastaout = fastaout,
                                   fastaout_rev = fastaout_rev),
               "When 'reverse' is specified and file_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  fastaout <- "some_file.fa"
  fastaout_rev <- NULL

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                   reverse = R2,
                                   file_format = file_format,
                                   fastqout = fastaout,
                                   fastqout_rev = fastaout_rev),
               "When 'reverse' is specified and file_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
})

test_that("error when fastx_input has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastx_trim_filt(fastx_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_fastx_trim_filt(fastx_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is fasta tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fasta_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fasta_dataframe.rds"))

  file_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  file_format = file_format),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when reverse has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fasta_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fasta_dataframe.rds")) %>%
    dplyr::select(-Header)

  file_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  file_format = file_format),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist when file format is fastq", {

  fastx_input <- "some_file.fq"
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  file_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                   reverse = reverse,
                                  file_format = file_format),
               paste("Cannot find input FASTQ file:", fastx_input))
})

test_that("error when reverse file does not exist when file format is fastq", {

  fastx_input <- test_path("testdata", "sample1", "R2_sample1.fq")
  reverse <- "some_file.fq"
  file_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                   reverse = reverse,
                                  file_format = file_format),
               paste("Cannot find reverse FASTQ file:", reverse))
})

test_that("error when input file does not exist when file format is fasta", {

  fastx_input <- "some_file.fa"
  reverse <- test_path("testdata", "sample1", "R2_sample1.fa")
  file_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  file_format = file_format),
               paste("Cannot find input FASTA file:", fastx_input))
})

test_that("error when reverse file does not exist when file format is fasta", {

  fastx_input <- test_path("testdata", "sample1", "R2_sample1.fa")
  reverse <- "some_file.fa"
  file_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  file_format = file_format),
               paste("Cannot find reverse FASTA file:", reverse))
})


test_that("filter fastq sequences from two files, and return two fastq files", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- withr::local_tempfile()
  fastqout_rev <- withr::local_tempfile()
  file_format <- "fastq"
  maxee_rate <- 0.01
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1
  log_file <- withr::local_tempfile()

  return_value <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                     reverse = reverse,
                                     fastqout = fastqout,
                                     fastqout_rev = fastqout_rev,
                                     file_format = file_format,
                                     maxee_rate = maxee_rate,
                                     minlen = minlen,
                                     trunclen = trunclen,
                                     log_file = log_file,
                                     threads = threads)

  expect_null(return_value)

  expect_true(file.exists(log_file))

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "R1_trim_filt_sample1.fq")))

  expect_equal(microseq::readFastq(fastqout_rev),
               microseq::readFastq(test_path("testdata", "output", "R2_trim_filt_sample1.fq")))
})

test_that("filter fastq sequences from two files, and return fastq tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- NULL
  fastqout_rev <- NULL
  file_format <- "fastq"
  maxee_rate <- 0.01
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1

  trim_filt_sample1 <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                          reverse = reverse,
                                          fastqout = fastqout,
                                          fastqout_rev = fastqout_rev,
                                          file_format = file_format,
                                          maxee_rate = maxee_rate,
                                          minlen = minlen,
                                          trunclen = trunclen,
                                          log_file = log_file,
                                          threads = threads)

  expect_equal(trim_filt_sample1,
               readRDS(test_path("testdata", "output", "trim_filt_sample1_fastq_files.rds")))

})

test_that("filter fasta sequences from two files, and return two fasta files", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fa")
  fastaout <- withr::local_tempfile()
  fastaout_rev <- withr::local_tempfile()
  file_format <- "fasta"
  maxlen <- 1000
  maxns <- 100
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                     reverse = reverse,
                                     fastaout = fastaout,
                                     fastaout_rev = fastaout_rev,
                                     file_format = file_format,
                                     maxlen = maxlen,
                                     minlen = minlen,
                                     maxns = maxns,
                                     trunclen = trunclen,
                                     log_file = log_file,
                                     threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "R1_trim_filt_sample1.fa")))

  expect_equal(microseq::readFasta(fastaout_rev),
               microseq::readFasta(test_path("testdata", "output", "R2_trim_filt_sample1.fa")))
})

test_that("filter fastq sequences from two tibbles, and return fastq tibble", {

  fastx_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastqout <- NULL
  fastqout_rev <- NULL
  file_format <- "fastq"
  maxee_rate <- 0.01
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1

  trim_filt_sample1 <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                          reverse = reverse,
                                          fastqout = fastqout,
                                          fastqout_rev = fastqout_rev,
                                          file_format = file_format,
                                          maxee_rate = maxee_rate,
                                          minlen = minlen,
                                          trunclen = trunclen,
                                          log_file = log_file,
                                          threads = threads)


  expect_equal(trim_filt_sample1,
               readRDS(test_path("testdata", "output", "trim_filt_sample1_fastq_tibbles.rds")))
})

test_that("filter fasta sequences from two tibbles, and return fasta tibble", {

  fastx_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  reverse <- microseq::readFasta(test_path("testdata", "sample1", "R2_sample1.fa"))
  fastaout <- NULL
  fastaout_rev <- NULL
  file_format <- "fasta"
  maxlen <- 1000
  maxns <- 100
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1

  trim_filt_sample1 <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                          reverse = reverse,
                                          fastaout = fastaout,
                                          fastaout_rev = fastaout_rev,
                                          file_format = file_format,
                                          maxlen = maxlen,
                                          minlen = minlen,
                                          maxns = maxns,
                                          trunclen = trunclen,
                                          log_file = log_file,
                                          threads = threads)


  expect_equal(trim_filt_sample1,
               readRDS(test_path("testdata", "output", "trim_filt_sample1_fasta_tibbles.rds")))
})

test_that("filter fastq sequences from one file, and return fastq tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- NULL
  fastqout <- NULL
  fastqout_rev <- NULL
  file_format <- "fastq"
  maxee_rate <- 0.01
  trunclen <- 150
  minlen <- 0
  log_file <- NULL
  threads <- 1

  trim_filt_sample1 <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                          reverse = reverse,
                                          fastqout = fastqout,
                                          fastqout_rev = fastqout_rev,
                                          file_format = file_format,
                                          maxee_rate = maxee_rate,
                                          minlen = minlen,
                                          trunclen = trunclen,
                                          log_file = log_file,
                                          threads = threads)

  expect_equal(trim_filt_sample1,
               readRDS(test_path("testdata", "output", "trim_filt_sample1_R1_fastq_file.rds")))

})
