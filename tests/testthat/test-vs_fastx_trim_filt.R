test_that("error when wrong output_format", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fastx"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when output_format is fastq, and fastaout and fastaout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fastq"
  fastaout <- "some_file.fa"
  fastaout_rev <- "some_other_file.fa"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format,
                                  fastaout = fastaout,
                                  fastaout_rev = fastaout_rev),
               "When output_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
})

test_that("error when input is fasta and output_format is fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fasta_df.rds"))
  output_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  output_format = output_format),
               "Invalid output_format when input tibble is of type 'fasta'")
})

test_that("error when output_format is 'fasta', and fastqout and fastqout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fasta"
  fastqout <- "some_file.fq"
  fastqout_rev <- "some_other_file.fq"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format,
                                  fastqout = fastqout,
                                  fastqout_rev = fastqout_rev),
               "When output_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "R1_fasta_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fasta_df.rds"))
  output_format <- "fasta"
  fastaout <- "some_file.fa"
  fastaout_rev <- NULL

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format,
                                  fastaout = fastaout,
                                  fastaout_rev = fastaout_rev),
               "When 'reverse' is specified and output_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fastq"
  fastqout <- "some_file.fq"
  fastqout_rev <- NULL

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format,
                                  fastqout = fastqout,
                                  fastqout_rev = fastqout_rev),
               "When 'reverse' is specified and output_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
})

test_that("error when fastx_input has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))

  expect_error(vs_fastx_trim_filt(fastx_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(-Header)

  expect_error(vs_fastx_trim_filt(fastx_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is fasta tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fasta_df.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_fasta_df.rds"))

  output_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when reverse has incorrect columns if input is fastq tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fasta_df.rds"))

  R2 <- readRDS(test_path("testdata", "R2_fasta_df.rds")) |>
    dplyr::select(-Header)

  output_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format,
                                  maxee_rate = NULL,
                                  truncqual = NULL),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist when output_format is fastq", {

  fastx_input <- "some_file.fq"
  reverse <- test_path("testdata", "R2.fastqq")
  output_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format),
               paste("Cannot find input file:", fastx_input))
})

test_that("error when reverse file does not exist when output_format is fastq", {

  fastx_input <- test_path("testdata", "R2.fastq")
  reverse <- "some_file.fq"
  output_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format),
               paste("Cannot find reverse file:", reverse))
})

test_that("error when input file does not exist when output_format is fasta", {

  fastx_input <- "some_file.fa"
  reverse <- test_path("testdata", "R2.fasta")
  output_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format),
               paste("Cannot find input file:", fastx_input))
})

test_that("error when reverse file does not exist when output_format is fasta", {

  fastx_input <- test_path("testdata", "R2.fasta")
  reverse <- "some_file.fa"
  output_format <- "fasta"

  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format,
                                  maxee_rate = NULL,
                                  truncqual = NULL),
               paste("Cannot find reverse file:", reverse))
})

test_that("error if reverse input is of type fasta and output_format is fastq", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fasta_df.rds"))
  output_format <- "fastq"

  expect_error(vs_fastx_trim_filt(fastx_input = R1,
                                  reverse = R2,
                                  output_format = output_format),
               "Invalid output_format when input tibble is of type 'fasta'")
})

test_that("trim/filter fastq sequences from two files, and return two fastq files", {

  fastx_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  fastqout <- withr::local_tempfile()
  fastqout_rev <- withr::local_tempfile()
  output_format <- "fastq"
  truncee <- 0.01
  log_file <- withr::local_tempfile()
  vsearch_options <- c("")
  relabel <- "OTU"

  return_value <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                     reverse = reverse,
                                     fastqout = fastqout,
                                     fastqout_rev = fastqout_rev,
                                     output_format = output_format,
                                     truncee = truncee,
                                     log_file = log_file,
                                     vsearch_options = vsearch_options,
                                     relabel= relabel)

  expect_null(return_value)

  expect_true(file.exists(log_file))

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "R1_trim_filt_fq_files.fastq")))

  expect_equal(microseq::readFastq(fastqout_rev),
               microseq::readFastq(test_path("testdata", "output", "R2_trim_filt_fq_files.fastq")))
})

test_that("trim/filter fastq sequences from two files, and return fastq tibble", {

  fastx_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  output_format <- "fastq"
  truncee <- 0.01
  sample <- "sample1"

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format,
                                  truncee = truncee,
                                  sample = sample,
                                  relabel_sha1 = TRUE)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_fq_files.rds")))

})

test_that("trim/filter fastq sequences from one file, and return fastq tibble", {

  fastx_input <- test_path("testdata", "R1.fastq")
  output_format <- "fastq"
  truncee <- 0.01

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  output_format = output_format,
                                  truncee = truncee)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_R1_fq_file.rds")))

})

test_that("trim/filter fastq sequences from two tibbles, and return fastq tibble", {

  fastx_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))
  output_format <- "fastq"
  truncee <- 0.01

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format,
                                  truncee = truncee)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_fq_tibbles.rds")))
})

test_that("trim/filter fasta sequences from two files, and return two fasta files", {

  fastx_input <- test_path("testdata", "R1.fasta")
  reverse <- test_path("testdata", "R2.fasta")
  fastaout <- withr::local_tempfile()
  fastaout_rev <- withr::local_tempfile()
  output_format <- "fasta"
  maxee_rate <- NULL
  truncqual <- NULL
  truncee <- NULL
  maxlen <- 1000
  trunclen <- 150

  return_value <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                     reverse = reverse,
                                     fastaout = fastaout,
                                     fastaout_rev = fastaout_rev,
                                     output_format = output_format,
                                     maxee_rate = maxee_rate,
                                     truncqual = truncqual,
                                     truncee = truncee,
                                     maxlen = maxlen,
                                     trunclen = trunclen)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "R1_trim_filt_fa_files.fasta")))

  expect_equal(microseq::readFasta(fastaout_rev),
               microseq::readFasta(test_path("testdata", "output", "R2_trim_filt_fa_files.fasta")))
})

test_that("trim/filter fasta sequences from two tibbles, and return fasta tibble", {

  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  reverse <- microseq::readFasta(test_path("testdata", "R2.fasta"))
  output_format <- "fasta"

  maxee_rate <- NULL
  truncqual <- NULL
  truncee <- NULL
  maxlen <- 1000
  trunclen <- 150

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format,
                                  maxee_rate = maxee_rate,
                                  truncqual = truncqual,
                                  truncee = truncee,
                                  maxlen = maxlen,
                                  trunclen = trunclen)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_fa_tibbles.rds")))
})

test_that("trim/filter fastq sequences from two files, and return fastq tibble with stripping", {

  fastx_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  output_format <- "fastq"
  truncee <- 0.01
  stripright <- 10
  stripleft <- 10

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  reverse = reverse,
                                  output_format = output_format,
                                  truncee = truncee,
                                  stripright = stripright,
                                  stripleft = stripleft)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_fq_files_strip.rds")))

})

test_that("trim/filter fastq sequences from one file with size values, and return fastq tibble", {

  fastx_input <- test_path("testdata", "output", "R1_derep.fastq")
  output_format <- "fastq"
  truncee <- 0.01
  minsize <- 1
  maxsize <- 10
  minqual <- 1

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input,
                                  output_format = output_format,
                                  truncee = truncee,
                                  minsize = minsize,
                                  maxsize = maxsize,
                                  minqual = minqual)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_fq_file_size.rds")))
})

test_that("trim/filter fastq sequences from pe_df data frame", {

  fastx_input <- readRDS(test_path("testdata", "pe_df.rds"))

  trim_filt <- vs_fastx_trim_filt(fastx_input = fastx_input)

  expect_equal(trim_filt,
               readRDS(test_path("testdata", "output", "trim_filt_pe_df.rds")))

  attr(fastx_input, "reverse") <- NULL
  expect_error(vs_fastx_trim_filt(fastx_input = fastx_input),
               "fastx_input has class 'pe_df' but no 'reverse' attribute found.")
})

test_that("trim/filter fasta sequences with wrong parameters", {

  fastx_input <- test_path("testdata", "R1.fasta")
  reverse <- test_path("testdata", "R2.fasta")
  output_format <- "fasta"
  maxee_rate <- 0.5
  truncqual <- 2
  truncee <- 100
  truncee_rate <- 0.8

  warnings <- capture_warnings(
    vs_fastx_trim_filt(
      fastx_input = fastx_input,
      reverse = reverse,
      output_format = output_format,
      maxee_rate = maxee_rate,
      truncqual = truncqual,
      truncee = truncee,
      truncee_rate = truncee_rate
    )
  )

  expect_true(any(grepl("maxee_rate is ignored for FASTA input", warnings)))
  expect_true(any(grepl("truncqual is ignored for FASTA input", warnings)))
  expect_true(any(grepl("truncee is ignored for FASTA input", warnings)))
  expect_true(any(grepl("truncee_rate is ignored for FASTA input", warnings)))
})
