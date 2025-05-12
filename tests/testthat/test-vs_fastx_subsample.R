test_that("error when wrong output_format", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  output_format <- "fastx"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when input is fasta and output_format is fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fasta_df.rds"))
  output_format <- "fastq"
  sample_size <- 20

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size),
               "Invalid output_format when input tibble is of type 'fasta'")
})

test_that("error if neither sample_size or sample_pct is specified", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  output_format <- "fastq"
  sample_size <- NULL
  sample_pct <- NULL

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size,
                                  sample_pct = sample_pct),
               "Either sample_size or sample_pct must be specified.")
})

test_that("error if both sample_size and sample_pct are specified", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  output_format <- "fastq"
  sample_size <- 10
  sample_pct <- 10.0

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size,
                                  sample_pct = sample_pct),
               "Only specify one of the following parameters, not both: sample_size, sample_pct ")
})

test_that("error when fastx_input has incorrect columns if input is fastq tibble", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  output_format <- "fastq"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is fasta tibble", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Header)

  output_format <- "fasta"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  output_format = output_format,
                                  sample_size = sample_size),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist", {

  fastx_input <- "some_file.fq"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  sample_size = sample_size),
               paste("Cannot find input file:", fastx_input))
})

test_that("subsample fastq file with size, and return fastq file", {

  fastx_input <- test_path("testdata", "R1.fastq")
  fastx_output <- withr::local_tempfile()
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     output_format = output_format,
                                     sample_size = sample_size,
                                     randseed = randseed
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_fq_file.fastq")))
})

test_that("subsample fastq file with size, and return fastq tibble", {

  fastx_input <- test_path("testdata", "R1.fastq")
  fastx_output <- NULL
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             output_format = output_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_fq_file_fq.rds")))
})

test_that("subsample fasta file with pct, and return fasta file", {

  fastx_input <- test_path("testdata", "R1.fasta")
  fastx_output <- withr::local_tempfile()
  output_format <- "fasta"
  sample_pct <- 10.0
  randseed <- 1

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     output_format = output_format,
                                     sample_pct = sample_pct,
                                     randseed = randseed
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_fa_file.fasta")))
})

test_that("subsample fasta file with pct, and return fasta tibble", {

  fastx_input <- test_path("testdata", "R1.fasta")
  fastx_output <- NULL
  output_format <- "fasta"
  sample_pct <- 10.0
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             output_format = output_format,
                                             sample_pct = sample_pct,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_fa_file_fa.rds")))
})

test_that("subsample fastq file with size, and return fastq file with relabeling", {

  fastx_input <- test_path("testdata", "R1.fastq")
  fastx_output <- withr::local_tempfile()
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1
  relabel <- "OTU"

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     output_format = output_format,
                                     sample_size = sample_size,
                                     randseed = randseed,
                                     relabel = relabel
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_fq_file_relabel.fastq")))
})

test_that("subsample fastq file with size, and return fastq file with sha1 relabeling", {

  fastx_input <- test_path("testdata", "R1.fastq")
  fastx_output <- withr::local_tempfile()
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1
  relabel_sha1 <- TRUE

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     output_format = output_format,
                                     sample_size = sample_size,
                                     randseed = randseed,
                                     relabel_sha1 = relabel_sha1
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_fq_file_sha1.fastq")))
})

test_that("subsample fastq tibble with size, and return fastq tibble", {

  fastx_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  fastx_output <- NULL
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             output_format = output_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_fq_file_size_fq.rds")))
})

test_that("subsample fasta tibble with size, and return fasta tibble", {

  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  fastx_output <- NULL
  output_format <- "fasta"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             output_format = output_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_fa_tibble_fa.rds")))
})

test_that("subsample fastq file with size, and return fastq tibble with vsearch_options", {

  fastx_input <- test_path("testdata", "R1.fastq")
  fastx_output <- NULL
  output_format <- "fastq"
  sample_size <- 100
  randseed <- 1
  sample <- "sample1"
  vsearch_options <- c("--relabel", "OTU")

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             output_format = output_format,
                                             sample_size = sample_size,
                                             randseed = randseed,
                                             sample = sample,
                                             vsearch_options = vsearch_options
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_fq_file_vs_fq.rds")))
})
