test_that("error when wrong file_format", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastx"

  expect_error(fastx_synchronize(file1 = file1,
                        file2 = file2,
                        file_format = file_format),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error if outputfiles are incorrectly specified", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  file1_out <- "output1.fq"
  file2_out <- NULL

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format,
                                 file1_out = file1_out,
                                 file2_out = file2_out),
               "Either both file1_out and file2_out must be NULL, or both must be specified.")

})

test_that("error if input is neither character or NULL", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  file1_out <- 1
  file2_out <- 1

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format,
                                 file1_out = file1_out,
                                 file2_out = file2_out),
               "file1_out must be a character string specifying the output file path.")

  file1_out <- "something"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format,
                                 file1_out = file1_out,
                                 file2_out = file2_out),
               "file2_out must be a character string specifying the output file path.")
})

test_that("error when file1 has incorrect columns if input is tibble and file_format = 'fastq'", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file1 FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when file2 has incorrect columns if input is tibble and file_format = 'fastq'", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file2 FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when file1 has incorrect columns if input is tibble and file_format = 'fasta'", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(Quality)
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  file_format <- "fasta"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file1 FASTA object must contain columns: Header and Sequence")
})

test_that("error when file2 has incorrect columns if input is tibble and file_format = 'fasta'", {

  file1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(Quality)

  file_format <- "fasta"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file2 FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file1 does not exist", {

  file1 <- "some_file.fq"
  file2 <- test_path("testdata", "sample1", "R2_sample1.fq")
  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               paste("Cannot find input file:", file1))
})

test_that("error when input file2 does not exist", {

  file1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  file2 <- "some_file.fq"
  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               paste("Cannot find input file:", file2))
})

test_that("two fasta files can be synchronized, and return fasta tibble with attribute", {

  file1 <- test_path("testdata", "sample1", "R1_sample1.fa")
  file2 <- test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fa")
  file_format <- "fasta"
  file1_out <- NULL
  file2_out <- NULL


  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "sync_file2")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_sample1_fasta.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_sample1_fasta.rds")), "sync_file2"))
})

test_that("two fastq files can be synchronized, and return fastq tibble with attribute", {

  file1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  file2 <- test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fq")
  file_format <- "fastq"
  file1_out <- NULL
  file2_out <- NULL


  sync_file1 <- fastx_synchronize(file1 = file1,
                                    file2 = file2,
                                    file_format = file_format,
                                    file1_out = file1_out,
                                    file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "sync_file2")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_sample1_fastq.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_sample1_fastq.rds")), "sync_file2"))
})

test_that("two fasta files can be synchronized, and return two fasta files", {

  file1 <- test_path("testdata", "sample1", "R1_sample1.fa")
  file2 <- test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fa")
  file_format <- "fasta"
  file1_out <- withr::local_tempfile()
  file2_out <- withr::local_tempfile()


  return_value <- fastx_synchronize(file1 = file1,
                                    file2 = file2,
                                    file_format = file_format,
                                    file1_out = file1_out,
                                    file2_out = file2_out)

  expect_null(return_value)

  expect_equal(microseq::readFasta(file1_out),
               microseq::readFasta(test_path("testdata", "output", "R1_sample1_sync.fa")))

  expect_equal(microseq::readFasta(file2_out),
               microseq::readFasta(test_path("testdata", "output", "R2_sample1_sync.fa")))
})

test_that("two fastq files can be synchronized, and return two fastq files", {

  file1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  file2 <- test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fq")
  file_format <- "fastq"
  file1_out <- withr::local_tempfile()
  file2_out <- withr::local_tempfile()


  return_value <- fastx_synchronize(file1 = file1,
                                    file2 = file2,
                                    file_format = file_format,
                                    file1_out = file1_out,
                                    file2_out = file2_out)

  expect_null(return_value)

  expect_equal(microseq::readFastq(file1_out),
               microseq::readFastq(test_path("testdata", "output", "R1_sample1_sync.fq")))

  expect_equal(microseq::readFastq(file2_out),
               microseq::readFastq(test_path("testdata", "output", "R2_sample1_sync.fq")))
})

test_that("two fasta tibbles can be synchronized, and return two fasta files", {

  file1 <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  file2 <- microseq::readFasta(test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fa"))
  file_format <- "fasta"
  file1_out <- withr::local_tempfile()
  file2_out <- withr::local_tempfile()


  return_value <- fastx_synchronize(file1 = file1,
                                    file2 = file2,
                                    file_format = file_format,
                                    file1_out = file1_out,
                                    file2_out = file2_out)

  expect_null(return_value)

  expect_equal(microseq::readFasta(file1_out),
               microseq::readFasta(test_path("testdata", "output", "R1_sample1_sync.fa")))

  expect_equal(microseq::readFasta(file2_out),
               microseq::readFasta(test_path("testdata", "output", "R2_sample1_sync.fa")))
})

test_that("two fastq tibbles can be synchronized, and return two fastq files", {

  file1 <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  file2 <- microseq::readFastq(test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fq"))
  file_format <- "fastq"
  file1_out <- withr::local_tempfile()
  file2_out <- withr::local_tempfile()


  return_value <- fastx_synchronize(file1 = file1,
                                    file2 = file2,
                                    file_format = file_format,
                                    file1_out = file1_out,
                                    file2_out = file2_out)

  expect_null(return_value)

  expect_equal(microseq::readFastq(file1_out),
               microseq::readFastq(test_path("testdata", "output", "R1_sample1_sync.fq")))

  expect_equal(microseq::readFastq(file2_out),
               microseq::readFastq(test_path("testdata", "output", "R2_sample1_sync.fq")))
})

test_that("two fasta tibbles can be synchronized, and return fasta tibble with attribute", {

  file1 <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  file2 <- microseq::readFasta(test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fa"))
  file_format <- "fasta"
  file1_out <- NULL
  file2_out <- NULL

  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "sync_file2")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_sample1_fasta.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_sample1_fasta.rds")), "sync_file2"))

})

test_that("two fastq tibbles can be synchronized, and return fastq tibble with attribute", {

  file1 <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  file2 <- microseq::readFastq(test_path("testdata", "sample1", "reduced_files", "R2_sample1_reduced.fq"))
  file_format <- "fastq"
  file1_out <- NULL
  file2_out <- NULL

  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "sync_file2")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_sample1_fastq.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_sample1_fastq.rds")), "sync_file2"))

})
