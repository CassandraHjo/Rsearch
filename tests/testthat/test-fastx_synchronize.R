test_that("error when wrong file_format", {

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  file_format <- "fastx"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "Invalid file_format. The files must be a fasta or fastq.")
})

test_that("error if outputfiles are incorrectly specified", {

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
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

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
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

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))

  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file1 FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when file2 has incorrect columns if input is tibble and file_format = 'fastq'", {

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(-Header)

  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file2 FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when file1 has incorrect columns if input is tibble and file_format = 'fasta'", {

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Quality)
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))

  file_format <- "fasta"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file1 FASTA object must contain columns: Header and Sequence")
})

test_that("error when file2 has incorrect columns if input is tibble and file_format = 'fasta'", {

  file1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  file2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(Quality)

  file_format <- "fasta"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               "file2 FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file1 does not exist", {

  file1 <- "some_file.fq"
  file2 <- test_path("testdata", "R2.fastq")
  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               paste("Cannot find input file:", file1))
})

test_that("error when input file2 does not exist", {

  file1 <- test_path("testdata", "R1.fastq")
  file2 <- "some_file.fq"
  file_format <- "fastq"

  expect_error(fastx_synchronize(file1 = file1,
                                 file2 = file2,
                                 file_format = file_format),
               paste("Cannot find input file:", file2))
})

test_that("two fasta files can be synchronized, and return fasta tibble with attribute", {

  file1 <- test_path("testdata", "R1.fasta")
  file2 <- test_path("testdata", "R2_reduced.fasta")
  file_format <- "fasta"
  file1_out <- NULL
  file2_out <- NULL


  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "reverse")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_fa_files_fa_tibble.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_fa_files_fa_tibble.rds")), "reverse"))
})

test_that("two fastq files can be synchronized, and return fastq tibble with attribute", {

  file1 <- test_path("testdata", "R1.fastq")
  file2 <- test_path("testdata", "R2_reduced.fastq")
  file_format <- "fastq"
  file1_out <- NULL
  file2_out <- NULL


  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "reverse")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_fq_files_fq_tibble.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_fq_files_fq_tibble.rds")), "reverse"))
})

test_that("two fasta files can be synchronized, and return two fasta files", {

  file1 <- test_path("testdata", "R1.fasta")
  file2 <- test_path("testdata", "R2_reduced.fasta")
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
               microseq::readFasta(test_path("testdata", "output", "R1_sync.fasta")))

  expect_equal(microseq::readFasta(file2_out),
               microseq::readFasta(test_path("testdata", "output", "R2_sync.fasta")))
})

test_that("two fastq files can be synchronized, and return two fastq files", {

  file1 <- test_path("testdata", "R1.fastq")
  file2 <- test_path("testdata", "R2_reduced.fastq")
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
               microseq::readFastq(test_path("testdata", "output", "R1_sync.fastq")))

  expect_equal(microseq::readFastq(file2_out),
               microseq::readFastq(test_path("testdata", "output", "R2_sync.fastq")))
})

test_that("two fasta tibbles can be synchronized, and return two fasta files", {

  file1 <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  file2 <- microseq::readFasta(test_path("testdata", "R2_reduced.fasta"))
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
               microseq::readFasta(test_path("testdata", "output", "R1_sync.fasta")))

  expect_equal(microseq::readFasta(file2_out),
               microseq::readFasta(test_path("testdata", "output", "R2_sync.fasta")))
})

test_that("two fastq tibbles can be synchronized, and return two fastq files", {

  file1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  file2 <- microseq::readFastq(test_path("testdata", "R2_reduced.fastq"))
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
               microseq::readFastq(test_path("testdata", "output", "R1_sync.fastq")))

  expect_equal(microseq::readFastq(file2_out),
               microseq::readFastq(test_path("testdata", "output", "R2_sync.fastq")))
})

test_that("two fasta tibbles can be synchronized, and return fasta tibble with attribute", {

  file1 <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  file2 <- microseq::readFasta(test_path("testdata", "R2_reduced.fasta"))
  file_format <- "fasta"
  file1_out <- NULL
  file2_out <- NULL

  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "reverse")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_fa_files_fa_tibble.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_fa_files_fa_tibble.rds")), "reverse"))
})

test_that("two fastq tibbles can be synchronized, and return fastq tibble with attribute", {

  file1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  file2 <- microseq::readFastq(test_path("testdata", "R2_reduced.fastq"))
  file_format <- "fastq"
  file1_out <- NULL
  file2_out <- NULL

  sync_file1 <- fastx_synchronize(file1 = file1,
                                  file2 = file2,
                                  file_format = file_format,
                                  file1_out = file1_out,
                                  file2_out = file2_out)

  sync_file2 <- attr(sync_file1, "reverse")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_fq_files_fq_tibble.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_fq_files_fq_tibble.rds")), "reverse"))
})

test_that("pe_df tibble can be synchronized", {

  file1 <- readRDS(test_path("testdata", "pe_df.rds"))

  sync_file1 <- fastx_synchronize(file1 = file1)

  sync_file2 <- attr(sync_file1, "reverse")

  expect_equal(sync_file1,
               readRDS(test_path("testdata", "output", "sync_pe_df.rds")))

  expect_equal(sync_file2,
               attr(readRDS(test_path("testdata", "output", "sync_pe_df.rds")), "reverse"))

  attr(file1, "reverse") <- NULL
  expect_error(fastx_synchronize(file1 = file1),
               "file1 has class 'pe_df' but no 'reverse' attribute found.")
})
