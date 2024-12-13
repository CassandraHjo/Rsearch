test_that("error when wrong file_format", {

  files_dir <- test_path("testdata")
  file_format <- "fastx"

  expect_error(fastx_combine_files(files_dir = files_dir,
                                   output_file = NULL,
                                   file_ext = ".fa",
                                   file_format = file_format),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error when directory does not exist", {

  files_dir <- "some_directory"

  expect_error(fastx_combine_files(files_dir = files_dir),
               paste("Directory does not exist:", files_dir))
})

test_that("error when no files are found in specified folder", {

  files_dir <- test_path("testdata")
  file_ext <- ".fa"

  expect_error(fastx_combine_files(files_dir = files_dir),
               paste("No", file_ext, "files found in the specified folder:", files_dir))
})
