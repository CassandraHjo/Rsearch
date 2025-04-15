test_that("check_vsearch_status() does nothing if status is NULL", {

  fake_output <- c("All good", "Processing done")
  expect_silent(check_vsearch_status(fake_output, args = c("--arg1", "value")))
})

test_that("check_vsearch_status() stops on error by default", {

  fake_output <- c("Some VSEARCH output", "Fatal error: bad argument")
  attr(fake_output, "status") <- 1
  attr(fake_output, "errmsg") <- "Resource temporarily unavailable"

  expect_error(
    check_vsearch_status(fake_output, args = c("--badarg")),
    regexp = "VSEARCH execution failed"
  )
})
