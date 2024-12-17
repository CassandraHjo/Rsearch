test_that("error if vsearch_executable is wrong", {

  vsearch_executable <- "invalid_vsearch"

  expect_error(vsearch_available(vsearch_executable),
               paste("Cannot run", vsearch_executable, "from R, use set_vsearch_executable\\(\\) to set proper command to invoke vsearch"))
})
