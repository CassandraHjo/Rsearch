test_that("hello() prints Hello, world!", {
  expect_equal(hello(), c("Hello, world!"))
})

test_that("bye() prints Goodbye, world!", {
  expect_equal(bye(), c("Goodbye, world!"))
})
