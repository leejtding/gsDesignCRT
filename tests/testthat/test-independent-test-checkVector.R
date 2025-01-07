# Test checkVector function
#----------------------------------------------

test_that("Test checkVector for invalid value for length", code = {
  x <- c(1, 5, 2, 3)
  expect_error(checkVector(x, length = "e"),
    info = "Test checkVector for invalid value for length"
  )
})


test_that("Test checkVector length not matching vector length", code = {
  x <- c(1, 5, 2, 3)
  expect_error(checkVector(x, length = 2),
    info = "Test checkVector length not matching vector length"
  )
})


test_that("Test checkVector with correct value of the length parameter ", code = {
  x <- c(1, 5, 2, 3)
  expect_invisible(checkVector(x, length = 4))
})