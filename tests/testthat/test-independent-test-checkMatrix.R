# --------------------------------------------
# Test checkMatrix function
#----------------------------------------------

test_that("Test checkMatrix with invalid input", code = {
  M <- 5
  expect_error(checkMatrix(M, nrows = 3, ncols = 3),
    info = "Test checkMatrix with invalid input"
  )
})


test_that("Test checkMatrix function for invalid nrows", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(checkRUEMatrix(M, nrows = 3, ncols = 3),
    info = "Test checkMatrix function for invalid nrows"
  )
})

test_that("Test checkMatrix function for invalid ncols", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(checkMatrix(M, nrows = 4, ncols = 2),
    info = "Test checkMatrix function for invalid ncols"
  )
})


test_that("Test checkMatrix function for invalid isType", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(checkMatrix(M, isType = "character"),
    info = "Test checkMatrix function for invalid isType"
  )
})

test_that("Test checkMatrix function for invalid isType", code = {
  M <- matrix(c("a", "b", "cc", "aa", "cc", "dd"), nrow = 2, byrow = TRUE)
  expect_error(checkMatrix(M, isType = "numeric"),
    info = "Test checkMatrix function for invalid isType"
  )
})

test_that("Test checkMatrix function for mixed data ", code = {
  M <- matrix(c(32, 42, 54, 16, 7, "a"), nrow = 2, byrow = TRUE)
  expect_error(checkMatrix(M, isType = "numeric"),
    info = "Test checkMatrix function for mixed data"
  )
})

test_that("Test checkMatrix function for data in valid interval",
  code = {
    M <- matrix(c(32, 42, 54, 16, 7, 34), nrow = 2, byrow = TRUE)
    expect_error(checkMatrix(M, interval = c(21, 45), isType = "numeric"),
      info = "Test checkMatrix function for data in valid interval"
    )
  }
)


test_that("Test checkMatrix function for integer data",
  code = {
    M <- matrix(c(32L, 42.5, 54L, 16L), nrow = 2, byrow = TRUE)
    expect_error(checkMatrix(M, isType = "integer"),
      info = "Test checkMatrix function for integer data"
    )
  }
)


test_that("Test checkMatrix function for valid data", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_invisible(checkMatrix(M, nrows = 4, ncols = 3))
})