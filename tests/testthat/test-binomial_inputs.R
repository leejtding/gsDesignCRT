testthat::context("binomial inputs")

testthat::test_that("test.testBinomial.x2", {
  testthat::expect_error(gsDesignCRT::testBinomial(x2 = "abc", x1 = 2, n1 = 2, n2 = 2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x2 = 3, x1 = 2, n1 = 2, n2 = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x2 = c(-1), x1 = 2, n1 = 2, n2 = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x2 = rep(2, 2), x1 = 2, n1 = rep(2, 3), n2 = 2),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.testBinomial.x1", {
  testthat::expect_error(gsDesignCRT::testBinomial(x1 = "abc", x2 = 2, n1 = 2, n2 = 2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x1 = 3, x2 = 2, n1 = 2, n2 = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x1 = c(-1), x2 = 2, n1 = 2, n2 = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(x1 = rep(2, 2), x2 = 2, n1 = rep(2, 3), n2 = 2),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.testBinomial.tol", {
  testthat::expect_error(gsDesignCRT::testBinomial(
    tol = "abc", x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::testBinomial(
    tol = 0, x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for out-of-range variable value")
})

testthat::test_that("test.testBinomial.scale", {
  testthat::expect_error(gsDesignCRT::testBinomial(
    scale = 1, x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::testBinomial(
    scale = "abc", x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    scale = rep("RR", 2), x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.testBinomial.n2", {
  testthat::expect_error(gsDesignCRT::testBinomial(n2 = "abc", x1 = 2, x2 = 2, n1 = 2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(n2 = 0, x1 = 2, x2 = 0, n1 = 2),
    info = "Checking for out-of-range variable value"
  )
})

testthat::test_that("test.testBinomial.n1", {
  testthat::expect_error(gsDesignCRT::testBinomial(n1 = "abc", x1 = 2, x2 = 2, n2 = 2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::testBinomial(n1 = 0, x1 = 0, x2 = 2, n2 = 2),
    info = "Checking for out-of-range variable value"
  )
})

testthat::test_that("test.testBinomial.delta0", {
  testthat::expect_error(gsDesignCRT::testBinomial(
    delta0 = "abc", x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::testBinomial(
    delta0 = c(-1), x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    delta0 = 1, x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    delta0 = rep(0.1, 2), x1 = 2, x2 = 2,
    n1 = rep(2, 3), n2 = 2
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.testBinomial.chisq", {
  testthat::expect_error(gsDesignCRT::testBinomial(
    chisq = "abc", x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::testBinomial(
    chisq = c(-1), x1 = 2, x2 = 2,
    n1 = 2, n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    chisq = 2, x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    chisq = rep(1, 2), x1 = 2, x2 = 2,
    n1 = rep(2, 3), n2 = 2
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.testBinomial.adj", {
  testthat::expect_error(gsDesignCRT::testBinomial(
    adj = "abc", x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::testBinomial(
    adj = c(-1), x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    adj = 2, x1 = 2, x2 = 2, n1 = 2,
    n2 = 2
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::testBinomial(
    adj = rep(1, 2), x1 = 2, x2 = 2,
    n1 = rep(2, 3), n2 = 2
  ), info = "Checking for incorrect variable length")
})
