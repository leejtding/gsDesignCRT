testthat::context("nbinomial inputs")

testthat::test_that("test.nBinomial.alpha", {
  testthat::expect_error(gsDesignCRT::nBinomial(alpha = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(alpha = 1, p1 = 0.1, p2 = 0.2, sided = 1),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(alpha = 0.5, p1 = 0.1, p2 = 0.2, sided = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(alpha = 0, p1 = 0.1, p2 = 0.2, sided = 2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(alpha = rep(0.1, 2), p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.nBinomial.beta", {
  testthat::expect_error(gsDesignCRT::nBinomial(beta = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(
    beta = 0.7, p1 = 0.1, p2 = 0.2, sided = 1,
    alpha = 0.3
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(beta = 0, p1 = 0.1, p2 = 0.2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(beta = rep(0.1, 2), p1 = 0.1, p2 = rep(
    0.2,
    3
  )), info = "Checking for incorrect variable length")
})

testthat::test_that("test.nBinomial.delta0", {
  testthat::expect_error(gsDesignCRT::nBinomial(delta0 = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(delta0 = 0, p1 = 0.1, p2 = 0.1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(delta0 = 0.1, p1 = 0.2, p2 = 0.1),
    info = "Checking for out-of-range variable value"
  )
})

testthat::test_that("test.nBinomial.outtype", {
  testthat::expect_error(gsDesignCRT::nBinomial(outtype = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(outtype = 0, p1 = 0.1, p2 = 0.2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(outtype = 4, p1 = 0.1, p2 = 0.2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(outtype = rep(1, 2), p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.nBinomial.p1", {
  testthat::expect_error(gsDesignCRT::nBinomial(p1 = "abc", p2 = 0.1), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::nBinomial(p1 = 0, p2 = 0.1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(p1 = 1, p2 = 0.1), info = "Checking for out-of-range variable value")
})

testthat::test_that("test.nBinomial.p2", {
  testthat::expect_error(gsDesignCRT::nBinomial(p2 = "abc", p1 = 0.1), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::nBinomial(p2 = 0, p1 = 0.1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(p2 = 1, p1 = 0.1), info = "Checking for out-of-range variable value")
})

testthat::test_that("test.nBinomial.ratio", {
  testthat::expect_error(gsDesignCRT::nBinomial(ratio = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(ratio = 0, p1 = 0.1, p2 = 0.2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(ratio = rep(0.5, 2), p1 = 0.1, p2 = rep(
    0.2,
    3
  )), info = "Checking for incorrect variable length")
})

testthat::test_that("test.nBinomial.scale", {
  testthat::expect_error(gsDesignCRT::nBinomial(scale = 1, p1 = 0.1, p2 = 0.2), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::nBinomial(scale = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(scale = rep("RR", 2), p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.nBinomial.sided", {
  testthat::expect_error(gsDesignCRT::nBinomial(sided = "abc", p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsDesignCRT::nBinomial(sided = 0, p1 = 0.1, p2 = 0.2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(sided = 3, p1 = 0.1, p2 = 0.2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::nBinomial(sided = rep(1, 2), p1 = 0.1, p2 = 0.2),
    info = "Checking for incorrect variable length"
  )
})
