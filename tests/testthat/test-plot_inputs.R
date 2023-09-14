testthat::context("plot inputs")

testthat::test_that("test.plot.gsDesignCRT.plottype", {
  testthat::expect_error(gsDesignCRT::plot.gsDesignCRT(plottype = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::plot.gsDesignCRT(plottype = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::plot.gsDesignCRT(plottype = 8), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::plot.gsDesignCRT(plottype = rep(2, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.plot.gsProbability.plottype", {
  testthat::expect_error(gsDesignCRT::plot.gsProbability(plottype = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::plot.gsProbability(plottype = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::plot.gsProbability(plottype = 8), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::plot.gsProbability(plottype = rep(2, 2)), info = "Checking for incorrect variable length")
})
