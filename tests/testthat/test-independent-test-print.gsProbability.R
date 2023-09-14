#testing print.gsProbability for gsDesignCRT and gsBinomialExact class.

testthat::test_that('Test: object of class gsDesignCRT - plot for y',
                    code = {
  x <- gsDesignCRT()
  y <- gsProbability(d = x, theta = x$delta * seq(0, 2, .25))
  local_edition(3)  #use the 3rd edition of the testthat package
  expect_snapshot_output(x = print.gsProbability(y))
})


testthat::test_that('Test: object of class gsDesignCRT  - plot for z',
                    code = {
  x <- gsDesignCRT()
  z <- gsProbability(k = 3, a = x$lower$bound, b = x$upper$bound, n.I = x$n.I, 
                     theta = x$delta * seq(0, 2, .25))
  local_edition(3)
  expect_snapshot_output(x = print.gsProbability(z))
})


testthat::test_that('Test: object of class gsBinomialExact',
                    code = {
  x <- gsBinomialExact()
  local_edition(3)
  expect_snapshot_output(x = print.gsProbability(x))
})
