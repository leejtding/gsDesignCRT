context("base tests")

### Utility functions for gsDesign stress tests

"alpha.beta.range.util" <- function(alpha, beta, test_type, sf) {
  no.err <- TRUE

  for (a in alpha) {
    for (b in beta) {
      if (b < 1 - a - 0.1) {
        res <- try(gsDesignCRT(test_type = test_type,
                               alpha = a, beta = b,
                               alpha_sf = sf, beta_sf = sf))

        if (inherits(res, "try-error")) {
          no.err <- FALSE
        }
      }
    }
  }

  no.err
}

"param.range.util" <- function(param, test_type, sf) {
  no.err <- TRUE

  for (p in param)
  {
    res <- try(gsDesignCRT(test_type = test_type,
                           alpha_sf = sf, alpha_sfpar = p,
                           beta_sf = sf, beta_sfpar = p))

    if (inherits(res, "try-error")) {
      no.err <- FALSE
    }
  }

  no.err
}

a1 <- round(seq(from = 0.05, to = 0.95, by = 0.05), 2)
a2 <- round(seq(from = 0.05, to = 0.45, by = 0.05), 2)
b <- round(seq(from = 0.05, to = 0.95, by = 0.05), 2)

# nu: sfExponential parameter
nu <- round(seq(from = 0.1, to = 1.5, by = 0.1), 1)

# rho: sfPower parameter
rho <- round(seq(from = 1, to = 15, by = 1), 0)

# gamma: sfHSD parameter
gamma <- round(seq(from = -5, to = 5, by = 1), 0)

######### gsUpperCRT tests #########

test_that("test.gsUpperCRT.theta", {
  expect_error(gsUpperCRT(theta = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(theta = rep(1, 2), I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gsUpperCRT(
    theta = -1, I = 1, falsepos = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsUpperCRT.I", {
  expect_error(gsUpperCRT(I = "abc", falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(I = 0, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsUpperCRT(I = rep(1, 2), falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gsUpperCRT.a", {
  expect_error(gsUpperCRT(a = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(a = rep(0.5, 2), I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gsUpperCRT(a = -0.5, I = 1, falsepos = 0.5, sides = 2),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsUpperCRT.falsepos", {
  expect_error(gsUpperCRT(falsepos = "abc", I = 1),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(falsepos = 1, I = 1),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsUpperCRT(falsepos = rep(0.5, 2), I = 1),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gsUpperCRT.sides", {
  expect_error(gsUpperCRT(sides = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(sides = 0, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsUpperCRT.tol", {
  expect_error(gsUpperCRT(tol = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(tol = 0, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsUpperCRT.r", {
  expect_error(gsUpperCRT(r = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsUpperCRT(r = 0, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsUpperCRT(r = 81, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsUpperCRT(r = rep(1, 2), I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
})

######## gsLowerCRT tests #########

test_that("test.gsLowerCRT.theta", {
  expect_error(gsLowerCRT(theta = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(theta = rep(1, 2), I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gsLowerCRT(
    theta = -1, I = 1, falseneg = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsLowerCRT.I", {
  expect_error(gsLowerCRT(I = "abc", falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(I = 0, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsLowerCRT(I = rep(1, 2), falseneg = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gsLowerCRT.b", {
  expect_error(gsLowerCRT(b = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(b = rep(0.5, 2), I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gsLowerCRT(b = -0.5, I = 1, falseneg = 0.5, sides = 2),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsLowerCRT.falseneg", {
  expect_error(gsLowerCRT(falseneg = "abc", I = 1),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(falseneg = 1, I = 1),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsLowerCRT(falseneg = rep(0.5, 2), I = 1),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gsLowerCRT.sides", {
  expect_error(gsLowerCRT(sides = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(sides = 0, I = 1, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsLowerCRT.binding", {
  expect_error(gsLowerCRT(binding = 2, I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(binding = "TRUE", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
})

test_that("test.gsLowerCRT.tol", {
  expect_error(gsLowerCRT(tol = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(tol = 0, I = 1, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsLowerCRT.r", {
  expect_error(gsLowerCRT(r = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsLowerCRT(r = 0, I = 1, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsLowerCRT(r = 81, I = 1, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsLowerCRT(r = rep(1, 2), I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gsBoundsCRT.theta", {
  expect_error(gsBoundsCRT(
    theta = "abc", I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(
    theta = rep(1, 2), I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable length"
  )
  expect_error(gsBoundsCRT(
    theta = -1, I = 1, falseneg = 0.5, falsepos = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

######## gsBoundsCRT tests #########

test_that("test.gsBoundsCRT.I", {
  expect_error(gsBoundsCRT(I = "abc", falseneg = 0.5, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(I = 0, falseneg = 0.5, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsBoundsCRT(
    I = rep(1, 2), falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gsBoundsCRT.falsepos", {
  expect_error(gsBoundsCRT(falsepos = "abc", I = 1, falseneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(falsepos = 1, I = 1, falseneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsBoundsCRT(
    falsepos = rep(0.5, 2), I = 1, falseneg = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gsBoundsCRT.falseneg", {
  expect_error(gsBoundsCRT(falseneg = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(falseneg = 1, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gsBoundsCRT(
    falseneg = rep(0.5, 2), I = 1, falsepos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gsBoundsCRT.sides", {
  expect_error(gsBoundsCRT(
    sides = "abc", I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(
    sides = 0, I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsBoundsCRT.binding", {
  expect_error(gsBoundsCRT(
    binding = 2, I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(
    binding = "TRUE", I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
})

test_that("test.gsBoundsCRT.tol", {
  expect_error(gsBoundsCRT(
    tol = "abc", I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(
    tol = 0, I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gsBoundsCRT.r", {
  expect_error(gsBoundsCRT(
    r = "abc", I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gsBoundsCRT(
    r = 0, I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
  expect_error(gsBoundsCRT(
    r = 81, I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
  expect_error(gsBoundsCRT(
    r = rep(1, 2), I = 1, falseneg = 0.5, falsepos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

######### gsDesignCRT tests #########

test_that("test.gsDesignCRT.k", {
  expect_error(gsDesignCRT(k = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(k = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(k = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(k = seq(2)),
               info = "Checking for incorrect variable length")
  expect_error(gsDesignCRT(
    k = 3, alpha_sf = sfPoints, alpha_sfpar = c(0.05, 0.1, 0.15, 0.2, 1)
  ), info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.outcome.type", {
  expect_error(gsDesignCRT(outcome_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(outcome_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(outcome_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(outcome_type = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(outcome_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.test.type", {
  expect_error(gsDesignCRT(test_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(test_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(test_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(test_type = 6),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(test_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.test.sides", {
  expect_error(gsDesignCRT(test_sides = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(test_sides = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(test_sides = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(test_sides = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(test_sides = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.size.type", {
  expect_error(gsDesignCRT(size_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(size_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(size_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(size_type = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(size_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.timing.type", {
  expect_error(gsDesignCRT(timing_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(timing_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(timing_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 4),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.recruit.type", {
  expect_error(gsDesignCRT(recruit_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(recruit_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(recruit_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(recruit_type = 4),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(recruit_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.delta", {
  expect_error(gsDesignCRT(delta = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(delta = -1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(delta = 2, outcome_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(delta = rep(0.1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.sigma_vec", {
  expect_error(gsDesignCRT(outcome_type = 1, sigma_vec = c("abc", 1)),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(outcome_type = 1, sigma_vec = 1),
               info = "Checking for incorrect variable length")
  expect_error(gsDesignCRT(outcome_type = 1, sigma_vec = seq(3)),
               info = "Checking for incorrect variable length")
  expect_error(gsDesignCRT(outcome_type = 1, sigma_vec = c(1, -1)),
               info = "Checking for out-of-range variable value")
})

test_that("test.gsDesignCRT.p_vec", {
  expect_error(gsDesignCRT(outcome_type = 2, p_vec = c("abc", 0.5)),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(outcome_type = 2, p_vec = 0.5),
               info = "Checking for incorrect variable length")
  expect_error(gsDesignCRT(outcome_type = 2, p_vec = rep(0.5, 3)),
               info = "Checking for incorrect variable length")
  expect_error(gsDesignCRT(outcome_type = 2, p_vec = c(0.5, -0.5)),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(outcome_type = 2, p_vec = c(1.5, 0.5)),
               info = "Checking for out-of-range variable value")
})

test_that("test.gsDesignCRT.rho", {
  expect_error(gsDesignCRT(rho = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(rho = -0.5),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(rho = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(rho = rep(1, 3)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.alpha", {
  expect_error(gsDesignCRT(alpha = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(alpha = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(alpha = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(alpha = rep(0.5, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.beta", {
  expect_error(gsDesignCRT(beta = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(beta = 0.5, alpha = 0.5),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(beta = 1, alpha = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(beta = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(beta = rep(0.1, 2), alpha = 0.5),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.m", {
  expect_error(gsDesignCRT(m = "abc", size_type = 2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(m = -1, size_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(m = rep(2, 3), size_type = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.m.alloc", {
  expect_error(gsDesignCRT(m_alloc = "abc", size_type = 2),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(m_alloc = c(0.5, 0.6), size_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(m_alloc = rep(0.5, 3), size_type = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.n", {
  expect_error(gsDesignCRT(n = "abc", size_type = 1),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(n = -1, size_type = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(n = rep(2, 3), size_type = 1),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.n.cv", {
  expect_error(gsDesignCRT(n_cv = "abc", size_type = 1),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(n_cv = -1, size_type = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(n_cv = rep(2, 3), size_type = 1),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.info.timing", {
  expect_error(gsDesignCRT(timing_type = 1, info_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(timing_type = 1, info_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 1, info_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 1, info_timing = c(0.1, 1.1), k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 1, info_timing = c(0.5, 0.1), k = 2),
               info = "Checking for incorrect variable specification")
  expect_error(gsDesignCRT(timing_type = 1, info_timing = c(0.1, 0.5, 1),
                           k = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.m.timing", {
  expect_error(gsDesignCRT(timing_type = 3, m_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(timing_type = 3, m_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3, m_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           m_timing = matrix(c(0.5, 0.5, 1, 2), nrow = 2),
                           k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           m_timing = matrix(c(0.5, 0.5, 0.25, 0.25, 1, 1),
                                             nrow = 2),
                           k = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           m_timing = matrix(c(0.25, 0.25, 0.5, 0.5, 1, 2),
                                             nrow = 2),
                           k = 2),
               info = "Checking for incorrect variable specification")
})

test_that("test.gsDesignCRT.n.timing", {
  expect_error(gsDesignCRT(timing_type = 3, n_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(timing_type = 3, n_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3, n_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           n_timing = matrix(c(0.5, 0.5, 1, 2), nrow = 2),
                           k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           n_timing = matrix(c(0.5, 0.5, 0.25, 0.25, 1, 1),
                                             nrow = 2),
                           k = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(timing_type = 3,
                           n_timing = matrix(c(0.25, 0.25, 0.5, 0.5, 1, 2),
                                             nrow = 2),
                           k = 2),
               info = "Checking for incorrect variable specification")
})

test_that("test.gsDesignCRT.alpha_sf", {
  expect_error(gsDesignCRT(alpha_sf = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(alpha_sf = rep(sfLDOF, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.beta_sf", {
  expect_error(gsDesignCRT(beta_sf = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(beta_sf = rep(sfLDOF, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.r", {
  expect_error(gsDesignCRT(r = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(r = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(r = 81),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(r = rep(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsDesignCRT.tol", {
  expect_error(gsDesignCRT(tol = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsDesignCRT(tol = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(tol = 0.10000001),
               info = "Checking for out-of-range variable value")
  expect_error(gsDesignCRT(tol = rep(0.1, 2)),
               info = "Checking for incorrect variable length")
})

######### Stress tests for spending functions #########

test_that("test.stress.sfExp.type1", {
  no_errors <- param.range.util(param = nu, test_type = 1, sf = sfExponential)
  expect_true(no_errors, info = "Type 1 sfExponential stress test")
})

test_that("test.stress.sfExp.type2", {
  no_errors <- param.range.util(param = nu, test_type = 2, sf = sfExponential)
  expect_true(no_errors, info = "Type 2 sfExponential stress test")
})

test_that("test.stress.sfExp.type3", {
  no_errors <- param.range.util(param = nu, test_type = 3, sf = sfExponential)
  expect_true(no_errors, info = "Type 3 sfExponential stress test")
})

test_that("test.stress.sfExp.type4", {
  no_errors <- param.range.util(param = nu, test_type = 4, sf = sfExponential)
  expect_true(no_errors, info = "Type 4 sfExponential stress test")
})

test_that("test.stress.sfExp.type5", {
  no_errors <- param.range.util(param = nu, test_type = 5, sf = sfExponential)
  expect_true(no_errors, info = "Type 5 sfExponential stress test")
})

test_that("test.stress.sfHSD.type1", {
  no_errors <- param.range.util(param = gamma, test_type = 1, sf = sfHSD)
  expect_true(no_errors, info = "Type 1 sfHSD stress test")
})

test_that("test.stress.sfHSD.type2", {
  no_errors <- param.range.util(param = gamma, test_type = 2, sf = sfHSD)
  expect_true(no_errors, info = "Type 2 sfHSD stress test")
})

test_that("test.stress.sfHSD.type3", {
  no_errors <- param.range.util(param = gamma, test_type = 3, sf = sfHSD)
  expect_true(no_errors, info = "Type 3 sfHSD stress test")
})

test_that("test.stress.sfHSD.type4", {
  no_errors <- param.range.util(param = gamma, test_type = 4, sf = sfHSD)
  expect_true(no_errors, info = "Type 4 sfHSD stress test")
})

test_that("test.stress.sfHSD.type5", {
  no_errors <- param.range.util(param = gamma, test_type = 5, sf = sfHSD)
  expect_true(no_errors, info = "Type 5 sfHSD stress test")
})

test_that("test.stress.sfLDOF.type1", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 1, sf = sfLDOF
  )
  expect_true(no_errors, info = "Type 1 LDOF stress test")
})

test_that("test.stress.sfLDOF.type2", {
  no_errors <- alpha.beta.range.util(
    alpha = a2, beta = b,
    test_type = 2, sf = sfLDOF
  )
  expect_true(no_errors, info = "Type 2 LDOF stress test")
})

test_that("test.stress.sfLDOF.type3", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 3, sf = sfLDOF
  )
  expect_true(no_errors, info = "Type 3 LDOF stress test")
})

test_that("test.stress.sfLDOF.type4", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 4, sf = sfLDOF
  )
  expect_true(no_errors, info = "Type 4 LDOF stress test")
})

test_that("test.stress.sfLDOF.type5", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 5, sf = sfLDOF
  )
  expect_true(no_errors, info = "Type 5 LDOF stress test")
})

test_that("test.stress.sfLDPocock.type1", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 1, sf = sfLDPocock
  )
  expect_true(no_errors, info = "Type 1 LDPocock stress test")
})

test_that("test.stress.sfLDPocock.type2", {
  no_errors <- alpha.beta.range.util(
    alpha = a2, beta = b,
    test_type = 2, sf = sfLDPocock
  )
  expect_true(no_errors, info = "Type 2 LDPocock stress test")
})

test_that("test.stress.sfLDPocock.type3", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 3, sf = sfLDPocock
  )
  expect_true(no_errors, info = "Type 3 LDPocock stress test")
})

test_that("test.stress.sfLDPocock.type4", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 4, sf = sfLDPocock
  )
  expect_true(no_errors, info = "Type 4 LDPocock stress test")
})

test_that("test.stress.sfLDPocock.type5", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 5, sf = sfLDPocock
  )
  expect_true(no_errors, info = "Type 5 LDPocock stress test")
})

test_that("test.stress.sfPower.type1", {
  no_errors <- param.range.util(param = rho, test_type = 1, sf = sfPower)
  expect_true(no_errors, info = "Type 1 sfPower stress test")
})

test_that("test.stress.sfPower.type2", {
  no_errors <- param.range.util(param = rho, test_type = 2, sf = sfPower)
  expect_true(no_errors, info = "Type 2 sfPower stress test")
})

test_that("test.stress.sfPower.type3", {
  no_errors <- param.range.util(param = rho, test_type = 3, sf = sfPower)
  expect_true(no_errors, info = "Type 3 sfPower stress test")
})

test_that("test.stress.sfPower.type4", {
  no_errors <- param.range.util(param = rho, test_type = 4, sf = sfPower)
  expect_true(no_errors, info = "Type 4 sfPower stress test")
})

test_that("test.stress.sfPower.type5", {
  no_errors <- param.range.util(param = rho, test_type = 5, sf = sfPower)
  expect_true(no_errors, info = "Type 5 sfPower stress test")
})

######## gsProbabilityCRT tests #########

test_that("test.gsProbabilityCRT.theta", {
  expect_error(gsProbabilityCRT(theta = "abc"),
               info = "Checking for incorrect variable type")
})

test_that("test.gsProbabilityCRT.I", {
  expect_error(gsProbabilityCRT(I = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsProbabilityCRT(I = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsProbabilityCRT(I = c(2, 1), a = rep(0, 2), b = rep(1, 2)),
               info = "Checking for out-of-order input sequence")
  expect_error(gsProbabilityCRT(I = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsProbabilityCRT.a", {
  expect_error(gsProbabilityCRT(a = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsProbabilityCRT(a = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsProbabilityCRT.b", {
  expect_error(gsProbabilityCRT(b = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsProbabilityCRT(b = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gsProbabilityCRT.sides", {
  expect_error(gsProbabilityCRT(sides = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsProbabilityCRT(sides = 0),
               info = "Checking for out-of-range variable value")
})

test_that("test.gsProbability.r", {
  expect_error(gsProbabilityCRT(r = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gsProbabilityCRT(r = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gsProbabilityCRT(r = 81),
               info = "Checking for out-of-range variable value")
  expect_error(gsProbabilityCRT(r = rep(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.Deming.gsProb", {
  w <- sum(gsProbabilityCRT(theta = 0, I = 1:2,
                            a = stats::qnorm(0.025) * c(1, 1),
                            b = stats::qnorm(0.975) * c(1, 1))$upper$prob)
  x <- sum(gsProbabilityCRT(theta = 0, I = 1:4,
                            a = stats::qnorm(0.025) * rep(1, 4),
                            b = stats::qnorm(0.975) * rep(1, 4))$upper$prob)
  y <- sum(gsProbabilityCRT(theta = 0, I = 1:10,
                            a = stats::qnorm(0.025) * array(1, 10),
                            b = stats::qnorm(0.975) * array(1, 10))$upper$prob)
  z <- sum(gsProbabilityCRT(theta = 0, I = 1:20,
                            a = stats::qnorm(0.025) * array(1, 20),
                            b = stats::qnorm(0.975) * array(1, 20))$upper$prob)
  expect_equal(0.042, round(w, 3), info = "Checking Type I error, k = 2")
  expect_equal(0.063, round(x, 3), info = "Checking Type I error, k = 4")
  expect_equal(0.097, round(y, 3), info = "Checking Type I error, k = 10")
  expect_equal(0.124, round(z, 3), info = "Checking Type I error, k = 20")
})

######### Spending function parameter tests #########

test_that("test.sfcauchy.param", {
  expect_error(sfcauchy(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sfcauchy(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sfcauchy.param ", {
  expect_error(sfcauchy(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfcauchy(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfexp.param", {
  expect_error(sfexp(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfexp(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sfexp(param = 0),
               info = "Checking for out-of-range variable value")
  expect_error(sfexp(param = 11),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfHSD.param", {
  expect_error(sfHSD(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfHSD(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sfHSD(param = -41),
               info = "Checking for out-of-range variable value")
  expect_error(sfHSD(param = 41),
               info = "Checking for out-of-range variable value")
})

test_that("test.sflogistic.param", {
  expect_error(sflogistic(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sflogistic(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sflogistic.param ", {
  expect_error(sflogistic(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sflogistic(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})


test_that("test.sfnorm.param", {
  expect_error(sfnorm(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sfnorm(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sfnorm.param ", {
  expect_error(sfnorm(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfnorm(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfpower.param", {
  expect_error(sfpower(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfpower(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sfpower(param = -1),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfTDist.param", {
  expect_error(sfTDist(param = rep(1, 4)),
               info = "Checking for incorrect variable length")
  expect_error(sfTDist(param = c(1, 0, 1)),
               info = "Checking for out-of-range variable value")
  expect_error(sfTDist(param = c(1, 1, 0.5)),
               info = "Checking for out-of-range variable value")
  expect_error(sfTDist(param = 1, 1:3 / 4, c(
    0.25, 0.5, 0.75,
    0.1, 0.2, 0.3
  )), info = "Checking for out-of-range variable value")
})

test_that("test.sfTDist.param ", {
  expect_error(sfTDist(param = "abc"),
               info = "Checking for incorrect variable type")
})

######### analysis schedule setting tests #########
test_that("mSetSchedule handles scalar (0/1) inputs correctly", {
  k <- 3

  # Both arms: pattern 0 -> constant max clusters at each analysis
  mt0 <- mSetSchedule(c(0, 0), k = k)
  expect_equal(dim(mt0), c(2L, k))
  expect_equal(mt0[1, ], rep(1, k))
  expect_equal(mt0[2, ], rep(1, k))

  # First arm constant, second arm linear accrual
  mt01 <- mSetSchedule(c(0, 1), k = k)
  expect_equal(dim(mt01), c(2L, k))
  expect_equal(mt01[1, ], rep(1, k))
  expect_equal(mt01[2, ], seq_len(k) / k)

  # Length-1 input recycled to both arms
  mt1 <- mSetSchedule(1, k = k)
  expect_equal(dim(mt1), c(2L, k))
  expect_equal(mt1[1, ], seq_len(k) / k)
  expect_equal(mt1[2, ], seq_len(k) / k)
})

test_that("nSetSchedule handles scalar (0/1) inputs correctly", {
  k <- 4

  # Both arms: pattern 1 -> linear accrual in cluster size
  nt11 <- nSetSchedule(c(1, 1), k = k)
  expect_equal(dim(nt11), c(2L, k))
  expect_equal(nt11[1, ], seq_len(k) / k)
  expect_equal(nt11[2, ], seq_len(k) / k)

  # First arm linear, second arm constant
  nt10 <- nSetSchedule(c(1, 0), k = k)
  expect_equal(dim(nt10), c(2L, k))
  expect_equal(nt10[1, ], seq_len(k) / k)
  expect_equal(nt10[2, ], rep(1, k))
})

test_that("mSetSchedule and nSetSchedule expand 2 x (k-1) matrices and enforce final timing = 1", {
  k <- 3

  # Valid 2 x (k-1) matrix: last column should be appended as 1s
  m_in <- matrix(c(0.3, 0.3,
                   0.7, 0.7), nrow = 2)  # 2 x 2, k = 3
  mt <- mSetSchedule(m_in, k = k)
  expect_equal(dim(mt), c(2L, k))
  expect_equal(mt[, 1:2], m_in)
  expect_equal(mt[, 3], rep(1, 2))

  n_in <- matrix(c(0.25, 0.25,
                   0.6,  0.6), nrow = 2) # 2 x 2, k = 3
  nt <- nSetSchedule(n_in, k = k)
  expect_equal(dim(nt), c(2L, k))
  expect_equal(nt[, 1:2], n_in)
  expect_equal(nt[, 3], rep(1, 2))

  # Invalid: full 2 x k matrix with final timing != 1 should error
  m_bad_final <- matrix(c(0.5, 0.5,
                          0.9, 0.9), nrow = 2) # k = 2
  expect_error(
    mSetSchedule(m_bad_final, k = 2),
    info = "final analysis is input, it must be 1"
  )

  n_bad_final <- matrix(c(0.4, 0.4,
                          0.8, 0.8), nrow = 2)
  expect_error(
    nSetSchedule(n_bad_final, k = 2),
    info = "final analysis is input, it must be 1"
  )
})

test_that("mSetSchedule and nSetSchedule enforce 2-row matrices and strictly increasing timings", {
  # Wrong number of rows
  m_bad_rows <- matrix(0.5, nrow = 4, ncol = 1)
  expect_error(
    mSetSchedule(m_bad_rows, k = 3),
    info = "if analysis timing for final analysis is input, it must be 1"
  )

  n_bad_rows <- matrix(0.5, nrow = 4, ncol = 1)
  expect_error(
    nSetSchedule(n_bad_rows, k = 3),
    info = "if analysis timing for final analysis is input, it must be 1"
  )

  # Not strictly increasing: flat timing should fail
  m_flat <- matrix(c(0.5, 0.5,
                     0.5, 0.5), nrow = 2)  # 2 x 2, k = 2
  expect_error(
    mSetSchedule(m_flat, k = 2),
    info = "Must be increasing strictly"
  )

  n_flat <- matrix(c(0.3, 0.3,
                     0.3, 0.3), nrow = 2)
  expect_error(
    nSetSchedule(n_flat, k = 2),
    info = "Must be increasing strictly"
  )
})

######## genClusterSizes tests #########

test_that("genClusterSizes with zero CV returns constant sizes", {
  sizes <- genClusterSizes(
    m     = 3,
    n     = 10,
    n_cv  = 0,
    n_min = rep(1, 3),
    n_max = rep(20, 3)
  )

  expect_equal(length(sizes), 3)
  expect_true(all(sizes == 10))
})

test_that("genClusterSizes errors for invalid negative binomial dispersion", {
  # n_cv^2 <= 1 / n should error
  expect_error(
    genClusterSizes(
      m     = 3,
      n     = 10,
      n_cv  = sqrt(1 / 10),  # exactly on the boundary
      n_min = rep(1, 3),
      n_max = rep(20, 3)
    ),
    info = "n_cv^2 must be greater"
  )
})

######## genContCRT tests #########

test_that("genContCRT generates continuous CRT data properly", {
  df <- genContCRT(
    m        = c(2, 2),
    m_alloc  = c(0.5, 0.5),
    n        = c(5, 5),
    n_cv     = c(0, 0),
    mu_vec   = c(0, 1),
    sigma_vec = c(1, 1),
    rho      = c(0.01, 0.01)
  )

  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("arm", "cluster", "individual", "response"))

  expect_true(all(df$arm %in% 0:1))
  expect_true(all(df$cluster >= 1))
  expect_true(is.numeric(df$response))
  expect_true(nrow(df) > 0)
})

######### genBinCRT tests #########

test_that("genBinCRT generates binary CRT data properly", {
  df <- genBinCRT(
    m        = c(2, 2),
    m_alloc  = c(0.5, 0.5),
    n        = c(5, 5),
    n_cv     = c(0, 0),
    p_vec    = c(0.3, 0.5),
    rho      = c(0.01, 0.01)
  )

  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("arm", "cluster", "individual", "response"))

  expect_true(all(df$arm %in% 0:1))
  expect_true(all(df$cluster >= 1))
  expect_true(all(df$response %in% c(0, 1)))
  expect_true(nrow(df) > 0)
})

######### gsSimContCRT tests #########

test_that("gsSimContCRT runs with known variance/ICC and returns expected structure", {
  # Simple continuous data set compatible with m_max / n_max
  m_max <- c(3, 3)
  n_max <- c(8, 8)
  data_cont <- genContCRT(
    m        = m_max,
    n        = n_max,
    n_cv     = c(0, 0),
    mu_vec   = c(0, 0.5),
    sigma_vec = c(1, 1),
    rho      = c(0.01, 0.01)
  )

  sim <- gsSimContCRT(
    k           = 3,
    data        = data_cont,
    test_type   = 1,
    test_sides  = 1,
    recruit_type = 1,
    stat_type   = 1,                 # Z_known
    delta       = 0.5,
    sigma_vec   = c(1, 1),
    rho         = c(0.01, 0.01),
    alpha       = 0.05,
    beta        = 0.1,
    i_max       = 1,                  # any positive value; just affects ifrac
    m_max       = m_max,
    n_max       = n_max,
    n_cv        = c(0, 0),
    schedule_m  = NULL,               # let function choose equal fractions
    schedule_n  = NULL,
    alpha_sf    = sfLDOF,
    beta_sf     = sfLDOF
  )

  expect_type(sim, "list")
  expect_true(all(c("reject", "t_i", "m_i", "n_i", "total_i",
                    "i_frac") %in% names(sim)))

  expect_true(is.logical(sim$reject))
  expect_true(length(sim$reject) == 1L)
  expect_true(is.numeric(sim$t_i))
  expect_true(sim$t_i >= 1)
  expect_true(is.numeric(sim$m_i))
  expect_true(is.numeric(sim$n_i))
  expect_true(is.numeric(sim$total_i))
})

######### gsSimBinCRT tests #########

test_that("gsSimBinCRT runs with known variance/ICC and returns expected structure", {
  set.seed(123)

  m_max <- c(3L, 3L)
  n_max <- c(8L, 8L)
  data_bin <- genBinCRT(
    m        = m_max,
    n        = n_max,
    n_cv     = c(0, 0),
    p_vec    = c(0.3, 0.5),
    rho      = c(0.01, 0.01)
  )

  sim <- gsSimBinCRT(
    k           = 3L,
    data        = data_bin,
    test_type   = 1L,
    test_sides  = 1L,
    recruit_type = 1L,
    stat_type   = 1L,                 # Z_known
    delta       = 0.2,
    p_vec       = c(0.3, 0.5),
    rho         = c(0.01, 0.01),
    alpha       = 0.05,
    beta        = 0.1,
    i_max       = 1,
    m_max       = m_max,
    n_max       = n_max,
    n_cv        = c(0, 0),
    schedule_m  = NULL,
    schedule_n  = NULL,
    alpha_sf    = sfLDOF,
    beta_sf     = sfLDOF
  )

  expect_type(sim, "list")
  expect_true(all(c("reject", "t_i", "m_i", "n_i", "total_i",
                    "i_frac") %in% names(sim)))

  expect_true(is.logical(sim$reject))
  expect_true(length(sim$reject) == 1L)
  expect_true(is.numeric(sim$t_i))
  expect_true(sim$t_i >= 1)
  expect_true(is.numeric(sim$m_i))
  expect_true(is.numeric(sim$n_i))
  expect_true(is.numeric(sim$total_i))
})