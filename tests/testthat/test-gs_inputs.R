testthat::context("gs inputs")

testthat::test_that("test.gsDesignCRT.alpha", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(alpha = "abc", test.type = 1), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(alpha = 0, test.type = 1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(alpha = 1, test.type = 1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(alpha = 0.51, test.type = 2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(alpha = rep(0.5, 2), test.type = 1),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsDesignCRT.astar", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(astar = "abc", test.type = 5), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(astar = 0.51, alpha = 0.5, test.type = 5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::gsDesignCRT(astar = 1, alpha = 0, test.type = 5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::gsDesignCRT(astar = -1, test.type = 6), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(astar = rep(0.1, 2), alpha = 0.5, test.type = 5),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsDesignCRT.beta", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(beta = "abc", test.type = 3), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(beta = 0.5, alpha = 0.5, test.type = 3),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::gsDesignCRT(beta = 1, alpha = 0, test.type = 3),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsDesignCRT::gsDesignCRT(beta = 0, test.type = 3), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(beta = rep(0.1, 2), alpha = 0.5, test.type = 3),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsDesignCRT.delta", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(delta = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(delta = -1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(delta = rep(0.1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.k", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = 1.2), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = 24, test.type = 4), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = seq(2)), info = "Checking for incorrect variable length")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(k = 3, sfu = sfpoints, sfupar = c(
    0.05,
    0.1, 0.15, 0.2, 1
  )), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.maxn.I", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(maxn.I = "abc"), info = "Checking for incorrect variable type")
})

testthat::test_that("test.gsDesignCRT.n.fix", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(n.fix = "abc", delta = 0), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(n.fix = -1, delta = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(n.fix = rep(2, 2), delta = 0), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.n.I", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(n.I = "abc"), info = "Checking for incorrect variable type")
})

testthat::test_that("test.gsDesignCRT.r", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(r = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(r = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(r = 81), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(r = rep(1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.sfl", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfl = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfl = rep(sfHSD, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.sflpar", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sflpar = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sflpar = rep(-2, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.sfu", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfu = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfu = rep(sfHSD, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.sfupar", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfupar = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(sfupar = rep(-4, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.test.type", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 1.2), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 7), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = seq(2)), info = "Checking for incorrect variable length")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 3, sfu = "WT"), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 4, sfu = "WT"), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 5, sfu = "WT"), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(test.type = 6, sfu = "WT"), info = "Checking for out-of-range variable value")
})

testthat::test_that("test.gsDesignCRT.timing", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = "abc", k = 1), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = -1, k = 1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = 2, k = 1), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = c(0.1, 1.1), k = 2), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = c(0.5, 0.1), k = 2), info = "NA")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(timing = c(0.1, 0.5, 1), k = 2), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsDesignCRT.tol", {
  testthat::expect_error(gsDesignCRT::gsDesignCRT(tol = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(tol = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(tol = 0.10000001), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsDesignCRT(tol = rep(0.1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsBoundCP.r", {
  testthat::expect_error(gsDesignCRT::gsBoundCP(r = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsBoundCP(r = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsBoundCP(r = 81), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsBoundCP(r = rep(1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsBoundCP.theta", {
  testthat::expect_error(gsDesignCRT::gsBoundCP(theta = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsBoundCP(theta = rep(1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsBoundCP.x", {
  testthat::expect_error(gsDesignCRT::gsBoundCP(x = "abc"), info = "Checking for incorrect variable type")
})

testthat::test_that("test.gsCP.r", {
  testthat::expect_error(gsDesignCRT::gsCP(r = "abc"), info = "Checking for incorrect variable type")
  testthat::expect_error(gsDesignCRT::gsCP(r = 0), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsCP(r = 81), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsDesignCRT::gsCP(r = rep(1, 2)), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsCP.x", {
  testthat::expect_error(gsDesignCRT::gsCP(x = "abc"), info = "Checking for incorrect variable type")
})

testthat::test_that("test.gsbound1.a", {
  testthat::expect_error(gsbound1(a = "abc", theta = 1, I = 1, probhi = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound1(
    a = rep(0.5, 2), theta = 1, I = 1,
    probhi = 0.5
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsbound1.I", {
  testthat::expect_error(gsbound1(I = "abc", theta = 1, a = 0, probhi = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound1(I = 0, theta = 1, a = 0, probhi = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound1(I = rep(1, 2), theta = 1, a = 0, probhi = 0.5),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsbound1.probhi", {
  testthat::expect_error(gsbound1(probhi = "abc", I = 1, a = 0, theta = 1),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound1(probhi = 0, I = 1, a = 0, theta = 1),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound1(probhi = 1, I = 1, a = 0, theta = 1),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound1(
    probhi = rep(0.5, 2), I = 1, a = 0,
    theta = 1
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsbound1.r", {
  testthat::expect_error(gsbound1(
    r = "abc", I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsbound1(
    r = 0, I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsbound1(
    r = 81, I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for out-of-range variable value")
  testthat::expect_error(gsbound1(
    r = rep(1, 2), I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsbound1.theta", {
  testthat::expect_error(gsbound1(theta = "abc", I = 1, a = 0, probhi = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound1(theta = rep(1, 2), I = 1, a = 0, probhi = 0.5),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsbound1.tol", {
  testthat::expect_error(gsbound1(
    tol = "abc", I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for incorrect variable type")
  testthat::expect_error(gsbound1(
    tol = 0, I = 1, a = 0, probhi = 0.5,
    theta = 1
  ), info = "Checking for out-of-range variable value")
})

testthat::test_that("test.gsbound.falsepos", {
  testthat::expect_error(gsbound(falsepos = "abc", I = 1, trueneg = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound(falsepos = 0, I = 1, trueneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(falsepos = 1, I = 1, trueneg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(falsepos = rep(0.5, 2), I = 1, trueneg = 0.5),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsbound.I", {
  testthat::expect_error(gsbound(I = "abc", trueneg = 0.5, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound(I = 0, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(I = rep(1, 2), trueneg = 0.5, falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
})

testthat::test_that("test.gsbound.r", {
  testthat::expect_error(gsbound(r = "abc", I = 1, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound(r = 0, I = 1, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(r = 81, I = 1, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(
    r = rep(1, 2), I = 1, trueneg = 0.5,
    falsepos = 0.5
  ), info = "Checking for incorrect variable length")
})

testthat::test_that("test.gsbound.tol", {
  testthat::expect_error(gsbound(tol = "abc", I = 1, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound(tol = 0, I = 1, trueneg = 0.5, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

testthat::test_that("test.gsbound.trueneg", {
  testthat::expect_error(gsbound(trueneg = "abc", I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable type"
  )
  testthat::expect_error(gsbound(trueneg = 0, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(trueneg = 1, I = 1, falsepos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  testthat::expect_error(gsbound(trueneg = rep(0.5, 2), I = 1, falsepos = 0.5),
    info = "Checking for incorrect variable length"
  )
})
