source('../gsDesign_independent_code.R')
# --------------------------------------------
# Test plot.gsDesignCRT function
## save_gg_plot() is used for storing plots created using ggplot2 package,
## whereas save_gr_plot() is used with plots created using graphics package.
#----------------------------------------------

test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 1 and base R plotting package set to TRUE", code = {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 1, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plotgsds1.png")
})

test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 1 and base R plotting package set to FALSE.", code = {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gg_plot(plot.gsDesignCRT(x, plottype = 1, base = FALSE))
  expect_snapshot_file(save_plot_obj, "plotgsds2.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype = power and base R plotting package set to FALSE.", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gg_plot(plot.gsDesignCRT(x, plottype = "power", base = FALSE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_3.png")
})

test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 2 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 2, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_4.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 3 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 3, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_5.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 4 and base R plotting package set to FALSE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gg_plot(plot.gsDesignCRT(x, plottype = 4, base = FALSE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_6.png")
})

test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 4 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 4, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_7.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 5 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 5, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_8.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 6 and base R plotting package set to FALSE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gg_plot(plot.gsDesignCRT(x, plottype = 6, base = FALSE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_9.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 6 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 6, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_11.png")
})

test_that("Test plot.gsDesignCRT graphs are correctly rendered for plottype 7 and base R plotting package set to TRUE", {
  x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(x, plottype = 7, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsds_12.png")
})


test_that("Test plot.gsDesignCRT graphs are correctly rendered for gsSurvR and base R plotting package set to TRUE", {
  z <- gsSurv(
    k = 3, sfl = sfPower, sflpar = .5,
    lambdaC = log(2) / 6, hr = .5, eta = log(2) / 40,
    gamma = 1, T = 36, minfup = 12
  )

  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(z, plottype = 2, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsdsa_3.png")
})


test_that("Test plotgsDesign graphs are correctly rendered for gSurv and base R plotting package set to False", {
  z <- gsSurv(
    k = 3, sfl = sfPower, sflpar = .5,
    lambdaC = log(2) / 6, hr = .5, eta = log(2) / 40,
    gamma = 1, T = 36, minfup = 12
  )
  save_plot_obj <- save_gg_plot(plot.gsDesignCRT(z, plottype = 2, base = FALSE))
  local_edition(3)
  expect_snapshot_file(save_plot_obj, "plot_plotgsdsa_4.png")
})

test_that("Test plot.gsDesignCRT graphs are correctly rendered for testtype 1 and plottype 2", {
  y <- gsDesignCRT(k = 5, test.type = 1, n.fix = 1)
  local_edition(3)
  save_plot_obj <- save_gr_plot(plot.gsDesignCRT(y, plottype = 2, base = TRUE))
  expect_snapshot_file(save_plot_obj, "plot_plotgsdsa_8.png")
})
