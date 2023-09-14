#Check LaTeX output for gsDesignCRT
#in xprint argument "comment = FALSE" to handle timestamp in snapshot.
#-------------------------------------------------------------------------------

test_that(
  desc = "test: Checking LaTeX output for gsDesignCRT",
  code = {
    n.fix <- nBinomial(p1 = .3, p2 = .15, scale = "OR")
    xOR <- gsDesignCRT(k = 2, n.fix = n.fix, delta1 = log(.15 / .3 / .85 * .7),
                    endpoint = "Binomial")
    local_edition(3)  #use the 3rd edition of the testthat package
    
    expect_snapshot_output(x = xprint(
      xtable::xtable(gsBoundSummary(xOR, deltaname = "OR", logdelta = TRUE)),
      comment = FALSE
    ))
  }
)
