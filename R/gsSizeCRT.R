# gsSizeCRT roxy [sinew] ----
#' @title Compute sample size for group sequential cluster-randomized trial
#' @description  \code{gsSizeCRT()} is used to trial size required for a
#' cluster-randomized group sequential design.
#'
#' @param t Number of analyses planned, including interim and final.
#' @param test_type \code{1=}one-sided test with early stopping for efficacy
#' \cr \code{2=}two-sided test with early stopping for efficacy
#' \cr \code{3=}one-sided test with efficacy or binding futility stopping
#' \cr \code{4=}two-sided test with efficacy or binding futility stopping
#' \cr \code{5=}one-sided test with efficacy or non-binding futility stopping
#' \cr \code{6=}two-sided test with efficacy or non-binding futility stopping
#' @param outcome_type \code{1=}continuous difference of means
#' \cr \code{2=}binary difference of proportions \cr \code{3=}binary odds ratio
#' @param size_type \code{1=}clusters \cr \code{2=}individuals per cluster
#' @param alpha Type I error, always one-sided. Default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#' @param sigma_vec Standard deviations for groups 1 and 2 (continuous case)
#' @param p_vec Probabilities of event for groups 1 and 2 (binary case)
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param delta Effect size for theta under alternative hypothesis.
#' @param m_fix Number of clusters; used to find maximum group sequential
#' size of each cluster
#' @param n_fix Mean size of each cluster; used to find maximum group sequential
#' number of clusters
#' @param timing Sets relative timing of interim analyses. Default of 1
#' produces equally spaced analyses.  Otherwise, this is a vector of length
#' \code{t} or \code{t-1}.  The values should satisfy \code{0 < timing[1] <
#' timing[2] < ... < timing[t-1] < timing[t]=1}.
#' @param alpha_sf A spending function or a character string indicating a
#' boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds, \dQuote{OF} for
#' O'Brien-Fleming bounds and \dQuote{Pocock} for Pocock bounds).  For
#' one-sided and symmetric two-sided testing is used to completely specify
#' spending (\code{test.type=1, 2}), \code{sfu}.  The default value is
#' \code{sfHSD} which is a Hwang-Shih-DeCani spending function.  See details,
#' \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param alpha_sfpar Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with the default
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions.  The parameter \code{sfupar} specifies any parameters needed for
#' the spending function specified by \code{sfu}; this will be ignored for
#' spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf Specifies the spending function for lower boundary crossing
#' probabilities when asymmetric, two-sided testing is performed
#' (\code{test.type = 3}, \code{4}, \code{5}, or \code{6}).  Unlike the upper
#' bound, only spending functions are used to specify the lower bound.  The
#' default value is \code{sfHSD} which is a Hwang-Shih-DeCani spending
#' function.  The parameter \code{sfl} is ignored for one-sided testing
#' (\code{test.type=1}) or symmetric 2-sided testing (\code{test.type=2}).  See
#' details, spending functions, manual and examples.
#' @param beta_sfpar Real value, default is \eqn{-2}, which, with the default
#' Hwang-Shih-DeCani spending function, specifies a less conservative spending
#' rate than the default for the upper bound.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#'
#' @return Object with maximum information sample size and expected sample size
#' under null and alternate hypotheses
#' @export
#' 
#' @useDynLib gsDesignCRT gsbounds1
#' @useDynLib gsDesignCRT gsbounds2
#' 
#' @name gsSizeCRT
# gsSizeCRT function [sinew] ----
gsSizeCRT <- function(t = 3, test_type = 1, outcome_type = 1, size_type = 1,
                      stat_type = 1, delta = 1, sigma_vec = NULL, p_vec = NULL,
                      rho = 0, alpha = 0.05, beta = 0.1, m_fix = 1, n_fix = 1,
                      timing = 1, alpha_sf = sfLDOF, alpha_sfpar = -4,
                      beta_sf = sfLDOF, beta_sfpar = -4,
                      tol = 0.000001, r = 42) {
  x <- list(t = t, test_type = test_type, outcome_type = outcome_type,
            size_type = size_type, stat_type = stat_type, delta = delta,
            sigma_vec = sigma_vec, p_vec = p_vec, rho = rho, alpha = alpha,
            beta = beta, m_fix = m_fix, n_fix = n_fix, timing = timing,
            tol = tol, r = r, i = 0, max_i = 0, i_frac = 0, e_i = 0,
            m = m_fix, max_m = m_fix, e_m = c(m_fix, m_fix),
            n = n_fix, max_n = n_fix, e_n = c(n_fix, n_fix),
            total = m_fix * n_fix, max_total = m_fix * n_fix,
            e_total = m_fix * n_fix, sufficient = 1, m_frac = 0, n_frac = 0)
  class(x) <- "gsDesignCRT"

  # Compute maximum and expected information
  x_info <- gsMaxInfoCRT(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar)
  x$i <- x_info$i
  x$max_i <- x_info$i[x$t]
  x$i_frac <- x$i / x$max_i
  x$e_i <- x_info$e_i
  x$lb <- x_info$lb
  x$ub <- x_info$ub
  x$falsepos <- x_info$falsepos

  # Convert maximum and expected information to sample sizes
  if (size_type == 1) {
    if (outcome_type == 1) {
      x$m <- mContDiff(x$i, x$n_fix, x$sigma_vec, x$rho)
      x$max_m <- mContDiff(x$max_i, x$n_fix, x$sigma_vec, x$rho)
      x$e_m <- mContDiff(x$e_i, x$n_fix, x$sigma_vec, x$rho)
    } else if (outcome_type == 2) {
      x$m <- mPropDiff(x$i, x$n_fix, x$p_vec, x$rho)
      x$max_m <- mPropDiff(x$max_i, x$n_fix, x$p_vec, x$rho)
      x$e_m <- mPropDiff(x$e_i, x$n_fix, x$p_vec, x$rho)
    } else {
      stop("invalid outcome type")
    }
  } else if (size_type == 2) {
    if (outcome_type == 1) {
      x$sufficient <- checkSuff(outcome_type = x$outcome_type,
                                m = x$m_fix,
                                max_i = x$max_i,
                                sigma_vec = x$sigma_vec,
                                rho = x$rho,
                                alpha = x$alpha,
                                beta = x$beta)

      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
      }
      x$n <- nContDiff(x$i, x$m_fix, x$sigma_vec, x$rho)
      x$max_n <- nContDiff(x$max_i, x$m_fix, x$sigma_vec, x$rho)
      x$e_n <- nContDiff(x$e_i, x$m_fix, x$sigma_vec, x$rho)
    } else if (outcome_type == 2) {
      x$sufficient <- checkSuff(outcome_type = x$outcome_type,
                                m = x$m_fix,
                                max_i = x$max_i,
                                p_vec = x$p_vec,
                                rho = x$rho,
                                alpha = x$alpha,
                                beta = x$beta)
      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
      }
      x$n <- nPropDiff(x$i, x$m_fix, x$p_vec, x$rho)
      x$max_n <- nPropDiff(x$max_i, x$m_fix, x$p_vec, x$rho)
      x$e_n <- nPropDiff(x$e_i, x$m_fix, x$p_vec, x$rho)
    } else {
      stop("invalid outcome type")
    }
  } else {
    stop("invalid variable")
  }

  # Compute sample sizes under information fractions
  x$m_frac <- sizeInfoFrac(1, x$i_frac, x$max_m, x$max_n, x$rho)
  x$n_frac <- sizeInfoFrac(2, x$i_frac, x$max_m, x$max_n, x$rho)

  # Compute total number per arm
  x$total <- x$m * x$n
  x$max_total <- x$max_m * x$max_n
  x$e_total <- x$e_m * x$e_n
  return(x)
}

# gsMaxInfoCRT function [sinew] ----
gsMaxInfoCRT <- function(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar) {
  # Check inputs
  checkScalar(x$t, "integer", c(1, Inf))
  checkScalar(x$test_type, "integer", c(1, 6))
  checkScalar(x$outcome_type, "integer", c(1, 2))
  checkScalar(x$size_type, "integer", c(1, 2))
  checkScalar(x$stat_type, "integer", c(1, 3))
  checkScalar(x$alpha, "numeric", 0:1, c(FALSE, FALSE))
  if (x$alpha > 0.5) {
    checkScalar(x$alpha, "numeric", c(0, 0.5), c(FALSE, TRUE))
  }
  checkScalar(x$beta, "numeric", c(0, 1 - x$alpha), c(FALSE, FALSE))
  checkScalar(x$rho, "numeric", c(0, 1))
  checkScalar(x$tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  checkScalar(x$r, "integer", c(1, 80))

  # Specify timing depending on provided information
  if (length(x$timing) < 1 ||
        (length(x$timing) == 1 &&
           (x$t > 2 || (x$t == 2 && (x$timing[1] <= 0 || x$timing[1] >= 1))))) {
    x$timing <- seq(x$t) / x$t
  } else if (x$timing == 1 && x$t == 1) {
    x$timing <- 1
  } else if (length(x$timing) == x$t - 1 || length(x$timing) == x$t) {
    if (length(x$timing) == x$t - 1) {
      x$timing <- c(x$timing, 1)
    } else if (x$timing[x$t] != 1) {
      stop("if analysis timing for final analysis is input, it must be 1")
    }
    if (min(x$timing - c(0, x$timing[1:(x$t - 1)])) <= 0) {
      stop("input timing of interim analyses must be increasing strictly
           between 0 and 1")
    }
  } else {
    stop("value input for timing must be length 1, k-1 or k")
  }

  # Get error spending at information levels
  if (is.character(alpha_sf) &&
        !is.element(alpha_sf, c("OF", "Pocock", "WT"))) {
    stop("Character specification of upper spending may only be WT, OF or
         Pocock")
  } else if (!is.function(alpha_sf)) {
    stop("Upper spending function mis-specified")
  }

  if (x$test_type == 1) {
    # Partition error probabilities based on error spending
    upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$t - 1])

    # Determine boundary values to achieve Type I error = alpha
    lower_bound <- rep(-30, x$t)
    bounds <- gsUpper1(0, x$timing, lower_bound, falsepos, x$tol, x$r)
    upper_bound <- bounds$b

    # Determine maximum information based on computed boundary values
    i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
    x$i <- stats::uniroot(gsbetadiff, lower = i0, upper = 10 * i0,
                          extendInt = "yes", theta = x$delta, beta = x$beta,
                          time = x$timing, a = lower_bound, b = upper_bound,
                          tol = x$tol, r = x$r)$root * x$timing

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprob1(theta, x$i, lower_bound, upper_bound, x$r)
    x$e_i <- y$eI
  } else if (x$test_type == 2) {
    # Partition error probabilities based on error spending
    upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$t - 1])

    # Determine boundary values to achieve Type I error = alpha
    lower_bound <- rep(0, x$t)
    bounds <- gsUpper2(0, x$timing, lower_bound, falsepos, x$tol, x$r)
    lower_bound <- bounds$a
    upper_bound <- bounds$b

    # Determine maximum information based on computed boundary values
    i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
    x$i <- stats::uniroot(gsbetadiff2, lower = i0, upper = 10 * i0,
                          extendInt = "yes", theta = x$delta, beta = x$beta,
                          time = x$timing, a = lower_bound, b = upper_bound,
                          tol = x$tol, r = x$r)$root * x$timing

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprob2(theta, x$i, lower_bound, upper_bound, x$r)
    x$e_i <- y$eI
  } else if (x$test_type == 3 || x$test_type == 5) {
    # Partition error probabilities based on error spending
    upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$t - 1])

    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$t - 1])

    # Determine maximum information that achieves a_T = b_T with error spending
    if (x$test_type == 3) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
    i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
    x$i <- stats::uniroot(gsbounddiff, lower = i0, upper = 10 * i0,
                          extendInt = "yes", theta = x$delta,
                          falseneg = falseneg, falsepos = falsepos,
                          time = x$timing, binding = binding,
                          tol = x$tol, r = x$r)$root * x$timing

    # Determine boundary values to achieve Type I error = alpha
    bounds <- gsBounds1(x$delta, x$i, falseneg, falsepos, binding, x$tol, x$r)
    lower_bound <- bounds$a
    upper_bound <- bounds$b

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprob1(theta, x$i, lower_bound, upper_bound, x$r)
    x$e_i <- y$eI
  } else if (x$test_type == 4 || x$test_type == 6) {
    # Partition error probabilities based on error spending
    upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$t - 1])

    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$t - 1])

    # Determine maximum information that achieves a_T = b_T with error spending
    if (x$test_type == 4) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
    i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
    x$i <- stats::uniroot(gsbounddiff2, lower = i0, upper = 10 * i0,
                          extendInt = "yes", theta = x$delta,
                          falseneg = falseneg, falsepos = falsepos,
                          time = x$timing, binding = binding,
                          tol = x$tol, r = x$r)$root * x$timing

    # Determine boundary values to achieve Type I error = alpha
    bounds <- gsBounds2(x$delta, x$i, falseneg, falsepos, binding, x$tol, x$r)
    lower_bound <- bounds$a
    upper_bound <- bounds$b

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprob2(theta, x$i, lower_bound, upper_bound, x$r)
    x$e_i <- y$eI
  } else {
    stop("invalid test type")
  }

  return(x)
}

# mContDiff roxy [sinew] ----
#' @title Convert information to number of clusters for difference of continuous
#' outcomes
#'
#' @param i Fisher information.
#' @param n Mean cluster size.
#' @param sigma_vec Standard deviations for groups 1 and 2 (continuous case).
#' @param rho Intraclass correlation coefficient.
#'
#' @return Number of clusters corresponding to inputs.
#' @export
#' @name mContDiff
# mContDiff function [sinew] ----
mContDiff <- function(i, n, sigma_vec, rho) {
  d_eff <- (1 + (n - 1) * rho) / n
  m <- i * (sigma_vec[1]^2 + sigma_vec[2]^2) * d_eff
  return(m)
}

# mPropDiff roxy [sinew] ----
#' @title Convert information to number of clusters for difference of
#' proportions
#'
#' @param i Fisher information.
#' @param n Mean cluster size.
#' @param p_vec Probability of success for groups 1 and 2.
#' @param rho Intraclass correlation coefficient.
#'
#' @return Number of clusters corresponding to inputs.
#' @export
#' @name mPropDiff
# mPropDiff function [sinew] ----
mPropDiff <- function(i, n, p_vec, rho) {
  p <- mean(p_vec)
  d_eff <- (1 + (n - 1) * rho) / n
  m <- i * (2 * p * (1 - p)) * d_eff
  return(m)
}

# nContDiff roxy [sinew] ----
#' @title Convert information to average cluster size for difference of
#' continuous outcomes
#'
#' @param i Fisher information.
#' @param m Number of clusters per arm.
#' @param sigma_vec Standard deviations for groups 1 and 2 (continuous case).
#' @param rho Intraclass correlation coefficient.
#'
#' @return Average cluster size corresponding to inputs.
#' @export
#' @name nContDiff
# nContDiff function [sinew] ----
nContDiff <- function(i, m, sigma_vec, rho) {
  num <- i * (sigma_vec[1]^2 + sigma_vec[2]^2) * (1 - rho) / m
  denom <- 1 - (i * (sigma_vec[1]^2 + sigma_vec[2]^2) * rho / m)
  n <- num / denom
  return(n)
}

# nPropDiff roxy [sinew] ----
#' @title Convert information to average cluster size for difference of
#' proportions
#'
#' @param i Fisher information.
#' @param m Number of clusters per arm.
#' @param p_vec Probability of success for groups 1 and 2.
#' @param rho Intraclass correlation coefficient.
#'
#' @return Average cluster size corresponding to inputs.
#' @export
#' @name nPropDiff
# nPropDiff function [sinew] ----
nPropDiff <- function(i, m, p_vec, rho) {
  p <- mean(p_vec)
  num <- i * (2 * p * (1 - p)) * (1 - rho) / m
  denom <- 1 - (i * (2 * p * (1 - p)) * rho / m)
  n <- num / denom
  return(n)
}

# checkSuff roxy [sinew] ----
#' @title Check whether number of clusters sufficient to achieve specified Type
#' I error and power
#'
#' @param test_type \code{1=}one-sided test with early stopping for efficacy
#' \cr \code{2=}two-sided symmetric with early stopping for efficacy
#' \cr \code{3=}one-sided with early stopping for efficacy or futility
#' \cr \code{4=}two-sided with early stopping for efficacy or futility
#' @param outcome_type \code{1=}continuous difference of means
#' \cr \code{2=}binary difference of proportions \cr \code{3=}binary odds ratio
#' @param m Number of clusters per arm.
#' @param sigma_vec Standard deviations for groups 1 and 2 (continuous case).
#' @param p_vec Probability of success for groups 1 and 2.
#' @param rho Intraclass correlation coefficient.
#' @param delta Effect size for theta under alternative hypothesis.
#' @param alpha Type I error, always one-sided. Default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#'
#' @return Binary variable corresponding to whether m clusters is sufficient to
#' achieve Type I error and power to detect an effect size of delta
#'
#' @export
#' @name checkSuff
# checkSuff function [sinew] ----
checkSuff <- function(outcome_type, m, max_i, sigma_vec = c(1, 1),
                      p_vec = c(0, 0), rho, alpha = 0.05, beta = 0.1) {
  # Compute individually randomized sample size for no interim looks
  if (outcome_type == 1) {
    ni <- max_i * (sigma_vec[1]^2 + sigma_vec[2]^2)
  } else if (outcome_type == 2) {
    p <- mean(p_vec)
    ni <- max_i * (2 * p * (1 - p))
  } else {
    stop("invalid outcome type")
  }

  # Conduct feasibility check
  check <- 0
  if (m > ni * rho) {
    check <- 1
  }
  return(check)
}

# sizeInfoFrac roxy [sinew] ----
#' @title Convert information fraction to sample size
#'
#' @param size_type \code{1=}clusters \cr \code{2=}individuals per cluster
#' @param ifrac information fraction I_{t} / I_{max}
#' @param n Mean cluster size.
#' @param p_vec Probability of success for groups 1 and 2.
#' @param rho Intraclass correlation coefficient.
#'
#' @return Sample size corresponding to inputs.
#' @export
#' @name sizeInfoFrac
# sizeInfoFrac function [sinew] ----
sizeInfoFrac <- function(size_type, ifrac, m_max = 1, n_max = 1, rho) {
  if (size_type == 1) {
    size <- ifrac * m_max
  } else if (size_type == 2) {
    size_num <- ifrac * (n_max / (1 + (n_max - 1) * rho)) * (1 - rho)
    size_denom <- 1 - ifrac * (n_max / (1 + (n_max - 1) * rho)) * rho
    size <- size_num / size_denom
  } else {
    stop("invalid outcome type")
  }
  return(size)
}