# gsDesignCRT roxy [sinew] ----
#' @title Compute stopping boundaries, maximum sample size, and expected sample
#'  sizes for a group sequential cluster randomized trial.
#'
#' @description  \code{gsDesignCRT()} is used to determine the maximum sample
#'  size needed for a specified parallel group sequential cluster randomized
#'  trial to detect a clinically meaningful effect size with some Type I error
#'  rate and power. Code adapted from gsDesign package.
#'
#' @param k Number of analyses planned, including interim and final.
#' @param outcome_type \code{1=}continuous difference of means
#' \cr \code{2=}binary difference of proportions
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param size_type \code{1=}clusters per arm \cr \code{2=}cluster size
#' @param recruit_type \code{1=}recruit clusters with fixed sizes
#' \cr \code{2=}recruit individuals into fixed number of clusters
#' @param delta Effect size for theta under alternative hypothesis. Must be > 0.
#' @param sigma_vec Standard deviations for control and treatment groups
#'  (continuous case).
#' @param p_vec Probabilities of event for control and treatment groups
#'  (binary case).
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Type I error, default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#' @param m_fix Number of clusters; used to find maximum size of each cluster.
#' @param n_fix Mean size of each cluster; used to find maximum number of
#'  clusters per arm.
#' @param timing Sets relative timing of interim analyses. Default of 1
#'  produces equally spaced analyses.  Otherwise, this is a vector of length
#'  \code{k} or \code{k-1}.  The values should satisfy \code{0 < timing[1] <
#'  timing[2] < ... < timing[k-1] < timing[k]=1}.
#' @param alpha_sf A spending function or a character string indicating an
#'  upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#'  \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#'  bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#'  O'Brien-Fleming spending function. See details,
#'  \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param alpha_sfpar Real value, default is \eqn{-4} which is an
#'  O'Brien-Fleming-like conservative bound when used with a
#'  Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#'  functions. The parameter \code{alpha_sfpar} specifies any parameters needed
#'  for the spending function specified by \code{alpha_sf}; this will be ignored
#'  for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#'  (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#'  lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#'  \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#'  bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#'  O'Brien-Fleming spending function. See details,
#'  \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param beta_sfpar Real value, default is \eqn{-4} which is an
#'  O'Brien-Fleming-like conservative bound when used with a
#'  Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#'  functions. The parameter \code{beta_sfpar} specifies any parameters needed
#'  for the spending function specified by \code{beta_sf}; this will be ignored
#'  for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#'  (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{k}{As
#'  input.} \item{test_type}{As input.} \item{outcome_type}{As input.} \item{
#'  size_type}{As input.} \item{size_type}{As input.} \item{recruit_type}{As
#'  input.} \item{delta}{As input.} \item{sigma_vec}{As input.} \item{p_vec}{As
#'  input.} \item{rho}{As input.} \item{alpha}{As input.} \item{beta}{As input.}
#'  \item{timing}{As input.} \item{i}{Fisher information at each planned interim
#'  analysis.} \item{fix_i}{Information for trial with no planned interim
#'  analyses.} \item{max_i}{Maximum information.} \item{i_frac}{Fraction of
#'  maximum information at each planned interim analysis} \item{e_i}{A vector of
#'  length 2 with expected information under the null and alternative
#'  hypotheses.} \item{m}{Number of clusters per arm at each planned interim
#'  analysis.} \item{max_m}{Maximum number of clusters per arm.} \item{e_m}{A
#'  vector of length 2 with expected number of clusters per arm under the null
#'  and alternative hypotheses.} \item{n}{Average cluster size at each planned
#'  interim analysis.} \item{max_n}{Maximum cluster size.} \item{e_m}{A vector
#'  of length 2 with expected cluster sizes under the null and alternative
#'  hypotheses.} \item{max_total}{Maximum number of individuals per arm.}
#'  \item{e_total}{A vector of length 2 with expected number of individuals per
#'  arm under the null and alternative hypotheses.} \item{sufficient}{Value
#'  denoting whether calculated sample size will be sufficient to achieve
#'  specified Type I error rate and power given the trial specifications.}
#'  \item{m_frac}{Number of clusters per arm at each interim analysis under
#'  equal information increments.} \item{n_frac}{Average cluster size at each
#'  interim analysis under equal information increments.}
#'  \item{lower_bound}{Calculated lower futility boundaries under equal
#'  information increments.} \item{upper_bound}{Calculated upper efficacy
#'  boundaries under equal information increments.} \item{tol}{As input.}
#'  \item{r}{As input.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#'
#' @useDynLib gsDesignCRT gsbounds1
#' @useDynLib gsDesignCRT gsbounds2
#'
#' @name gsDesignCRT
# gsDesignCRT function [sinew] ----
gsDesignCRT <- function(k = 3, outcome_type = 1, test_type = 1, test_sides = 1,
                        size_type = 1, recruit_type = 1,
                        delta = 1, sigma_vec = c(1, 1), p_vec = c(0.5, 0.5),
                        rho = 0, alpha = 0.05, beta = 0.1,
                        m_fix = 1, n_fix = 1, timing = 1,
                        alpha_sf = sfLDOF, alpha_sfpar = -4,
                        beta_sf = sfLDOF, beta_sfpar = -4,
                        tol = 0.000001, r = 18) {
  x <- list(k = k, outcome_type = outcome_type, test_type = test_type,
            test_sides = test_sides, size_type = size_type,
            recruit_type = recruit_type, delta = delta,
            sigma_vec = sigma_vec, p_vec = p_vec, rho = rho,
            alpha = alpha, beta = beta, timing = timing,
            i = 0, fix_i = 0, max_i = 0, i_frac = 0, e_i = 0,
            m = m_fix, max_m = m_fix, e_m = c(m_fix, m_fix),
            n = n_fix, max_n = n_fix, e_n = c(n_fix, n_fix),
            max_total = m_fix * n_fix, e_total = m_fix * n_fix,
            sufficient = 1, m_frac = 0, n_frac = 0,
            lower_bound = rep(0, k), upper_bound = rep(0, k), tol = tol, r = r)

  # Compute maximum information
  x_max <- gsMaxInfoCRT(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar)
  x$fix_i <- x_max$fix_i
  x$max_i <- x_max$i[x$k]
  x$i_frac <- x_max$i / x$max_i
  x$lower_bound <- x_max$lower_bound
  x$upper_bound <- x_max$upper_bound

  # Convert maximum information to sample sizes
  if (size_type == 1) {
    if (outcome_type == 1 && !is.null(x$sigma_vec)) {
      x$max_m <- mContDiff(x$max_i, x$max_n, x$sigma_vec, x$rho)
    } else if (outcome_type == 2 && !is.null(x$p_vec)) {
      x$max_m <- mPropDiff(x$max_i, x$max_n, x$p_vec, x$rho)
    } else {
      stop("invalid specification of outcome")
    }
  } else if (size_type == 2) {
    if (outcome_type == 1 && !is.null(x$sigma_vec)) {
      x$sufficient <- checkSuff(outcome_type = x$outcome_type,
                                m = x$max_m,
                                max_i = x$max_i,
                                sigma_vec = x$sigma_vec,
                                rho = x$rho,
                                alpha = x$alpha,
                                beta = x$beta)

      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
      }
      x$max_n <- nContDiff(x$max_i, x$max_m, x$sigma_vec, x$rho)
    } else if (outcome_type == 2 && !is.null(x$p_vec)) {
      x$sufficient <- checkSuff(outcome_type = x$outcome_type,
                                m = x$max_m,
                                max_i = x$max_i,
                                p_vec = x$p_vec,
                                rho = x$rho,
                                alpha = x$alpha,
                                beta = x$beta)
      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
      }
      x$max_n <- nPropDiff(x$max_i, x$max_m, x$p_vec, x$rho)
    } else {
      stop("invalid specification of outcome")
    }
  } else {
    stop("invalid sample size variable")
  }

  # Compute sample increments under equal information fraction increments
  x$m_frac <- sizeInfoFrac(1, x$i_frac, x$max_m, x$max_n, x$rho)
  x$n_frac <- sizeInfoFrac(2, x$i_frac, x$max_m, x$max_n, x$rho)

  # Compute information, expected information, and boundaries under equal
  # sample size increments if sufficient
  if (x$sufficient == 1) {
    x_expect <- gsExpectedInfoCRT(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar)
    x$i <- x_expect$i
    x$e_i <- x_expect$e_i

    # Convert expected information to sample sizes
    if (recruit_type == 1 || recruit_type == 3) {
      if (outcome_type == 1) {
        x$m <- mContDiff(x$i, x$max_n, x$sigma_vec, x$rho)
        x$e_m <- mContDiff(x$e_i, x$max_n, x$sigma_vec, x$rho)
        x$e_n <- x$max_n
      } else if (outcome_type == 2) {
        x$m <- mPropDiff(x$i, x$max_n, x$p_vec, x$rho)
        x$e_m <- mPropDiff(x$e_i, x$max_n, x$p_vec, x$rho)
        x$e_n <- x$max_n
      }
    } else if (recruit_type == 2) {
      if (outcome_type == 1) {
        x$n <- nContDiff(x$i, x$max_m, x$sigma_vec, x$rho)
        x$e_n <- nContDiff(x$e_i, x$max_m, x$sigma_vec, x$rho)
        x$e_m <- x$max_m
      } else if (outcome_type == 2) {
        x$n <- mPropDiff(x$i, x$max_m, x$p_vec, x$rho)
        x$e_n <- nPropDiff(x$e_i, x$max_m, x$p_vec, x$rho)
        x$e_m <- x$max_m
      }
    } else {
      stop("invalid recruitment strategy")
    }

    # Compute total number per arm
    x$total <- x$m * x$n
    x$max_total <- x$max_m * x$max_n
    x$e_total <- x$e_m * x$e_n
  }
  return(x)
}

###
# Hidden Functions
###

# gsMaxInfoCRT function [sinew] ----
gsMaxInfoCRT <- function(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar) {
  # Check inputs
  checkScalar(x$k, "integer", c(1, Inf))
  checkScalar(x$outcome_type, "integer", c(1, 2))
  checkScalar(x$test_type, "integer", c(1, 5))
  checkScalar(x$test_sides, "integer", c(1, 2))
  checkScalar(x$size_type, "integer", c(1, 2))
  checkScalar(x$recruit_type, "integer", c(1, 3))
  if (x$outcome_type == 1) {
    checkScalar(x$delta, "numeric", c(0, Inf), c(FALSE, FALSE))
    checkVector(x$sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
  } else if (x$outcome_type == 2) {
    checkScalar(x$delta, "numeric", c(0, 1), c(FALSE, TRUE))
    checkVector(x$p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
    checkScalar(abs(x$p_vec[2] - x$p_vec[1]), "numeric", x$delta)
  }
  checkScalar(x$alpha, "numeric", 0:1, c(FALSE, FALSE))
  checkScalar(x$beta, "numeric", c(0, 1 - x$alpha), c(FALSE, FALSE))
  if (x$size_type == 1) {
    checkScalar(x$n, "integer", c(0, Inf), c(FALSE, FALSE))
  } else if (x$size_type == 2) {
    checkScalar(x$m, "integer", c(0, Inf), c(FALSE, FALSE))
  }
  checkScalar(x$rho, "numeric", c(0, 1))
  checkScalar(x$tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  checkScalar(x$r, "integer", c(1, 80))

  # Specify timing depending on provided information
  if (length(x$timing) < 1 ||
        (length(x$timing) == 1 &&
           (x$k > 2 || (x$k == 2 && (x$timing[1] <= 0 || x$timing[1] >= 1))))) {
    x$timing <- seq(x$k) / x$k
  } else if (x$timing == 1 && x$k == 1) {
    x$timing <- 1
  } else if (length(x$timing) == x$k - 1 || length(x$timing) == x$k) {
    if (length(x$timing) == x$k - 1) {
      x$timing <- c(x$timing, 1)
    } else if (x$timing[x$k] != 1) {
      stop("if analysis timing for final analysis is input, it must be 1")
    }
    if (min(x$timing - c(0, x$timing[1:(x$k - 1)])) <= 0) {
      stop("input timing of interim analyses must be increasing strictly
           between 0 and 1")
    }
  } else {
    stop("value input for timing must be length 1, k-1 or k")
  }

  # Check error spending
  if (is.character(alpha_sf) &&
        !is.element(alpha_sf, c("OF", "Pocock", "WT"))) {
    stop("Character specification of upper spending may only be WT, OF or
         Pocock")
  } else if (!is.function(alpha_sf)) {
    stop("Upper spending function mis-specified")
  }
  if (is.character(beta_sf) &&
        !is.element(beta_sf, c("OF", "Pocock", "WT"))) {
    stop("Character specification of upper spending may only be WT, OF or
         Pocock")
  } else if (!is.function(beta_sf)) {
    stop("Upper spending function mis-specified")
  }

  if (x$test_type == 1) {
    # Partition error probabilities based on error spending
    if (x$test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
      i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
    } else {
      upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
      i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
    }
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$k - 1])
    falseneg <- c(rep(1e-30, x$k - 1), x$beta)
    
    binding <- TRUE
  } else if (x$test_type == 2 || x$test_type == 3) {
    # Partition error probabilities based on error spending
    if (x$test_sides == 1) {
      i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
      falsepos <- c(rep(1e-30, x$k - 1), x$alpha)
    } else {
      i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
      falsepos <- c(rep(1e-30, x$k - 1), x$alpha / 2)
    }

    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$k - 1])

    if (x$test_type == 2) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
  } else if (x$test_type == 4 || x$test_type == 5) {
    # Partition error probabilities based on error spending
    if (x$test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
      i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
    } else {
      upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
      i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
    }
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$k - 1])

    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$k - 1])

    if (x$test_type == 4) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
  } else {
    stop("invalid test type")
  }

  # Determine maximum information that achieves a_K = b_K with error spending
  x$i <- stats::uniroot(gsbounddiffCRT, lower = i0, upper = 10 * i0,
                        extendInt = "yes", theta = x$delta,
                        falseneg = falseneg, falsepos = falsepos,
                        time = x$timing, sides = x$test_sides,
                        binding = binding, tol = x$tol,
                        r = x$r)$root * x$timing

  # Calculate corresponding boundaries based on maximum information
  bounds <- gsBoundsCRT(x$delta, x$i, falseneg, falsepos, x$test_sides,
                        binding, x$tol, x$r)
  x$lower_bound <- bounds$a
  x$upper_bound <- bounds$b

  # Save information for trial with no interim analyses
  x$fix_i <- i0

  return(x)
}

# gsExpectedInfoCRT function [sinew] ----
gsExpectedInfoCRT <- function(x, alpha_sf, alpha_sfpar, beta_sf, beta_sfpar) {
  # Check additional inputs
  checkScalar(x$recruit_type, "integer", c(1, 3))

  # Specify timing depending on provided information
  if (x$recruit_type == 1) {
    x$m <- ((1:x$k) / x$k) * x$max_m
    if (x$outcome_type == 1) {
      x$i <- iContDiff(x$m, x$max_n, x$sigma_vec, x$rho)
    } else if (x$outcome_type == 2) {
      x$i <- iPropDiff(x$m, x$max_n, x$p_vec, x$rho)
    }
  } else if (x$recruit_type == 2) {
    x$n <- ((1:x$k) / x$k) * x$max_n
    if (x$outcome_type == 1) {
      x$i <- iContDiff(x$max_m, x$n, x$sigma_vec, x$rho)
    } else if (x$outcome_type == 2) {
      x$i <- iPropDiff(x$max_m, x$n, x$p_vec, x$rho)
    }
  } else if (x$recruit_type == 3) {
    x$m <- ((1:x$k) / x$k) * x$max_m
    x$n <- ((1:x$k) / x$k) * x$max_n
    if (x$outcome_type == 1) {
      x$i <- iContDiff(x$m, x$n, x$sigma_vec, x$rho)
    } else if (x$outcome_type == 2) {
      x$i <- iPropDiff(x$m, x$n, x$p_vec, x$rho)
    }
  }

  x$timing <- x$i / x$i[x$k]

  if (x$test_type == 1) {
    # Partition error probabilities based on error spending
    if (x$test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
      x$lower_bound <- rep(-30, x$k)
    } else {
      upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
      x$lower_bound <- rep(0, x$k)
    }
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$k - 1])

    # Determine boundary values to achieve Type I error = alpha
    bounds <- gsUpperCRT(0, x$i, x$lower_bound, falsepos, x$test_sides,
                         x$tol, x$r)
    x$upper_bound <- bounds$b

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprobCRT(theta, x$i, x$lower_bound, x$upper_bound, x$test_sides, x$r)
    x$e_i <- y$eI
  } else if (x$test_type == 2 || x$test_type == 3) {
    # Partition error probabilities based on error spending
    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$k - 1])

    # Determine boundary values to achieve Type II error = beta
    if (x$test_type == 3) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
    x$upper_bound <- rep(30, x$k)
    bounds <- gsLowerCRT(x$delta, x$i, falseneg, x$upper_bound, x$test_sides,
                         binding, x$tol, x$r)
    x$lower_bound <- bounds$a

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprobCRT(theta, x$i, x$lower_bound, x$upper_bound, x$test_sides, x$r)
    x$e_i <- y$eI
  } else if (x$test_type == 4 || x$test_type == 5) {
    # Partition error probabilities based on error spending
    if (x$test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$timing, alpha_sfpar)
    } else {
      upper <- alpha_sf(x$alpha / 2, x$timing, alpha_sfpar)
    }
    falsepos <- upper$spend
    falsepos <- falsepos - c(0, falsepos[1:x$k - 1])

    lower <- beta_sf(x$beta, x$timing, beta_sfpar)
    falseneg <- lower$spend
    falseneg <- falseneg - c(0, falseneg[1:x$k - 1])

    # Determine boundary values to achieve Type I error = alpha
    if (x$test_type == 5) {
      binding <- TRUE
    } else {
      binding <- FALSE
    }
    bounds <- gsBoundsCRT(x$delta, x$i, falseneg, falsepos, x$test_sides,
                          binding, x$tol, x$r)
    x$lower_bound <- bounds$a
    x$upper_bound <- bounds$b

    # Compute crossing probabilities and expected information
    theta <- as.double(c(0, x$delta)) # For H0 and H1
    y <- gsprobCRT(theta, x$i, x$lower_bound, x$upper_bound, x$test_sides, x$r)
    x$e_i <- y$eI
  } else {
    stop("invalid test type")
  }

  return(x)
}

# iContDiff function [sinew] ----
iContDiff <- function(m, n, sigma_vec, rho) {
  num <- m * n
  denom <- (sigma_vec[1]^2 + sigma_vec[2]^2) * (1 + (n - 1) * rho)
  i <- num / denom
  return(i)
}

# iPropDiff function [sinew] ----
iPropDiff <- function(m, n, p_vec, rho) {
  p <- mean(p_vec)
  num <- m * n
  denom <- (2 * p * (1 - p)) * (1 + (n - 1) * rho)
  i <- num / denom
  return(i)
}

# mContDiff function [sinew] ----
mContDiff <- function(i, n, sigma_vec, rho) {
  d_eff <- (1 + (n - 1) * rho) / n
  m <- i * (sigma_vec[1]^2 + sigma_vec[2]^2) * d_eff
  return(m)
}

# mPropDiff function [sinew] ----
mPropDiff <- function(i, n, p_vec, rho) {
  p <- mean(p_vec)
  d_eff <- (1 + (n - 1) * rho) / n
  m <- i * (2 * p * (1 - p)) * d_eff
  return(m)
}

# nContDiff function [sinew] ----
nContDiff <- function(i, m, sigma_vec, rho) {
  num <- i * (sigma_vec[1]^2 + sigma_vec[2]^2) * (1 - rho) / m
  denom <- 1 - (i * (sigma_vec[1]^2 + sigma_vec[2]^2) * rho / m)
  n <- num / denom
  return(n)
}

# nPropDiff function [sinew] ----
nPropDiff <- function(i, m, p_vec, rho) {
  p <- mean(p_vec)
  num <- i * (2 * p * (1 - p)) * (1 - rho) / m
  denom <- 1 - (i * (2 * p * (1 - p)) * rho / m)
  n <- num / denom
  return(n)
}

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