# gsSimContCRT roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial with continuous
#' outcomes.
#'
#' @param k Number of analyses planned, including interim and final.
#' @param data Simulated continuous outcomes. Should be an n x 4 matrix with the
#' columns encoding the treatment arm, cluster, individual, and response.
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param recruit_type \code{1=}by clusters with all individuals recruited
#' \cr \code{2=}by individuals in recruited cluster
#' \cr \code{3=}by both clusters and individuals in clusters
#' @param stat_type \code{1=} Z-test with known variance and ICC
#' \cr \code{2=} Z-test with re-estimated variance and ICC
#' \cr \code{3=} t-test with re-estimated variance and ICC
#' @param balance_size \code{1=} exact sample increments according to the
#' scheduled interim analyses.
#' \cr \code{2=} randomized sample increments from multinomial distribution
#' according to the scheduled interim analyses.
#' @param precompute Use pre-computed stopping boundaries if true.
#' @param delta Effect size for theta under alternative hypothesis.
#' @param sigma_vec Standard deviations for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Desired Type I error, always one-sided. Default value is 0.05.
#' @param beta Desired Type II error, default value is 0.1 (90\% power).
#' @param lower_bound Pre-computed lower futility boundaries at the specified
#' interim analyses. Must be specified if precompute is TRUE. NULL otherwise.
#' @param upper_bound Pre-computed upper efficacy boundaries at the specified
#' interim analyses. Must be specified if precompute is TRUE. NULL otherwise.
#' @param m_max Number of clusters.
#' @param n_max Mean size of each cluster.
#' @param schedule_m Number of clusters at each interim look. Interim analyses
#' will be conducted according to the information levels in \code{schedule_m}
#' and \code{schedule_n} if provided. Otherwise, interim analyses will be
#' conducted at equal-sized sample increments according to \code{recruit_type}.
#' @param schedule_n Average cluster size at each interim look. Interim
#' analyses will be conducted according to the information levels in
#' \code{schedule_m} and \code{schedule_n} if provided. Otherwise, interim
#' analyses will be conducted at equal-sized sample increments according to
#' \code{recruit_type}.
#' @param alpha_sf A spending function or a character string indicating an
#' upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param alpha_sfpar Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{alpha_sfpar} specifies any parameters needed
#' for the spending function specified by \code{alpha_sf}; this will be ignored
#' for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#' lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param beta_sfpar Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{beta_sfpar} specifies any parameters needed
#' for the spending function specified by \code{beta_sf}; this will be ignored
#' for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{reject}{Whether the
#' null hypothesis was rejected in the simulated trial.} \item{k_i}{Interim
#' analysis at which simulated trial was stopped.} \item{m_i}{Number of clusters
#' per arm when the simulated trial was stopped.} \item{n_i}{Average number of
#' individuals per cluster when the simulated trial was stopped.}
#' \item{total_i}{Total number of individuals per arm when the simulated trial
#' was stopped.} \item{i_frac}{Information fraction when the simulated trial
#' was stopped.}
#'
#' @importFrom stats rmultinom
#' @importFrom lme4 lmer
#' @importFrom performance icc
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gsSimContCRT
# gsSimContCRT function [sinew] ----
gsSimContCRT <- function(k, data, test_type, test_sides, recruit_type,
                         stat_type, balance_size, precompute = FALSE,
                         delta, sigma_vec = c(1, 1), rho,
                         alpha = 0.05, beta = 0.1,
                         lower_bound = NULL, upper_bound = NULL,
                         m_max = 1, n_max = 1,
                         schedule_m = NULL, schedule_n = NULL,
                         alpha_sf, alpha_sfpar = -4, beta_sf, beta_sfpar = -4,
                         tol = 0.000001, r = 18) {
  # Check inputs
  checkScalar(k, "integer", c(1, Inf))
  checkScalar(test_type, "integer", c(1, 5))
  checkScalar(test_sides, "integer", c(1, 2))
  checkScalar(recruit_type, "integer", c(1, 3))
  checkScalar(stat_type, "integer", c(1, 3))
  checkScalar(balance_size, "integer", c(1, 2))
  checkScalar(delta, "numeric", c(0, Inf), c(FALSE, FALSE))
  checkVector(sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
  checkScalar(rho, "numeric", c(0, 1))
  checkScalar(alpha, "numeric", 0:1, c(FALSE, FALSE))
  checkScalar(beta, "numeric", c(0, 1 - alpha), c(FALSE, FALSE))
  if (precompute) {
    checkLengths(lower_bound, upper_bound)
  }
  checkScalar(m_max, "integer", c(0, Inf), c(FALSE, FALSE))
  checkScalar(n_max, "integer", c(0, Inf), c(FALSE, FALSE))
  if (!is.null(schedule_m) && !is.null(schedule_n)) {
    checkLengths(schedule_m, schedule_n)
    if (precompute) {
      checkLengths(schedule_m, schedule_n, lower_bound, upper_bound)
    }
  }
  checkScalar(tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))

  # Set interim look parameters
  ka_i <- 0
  kb_i <- 0

  info_vec <- c()
  info_vec_ret <- c()

  ifrac_prev <- 0
  ifrac_vec <- c()
  ifrac_vec_ret <- c()

  falsepos_vec <- c()
  falseneg_vec <- c()
  stop <- FALSE

  if (is.null(schedule_m)) {
    if (balance_size == 1) {
      schedule_m1 <- ((1:k) / k) * m_max
      schedule_m2 <- ((1:k) / k) * m_max
    } else {
      schedule_m1 <- cumsum(rep(1, k) + rmultinom(1, m_max - k, rep(1 / k, k)))
      schedule_m2 <- cumsum(rep(1, k) + rmultinom(1, m_max - k, rep(1 / k, k)))
    }
  } else {
    if (balance_size == 1) {
      schedule_m1 <- schedule_m
      schedule_m2 <- schedule_m
    } else {
      schedule_m1 <- cumsum(rep(1, k) + rmultinom(1, m_max - k,
                                                  schedule_m / m_max))
      schedule_m2 <- cumsum(rep(1, k) + rmultinom(1, m_max - k,
                                                  schedule_m / m_max))
    }
  }

  if (is.null(schedule_n)) {
    schedule_n1 <- ((1:k) / k) * n_max
    schedule_n2 <- ((1:k) / k) * n_max
  } else {
    schedule_n1 <- schedule_n
    schedule_n2 <- schedule_n
  }

  while (!stop && ka_i < k) {
    # Update number of looks
    ka_i <- ka_i + 1 # number of interim looks
    kb_i <- kb_i + 1 # number of crossing bounds (usually = ka_i)

    # "Recruit" clusters
    if (recruit_type == 1) {
      ifrac_size_total1 <- ceiling(schedule_m1) * n_max
      ifrac_size_total2 <- ceiling(schedule_m2) * n_max

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max
      } else {
        m1_i <- max(ceiling(schedule_m1[ka_i]), 1)
        m2_i <- max(ceiling(schedule_m2[ka_i]), 1)

        n1_i <- n_max
        n2_i <- n_max
      }

      size_vec1_i <- rep(n1_i, m1_i)
      size_vec2_i <- rep(n1_i, m1_i)

      data1_i <- data[(data$arm == 0 & data$cluster %in% 1:m1_i), ]
      data2_i <- data[(data$arm == 1 &
                         data$cluster %in% (m_max + 1):(m_max + m2_i)), ]
      data_i <- rbind(data1_i, data2_i)
    } else if (recruit_type == 2) {
      ifrac_size_total1 <- m_max * ceiling(schedule_n1)
      ifrac_size_total2 <- m_max * ceiling(schedule_n2)

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max

        size_vec1_i <- rep(n1_i, m1_i)
        size_vec2_i <- rep(n1_i, m1_i)

        data1_i <- data[(data$arm == 0 & data$individual %in% 1:n1_i), ]
        data2_i <- data[(data$arm == 1 & data$individual %in% 1:n2_i), ]
        data_i <- rbind(data1_i, data2_i)
      } else {
        if (balance_size == 1) {
          n1_vec_i <- rep(schedule_n1[ka_i], m_max)
          n2_vec_i <- rep(schedule_n2[ka_i], m_max)
        } else {
          n1_vec_i <- rep(2, m_max) + rmultinom(1,
                                            max(ifrac_size_total1[ka_i] -
                                                  (2 * m_max), 0),
                                            rep(1 / m_max, m_max))
          n2_vec_i <- rep(2, m_max) + rmultinom(1,
                                            max(ifrac_size_total2[ka_i] -
                                                  (2 * m_max), 0),
                                            rep(1 / m_max, m_max))
        }
        m1_i <- 1
        m2_i <- m_max + 1

        size_vec1_i <- c()
        size_vec2_i <- c()

        data_i <- data.frame()
        for (m_i in 1:m_max) {
          n1_i <- n1_vec_i[m_i]
          n2_i <- n2_vec_i[m_i]

          size_vec1_i <- c(size_vec1_i, n1_i)
          size_vec2_i <- c(size_vec2_i, n1_i)

          data1_i <- data[(data$arm == 0 &
                             data$cluster == m1_i &
                             data$individual %in% 1:n1_i), ]
          data2_i <- data[(data$arm == 1 &
                             data$cluster == m2_i &
                             data$individual %in% 1:n2_i), ]
          data_i <- rbind(data_i, data1_i, data2_i)

          m1_i <- m1_i + 1
          m2_i <- m2_i + 1
        }
      }
    } else if (recruit_type == 3) {
      ifrac_size_total1 <- ceiling(schedule_m1) * ceiling(schedule_n1)
      ifrac_size_total2 <- ceiling(schedule_m2) * ceiling(schedule_n2)

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max

        size_vec1_i <- rep(n1_i, m1_i)
        size_vec2_i <- rep(n1_i, m1_i)

        data1_i <- data[(data$arm == 0 & data$individual %in% 1:n1_i), ]
        data2_i <- data[(data$arm == 1 & data$individual %in% 1:n2_i), ]
        data_i <- rbind(data1_i, data2_i)
      } else {
        m1 <- max(ceiling(schedule_m1[ka_i]), 1)
        m2 <- max(ceiling(schedule_m2[ka_i]), 1)

        if (balance_size == 1) {
          n1_vec_i <- rep(schedule_n1[ka_i], m1)
          n2_vec_i <- rep(schedule_n2[ka_i], m2)
        } else {
          n1_vec_i <- rep(2, m1) + rmultinom(1,
                                             max(ifrac_size_total1[ka_i] -
                                                   (2 * m1), 0),
                                             rep(1 / m1, m1))
          n2_vec_i <- rep(2, m2) + rmultinom(1,
                                             max(ifrac_size_total2[ka_i] -
                                                   (2 * m2), 0),
                                             rep(1 / m2, m2))
        }

        data_i <- data.frame()

        m1_i <- 1
        size_vec1_i <- c()
        for (m_i in 1:m1) {
          n1_i <- n1_vec_i[m_i]
          size_vec1_i <- c(size_vec1_i, n1_i)

          data1_i <- data[(data$arm == 0 &
                             data$cluster == m1_i &
                             data$individual %in% 1:n1_i), ]
          data_i <- rbind(data_i, data1_i)

          m1_i <- m1_i + 1
        }

        m2_i <- m_max + 1
        size_vec2_i <- c()
        for (m_i in 1:m2) {
          n2_i <- n2_vec_i[m_i]
          size_vec2_i <- c(size_vec2_i, n2_i)

          data2_i <- data[(data$arm == 1 &
                             data$cluster == m2_i &
                             data$individual %in% 1:n2_i), ]
          data_i <- rbind(data_i, data2_i)

          m2_i <- m2_i + 1
        }
      }
    }

    # Get responses
    data1_i <- data_i[data_i$arm == 0, ]
    data2_i <- data_i[data_i$arm == 1, ]

    x1_i <- data1_i$response
    x2_i <- data2_i$response

    # Compute corresponding means
    x_bar1_i <- mean(x1_i)
    x_bar2_i <- mean(x2_i)

    # Compute test statistics
    if (stat_type == 1) {
      rho_i <- rho
      se <- seContDiff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, sigma_vec)
    } else {
      fit <- lmer(response ~ as.factor(arm) + (1 | cluster), data = data_i)
      icc_est_i <- icc(fit)
      if (length(icc_est_i) > 1) {
        rho_i <- icc_est_i$ICC_adjusted
      } else {
        rho_i <- 0
      }
      se <- seContDiff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, NULL)
    }

    z_i <- (x_bar2_i - x_bar1_i) / se
    if (test_sides == 2) {
      z_i <- abs(z_i)
    }

    # Conduct group sequential test
    info <- 1 / (se^2)
    info_vec_ret <- c(info_vec_ret, info)

    ifrac <- iFrac(size_vec1_i, size_vec2_i, m_max, n_max, rho_i)
    ifrac_vec_ret <- c(ifrac_vec_ret, ifrac)

    if (ifrac <= ifrac_prev) {
      kb_i <- kb_i - 1
    } else if (test_type == 1) {
      # Get stopping bounds depending on test
      if (precompute) {
        if (stat_type <= 2) {
          ub_i <- upper_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          ub_i <- qt(p = pnorm(upper_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Alpha spending
        if (test_sides == 1) {
          spend_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
          spend_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
          a <- rep(-30, kb_i)
        } else {
          spend_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
          spend_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
          a <- rep(0, kb_i)
        }
        falsepos_i <- c(falsepos_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        bounds_i <- gsUpperCRT(0, I = c(info_vec, info), a = a,
                               falsepos = falsepos_i, sides = test_sides,
                               tol = tol, r = r)
        if (stat_type <= 2) {
          ub_i <- bounds_i$b[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
      }

    } else if (test_type == 2 || test_type == 3) {
      if (precompute) {
        if (stat_type <= 2) {
          lb_i <- lower_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(lower_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Beta spending
        spend_i <- beta_sf(beta, ifrac, beta_sfpar)
        spend_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
        falseneg_i <- c(falseneg_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 3) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
        bounds_i <- gsLowerCRT(delta, I = c(info_vec, info), b = rep(30, kb_i),
                               falseneg = falseneg_i, sides = test_sides,
                               binding = binding, tol = tol, r = r)
        if (stat_type <= 2) {
          lb_i <- bounds_i$a[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
      }
    } else if (test_type == 4 || test_type == 5) {
      if (precompute) {
        if (stat_type <= 2) {
          lb_i <- lower_bound[ka_i]
          ub_i <- upper_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(lower_bound[ka_i]), df = df_i)
          ub_i <- qt(p = pnorm(upper_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- FALSE
          stop <- TRUE
        } else if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Alpha spending
        if (test_sides == 1) {
          spend1_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
          spend1_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
        } else {
          spend1_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
          spend1_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
        }
        falsepos_i <- c(falsepos_vec, spend1_i$spend - spend1_i_prev$spend)

        # Beta spending
        spend2_i <- beta_sf(beta, ifrac, beta_sfpar)
        spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
        falseneg_i <- c(falseneg_vec, spend2_i$spend - spend2_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 5) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
        bounds_i <- gsBoundsCRT(theta = delta, I = c(info_vec, info),
                                falseneg = falseneg_i, falsepos = falsepos_i,
                                sides = test_sides, binding = binding,
                                tol = tol, r = r)
        if (stat_type <= 2) {
          lb_i <- bounds_i$a[kb_i]
          ub_i <- bounds_i$b[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
          ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- FALSE
          stop <- TRUE
        } else if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
        falseneg_vec <- falseneg_i
      }
    }
  }

  # Return results
  m_i <- mean(length(size_vec1_i), length(size_vec2_i))
  n_i <- mean(c(size_vec1_i, size_vec2_i))
  total_i <- sum(size_vec1_i) + sum(size_vec2_i)
  out <- list("reject" = reject,
              "t_i" = ka_i,
              "m_i" = m_i,
              "n_i" = n_i,
              "total_i" = total_i,
              "i_frac" = ifrac_prev)
  return(out)
}

# gsSimBinCRT roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial with binary
#' outcomes
#'
#' @param k Number of analyses planned, including interim and final.
#' @param data Simulated binary outcomes. Should be an n x 4 matrix with the
#' columns encoding the treatment arm, cluster, individual, and response.
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param recruit_type \code{1=}by clusters with all individuals recruited
#' \cr \code{2=}by individuals in recruited cluster
#' \cr \code{3=}by both clusters and individuals in clusters
#' @param stat_type \code{1=} Z-test with known variance and ICC
#' \cr \code{2=} Z-test with re-estimated variance and ICC
#' \cr \code{3=} t-test with re-estimated variance and ICC
#' @param balance_size \code{1=} exact sample increments according to the
#' scheduled interim analyses.
#' \cr \code{2=} randomized sample increments from multinomial distribution
#' according to the scheduled interim analyses.
#' @param precompute Use pre-computed stopping boundaries if true.
#' @param delta Effect size for theta under alternative hypothesis.
#' @param p_vec Probabilities of event for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Type I error, always one-sided. Default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#' @param lower_bound Pre-computed lower futility boundaries at the specified
#' interim analyses. Must be specified if precompute is TRUE. NULL otherwise.
#' @param upper_bound Pre-computed upper efficacy boundaries at the specified
#' interim analyses. Must be specified if precompute is TRUE. NULL otherwise.
#' @param m_max Number of clusters.
#' @param n_max Mean size of each cluster.
#' @param schedule_m Number of clusters at each interim look. Interim analyses
#' will be conducted according to the information levels in \code{schedule_m}
#' and \code{schedule_n} if provided. Otherwise, interim analyses will be
#' conducted at equal-sized sample increments according to \code{recruit_type}.
#' @param schedule_n Average cluster size at each interim look. Interim
#' analyses will be conducted according to the information levels in
#' \code{schedule_m} and \code{schedule_n} if provided. Otherwise, interim
#' analyses will be conducted at equal-sized sample increments according to
#' \code{recruit_type}.
#' @param alpha_sf A spending function or a character string indicating an
#' upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param alpha_sfpar Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{alpha_sfpar} specifies any parameters needed
#' for the spending function specified by \code{alpha_sf}; this will be ignored
#' for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#' lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sfLDOF} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("SpendingFunctionOverview")}, manual and examples.
#' @param beta_sfpar Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{beta_sfpar} specifies any parameters needed
#' for the spending function specified by \code{beta_sf}; this will be ignored
#' for spending functions (\code{sfLDOF}, \code{sfLDPocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{reject}{Whether the
#' null hypothesis was rejected in the simulated trial.} \item{k_i}{Interim
#' analysis at which simulated trial was stopped.} \item{m_i}{Number of clusters
#' per arm when the simulated trial was stopped.} \item{n_i}{Average number of
#' individuals per cluster when the simulated trial was stopped.}
#' \item{total_i}{Total number of individuals per arm when the simulated trial
#' was stopped.} \item{i_frac}{Information fraction when the simulated trial
#' was stopped.}
#'
#' @importFrom stats rmultinom
#' @importFrom lme4 lmer
#' @importFrom performance icc
#'
#' @export
#' @name gsSimBinCRT
# gsSimBinCRT function [sinew] ----
gsSimBinCRT <- function(k, data, test_type, test_sides, recruit_type,
                        stat_type, balance_size, precompute = FALSE,
                        delta, p_vec, rho, alpha = 0.05, beta = 0.1,
                        lower_bound = NULL, upper_bound = NULL,
                        m_max = 1, n_max = 1,
                        schedule_m = NULL, schedule_n = NULL,
                        alpha_sf, alpha_sfpar = -4, beta_sf, beta_sfpar = -4,
                        tol = 0.000001, r = 18) {
  # Check inputs
  checkScalar(k, "integer", c(1, Inf))
  checkScalar(test_type, "integer", c(1, 5))
  checkScalar(test_sides, "integer", c(1, 2))
  checkScalar(recruit_type, "integer", c(1, 3))
  checkScalar(stat_type, "integer", c(1, 3))
  checkScalar(balance_size, "integer", c(1, 2))
  checkScalar(delta, "numeric", c(0, 1), c(FALSE, TRUE))
  checkVector(p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  checkScalar(rho, "numeric", c(0, 1))
  checkScalar(alpha, "numeric", 0:1, c(FALSE, FALSE))
  checkScalar(beta, "numeric", c(0, 1 - alpha), c(FALSE, FALSE))
  if (precompute) {
    checkLengths(lower_bound, upper_bound)
  }
  checkScalar(m_max, "integer", c(0, Inf), c(FALSE, FALSE))
  checkScalar(n_max, "integer", c(0, Inf), c(FALSE, FALSE))
  if (!is.null(schedule_m) && !is.null(schedule_n)) {
    checkLengths(schedule_m, schedule_n)
    if (precompute) {
      checkLengths(schedule_m, schedule_n, lower_bound, upper_bound)
    }
  }
  checkScalar(tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))

  # Set interim look parameters
  ka_i <- 0
  kb_i <- 0

  info_vec <- c()
  info_vec_ret <- c()

  ifrac_prev <- 0
  ifrac_vec <- c()
  ifrac_vec_ret <- c()

  falsepos_vec <- c()
  falseneg_vec <- c()
  stop <- FALSE

  if (is.null(schedule_m)) {
    if (balance_size == 1) {
      schedule_m1 <- ((1:k) / k) * m_max
      schedule_m2 <- ((1:k) / k) * m_max
    } else {
      schedule_m1 <- cumsum(rep(1, k) + rmultinom(1, m_max - k, rep(1 / k, k)))
      schedule_m2 <- cumsum(rep(1, k) + rmultinom(1, m_max - k, rep(1 / k, k)))
    }
  } else {
    if (balance_size == 1) {
      schedule_m1 <- schedule_m
      schedule_m2 <- schedule_m
    } else {
      schedule_m1 <- cumsum(rep(1, k) + rmultinom(1, m_max - k,
                                                  schedule_m / m_max))
      schedule_m2 <- cumsum(rep(1, k) + rmultinom(1, m_max - k,
                                                  schedule_m / m_max))
    }
  }

  if (is.null(schedule_n)) {
    schedule_n1 <- ((1:k) / k) * n_max
    schedule_n2 <- ((1:k) / k) * n_max
  } else {
    schedule_n1 <- schedule_n
    schedule_n2 <- schedule_n
  }

  while (!stop && ka_i < k) {
    # Update number of looks
    ka_i <- ka_i + 1 # number of interim looks
    kb_i <- kb_i + 1 # number of crossing bounds (usually = ka_i)

    # "Recruit" clusters
    if (recruit_type == 1) {
      ifrac_size_total1 <- ceiling(schedule_m1) * n_max
      ifrac_size_total2 <- ceiling(schedule_m2) * n_max

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max
      } else {
        m1_i <- max(ceiling(schedule_m1[ka_i]), 1)
        m2_i <- max(ceiling(schedule_m2[ka_i]), 1)

        n1_i <- n_max
        n2_i <- n_max
      }

      size_vec1_i <- rep(n1_i, m1_i)
      size_vec2_i <- rep(n1_i, m1_i)

      data1_i <- data[(data$arm == 0 & data$cluster %in% 1:m1_i), ]
      data2_i <- data[(data$arm == 1 &
                         data$cluster %in% (m_max + 1):(m_max + m2_i)), ]
      data_i <- rbind(data1_i, data2_i)
    } else if (recruit_type == 2) {
      ifrac_size_total1 <- m_max * ceiling(schedule_n1)
      ifrac_size_total2 <- m_max * ceiling(schedule_n2)

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max

        size_vec1_i <- rep(n1_i, m1_i)
        size_vec2_i <- rep(n1_i, m1_i)

        data1_i <- data[(data$arm == 0 & data$individual %in% 1:n1_i), ]
        data2_i <- data[(data$arm == 1 & data$individual %in% 1:n2_i), ]
        data_i <- rbind(data1_i, data2_i)
      } else {
        if (balance_size == 1) {
          n1_vec_i <- rep(schedule_n1[ka_i], m_max)
          n2_vec_i <- rep(schedule_n2[ka_i], m_max)
        } else {
          n1_vec_i <- rep(2, m_max) + rmultinom(1,
                                                max(ifrac_size_total1[ka_i] -
                                                      (2 * m_max), 0),
                                                rep(1 / m_max, m_max))
          n2_vec_i <- rep(2, m_max) + rmultinom(1,
                                                max(ifrac_size_total2[ka_i] -
                                                      (2 * m_max), 0),
                                                rep(1 / m_max, m_max))
        }
        m1_i <- 1
        m2_i <- m_max + 1

        size_vec1_i <- c()
        size_vec2_i <- c()

        data_i <- data.frame()
        for (m_i in 1:m_max) {
          n1_i <- n1_vec_i[m_i]
          n2_i <- n2_vec_i[m_i]

          size_vec1_i <- c(size_vec1_i, n1_i)
          size_vec2_i <- c(size_vec2_i, n1_i)

          data1_i <- data[(data$arm == 0 &
                             data$cluster == m1_i &
                             data$individual %in% 1:n1_i), ]
          data2_i <- data[(data$arm == 1 &
                             data$cluster == m2_i &
                             data$individual %in% 1:n2_i), ]
          data_i <- rbind(data_i, data1_i, data2_i)

          m1_i <- m1_i + 1
          m2_i <- m2_i + 1
        }
      }
    } else if (recruit_type == 3) {
      ifrac_size_total1 <- ceiling(schedule_m1) * ceiling(schedule_n1)
      ifrac_size_total2 <- ceiling(schedule_m2) * ceiling(schedule_n2)

      if (ka_i == k) {
        m1_i <- m_max
        m2_i <- m_max

        n1_i <- n_max
        n2_i <- n_max

        size_vec1_i <- rep(n1_i, m1_i)
        size_vec2_i <- rep(n1_i, m1_i)

        data1_i <- data[(data$arm == 0 & data$individual %in% 1:n1_i), ]
        data2_i <- data[(data$arm == 1 & data$individual %in% 1:n2_i), ]
        data_i <- rbind(data1_i, data2_i)
      } else {
        m1 <- max(ceiling(schedule_m1[ka_i]), 1)
        m2 <- max(ceiling(schedule_m2[ka_i]), 1)

        if (balance_size == 1) {
          n1_vec_i <- rep(schedule_n1[ka_i], m1)
          n2_vec_i <- rep(schedule_n2[ka_i], m2)
        } else {
          n1_vec_i <- rep(2, m1) + rmultinom(1,
                                             max(ifrac_size_total1[ka_i] -
                                                   (2 * m1), 0),
                                             rep(1 / m1, m1))
          n2_vec_i <- rep(2, m2) + rmultinom(1,
                                             max(ifrac_size_total2[ka_i] -
                                                   (2 * m2), 0),
                                             rep(1 / m2, m2))
        }

        data_i <- data.frame()

        m1_i <- 1
        size_vec1_i <- c()
        for (m_i in 1:m1) {
          n1_i <- n1_vec_i[m_i]
          size_vec1_i <- c(size_vec1_i, n1_i)

          data1_i <- data[(data$arm == 0 &
                             data$cluster == m1_i &
                             data$individual %in% 1:n1_i), ]
          data_i <- rbind(data_i, data1_i)

          m1_i <- m1_i + 1
        }

        m2_i <- m_max + 1
        size_vec2_i <- c()
        for (m_i in 1:m2) {
          n2_i <- n2_vec_i[m_i]
          size_vec2_i <- c(size_vec2_i, n2_i)

          data2_i <- data[(data$arm == 1 &
                             data$cluster == m2_i &
                             data$individual %in% 1:n2_i), ]
          data_i <- rbind(data_i, data2_i)

          m2_i <- m2_i + 1
        }
      }
    }

    # Get responses
    data1_i <- data_i[data_i$arm == 0, ]
    data2_i <- data_i[data_i$arm == 1, ]

    x1_i <- data1_i$response
    x2_i <- data2_i$response

    # Compute corresponding proportions
    p_hat1_i <- mean(x1_i)
    p_hat2_i <- mean(x2_i)

    # Compute test statistics
    if (stat_type == 1) {
      rho_i <- rho
      se <- seBinDiff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, p_vec)
    } else {
      fit <- lmer(response ~ as.factor(arm) + (1 | cluster), data = data_i)
      icc_est_i <- icc(fit)
      if (length(icc_est_i) > 1) {
        rho_i <- icc_est_i$ICC_adjusted
      } else {
        rho_i <- 0
      }
      se <- seBinDiff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, NULL)
    }

    z_i <- (p_hat2_i - p_hat1_i) / se
    if (test_sides == 2) {
      z_i <- abs(z_i)
    }

    # Conduct group sequential test
    info <- 1 / (se^2)
    info_vec_ret <- c(info_vec_ret, info)

    ifrac <- iFrac(size_vec1_i, size_vec2_i, m_max, n_max, rho_i)
    ifrac_vec_ret <- c(ifrac_vec_ret, ifrac)

    if (ifrac <= ifrac_prev) {
      kb_i <- kb_i - 1
    } else if (test_type == 1) {
      # Get stopping bounds depending on test
      if (precompute) {
        if (stat_type <= 2) {
          ub_i <- upper_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          ub_i <- qt(p = pnorm(upper_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Alpha spending
        if (test_sides == 1) {
          spend_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
          spend_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
          a <- rep(-30, kb_i)
        } else {
          spend_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
          spend_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
          a <- rep(0, kb_i)
        }
        falsepos_i <- c(falsepos_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        bounds_i <- gsUpperCRT(0, I = c(info_vec, info), a = a,
                               falsepos = falsepos_i, sides = test_sides,
                               tol = tol, r = r)
        if (stat_type <= 2) {
          ub_i <- bounds_i$b[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
      }

    } else if (test_type == 2 || test_type == 3) {
      if (precompute) {
        if (stat_type <= 2) {
          lb_i <- lower_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(lower_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Beta spending
        spend_i <- beta_sf(beta, ifrac, beta_sfpar)
        spend_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
        falseneg_i <- c(falseneg_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 3) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
        bounds_i <- gsLowerCRT(delta, I = c(info_vec, info), b = rep(30, kb_i),
                               falseneg = falseneg_i, sides = test_sides,
                               binding = binding, tol = tol, r = r)
        if (stat_type <= 2) {
          lb_i <- bounds_i$a[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
      }
    } else if (test_type == 4 || test_type == 5) {
      if (precompute) {
        if (stat_type <= 2) {
          lb_i <- lower_bound[ka_i]
          ub_i <- upper_bound[ka_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(lower_bound[ka_i]), df = df_i)
          ub_i <- qt(p = pnorm(upper_bound[ka_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- FALSE
          stop <- TRUE
        } else if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }
      } else {
        # Alpha spending
        if (test_sides == 1) {
          spend1_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
          spend1_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
        } else {
          spend1_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
          spend1_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
        }
        falsepos_i <- c(falsepos_vec, spend1_i$spend - spend1_i_prev$spend)

        # Beta spending
        spend2_i <- beta_sf(beta, ifrac, beta_sfpar)
        spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
        falseneg_i <- c(falseneg_vec, spend2_i$spend - spend2_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 5) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
        bounds_i <- gsBoundsCRT(theta = delta, I = c(info_vec, info),
                                falseneg = falseneg_i, falsepos = falsepos_i,
                                sides = test_sides, binding = binding,
                                tol = tol, r = r)
        if (stat_type <= 2) {
          lb_i <- bounds_i$a[kb_i]
          ub_i <- bounds_i$b[kb_i]
        } else {
          df_i <- dim(data_i)[1] - 2
          lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
          ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
        }

        # Conduct group sequential test for interim look
        if (z_i <= lb_i) {
          reject <- FALSE
          stop <- TRUE
        } else if (z_i >= ub_i) {
          reject <- TRUE
          stop <- TRUE
        } else {
          reject <- FALSE
        }

        # Update previous information fractions
        info_vec <- c(info_vec, info)
        ifrac_prev <- ifrac
        ifrac_vec <- c(ifrac_vec, ifrac)
        falsepos_vec <- falsepos_i
        falseneg_vec <- falseneg_i
      }
    }
  }

  # Return results
  m_i <- mean(length(size_vec1_i), length(size_vec2_i))
  n_i <- mean(c(size_vec1_i, size_vec2_i))
  total_i <- sum(size_vec1_i) + sum(size_vec2_i)
  out <- list("reject" = reject,
              "t_i" = ka_i,
              "m_i" = m_i,
              "n_i" = n_i,
              "total_i" = total_i,
              "i_frac" = ifrac_prev)
  return(out)
}

# genContCRT roxy [sinew] ----
#' @title Simulate cluster-randomized trial data with continuous outcomes
#'
#' @param m Number of clusters.
#' @param n Mean size of each cluster.
#' @param mu_vec Vector of means for control and treatment groups, respectively.
#' @param sigma_vec Vector of standard deviations for control and treatment
#'  groups, respectively.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#'
#' @return Simulated continuous outcomes represented as an n x 4 matrix with the
#'  columns encoding the treatment arm, cluster, individual, and response.
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @importFrom stats rnorm
#'
#' @export
#' @name genContCRT
# genContCRT function [sinew] ----
genContCRT <- function(m = 1, n = 1, mu_vec = c(0, 1), sigma_vec = c(1, 1),
                       rho = 0) {
  # Check inputs
  checkScalar(m, "integer", c(0, Inf), c(FALSE, FALSE))
  checkScalar(n, "integer", c(0, Inf), c(FALSE, FALSE))
  checkVector(mu_vec, "numeric", c(-Inf, Inf), c(FALSE, FALSE), length = 2)
  checkVector(sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
  checkScalar(rho, "numeric", c(0, 1))

  # Specify population means
  mu1 <- mu_vec[1]
  mu2 <- mu_vec[2]

  # Specify between and within-cluster variances
  sigma_b_1 <- sqrt(rho * sigma_vec[1]^2)
  sigma_w_1 <- sqrt(sigma_vec[1]^2 - sigma_b_1^2)
  sigma_b_2 <- sqrt(rho * sigma_vec[2]^2)
  sigma_w_2 <- sqrt(sigma_vec[2]^2 - sigma_b_2^2)

  # Generate samples by cluster
  a1 <- rep(0, m * n)
  r1 <- rep(0, m * n)
  c1 <- rep(1:m, each = n)
  i1 <- rep(1:n, m)

  a2 <- rep(1, m * n)
  r2 <- rep(0, m * n)
  c2 <- rep(m + 1:m, each = n)
  i2 <- rep(1:n, m)

  for (mi in 1:m) {
    re1 <- rnorm(1, mean = 0, sd = sigma_b_1)
    r1[(mi - 1) * n + 1:n] <- (mu1 + re1 + rnorm(n, mean = 0, sd = sigma_w_1))

    re2 <- rnorm(1, mean = 0, sd = sigma_b_2)
    r2[(mi - 1) * n + 1:n] <- (mu2 + re2 + rnorm(n, mean = 0, sd = sigma_w_2))
  }

  # Reshape samples into dataframe with cluster assignment
  df1 <- cbind.data.frame(a1, c1, i1, r1)
  colnames(df1) <- c("arm", "cluster", "individual", "response")

  df2 <- cbind.data.frame(a2, c2, i2, r2)
  colnames(df2) <- c("arm", "cluster", "individual", "response")

  # Combine samples
  df <- rbind(df1, df2)
  df <- df[order(df[, 1], df[, 2]), ]
  return(df)
}

# genBinCRT roxy [sinew] ----
#' @title Simulate cluster-randomized trial data with binary outcomes
#'
#' @param m Number of clusters.
#' @param n Mean size of each cluster.
#' @param p_vec Probabilities of event for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#'
#' @return Simulated binary outcomes represented as an n x 4 matrix with the
#'  columns encoding the treatment arm, cluster, individual, and response.
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name genBinCRT
# genBinCRT function [sinew] ----
genBinCRT <- function(m = 1, n = 1, p_vec = c(0.5, 0.5), rho = NULL) {
  # Check inputs
  checkScalar(m, "integer", c(0, Inf), c(FALSE, FALSE))
  checkScalar(n, "integer", c(0, Inf), c(FALSE, FALSE))
  checkVector(p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  checkScalar(rho, "numeric", c(0, 1))

  # Generate samples by cluster
  a1 <- rep(0, m * n)
  r1 <- rep(0, m * n)
  c1 <- rep(1:m, each = n)
  i1 <- rep(1:n, m)

  a2 <- rep(1, m * n)
  r2 <- rep(0, m * n)
  c2 <- rep(m + 1:m, each = n)
  i2 <- rep(1:n, m)

  # Generate data according to Qaqish 2003 paper
  b1 <- simB(p_vec[1], n, rho)
  b2 <- simB(p_vec[2], n, rho)

  for (mi in 1:m) {
    r1[(mi - 1) * n + 1:n] <- simResponse(p_vec[1], n, b1)
    r2[(mi - 1) * n + 1:n] <- simResponse(p_vec[2], n, b2)
  }

  # Reshape samples into dataframe with cluster assignment
  df1 <- cbind.data.frame(a1, c1, i1, r1)
  colnames(df1) <- c("arm", "cluster", "individual", "response")

  df2 <- cbind.data.frame(a2, c2, i2, r2)
  colnames(df2) <- c("arm", "cluster", "individual", "response")

  # Combine samples
  df <- rbind(df1, df2)
  df <- df[order(df[, 1], df[, 2]), ]
  return(df)
}

# simB function [sinew] ----
simB <- function(p, n, rho) {
  B <- NULL
  R <- (1 - rho) * diag(n) + rho * matrix(1, nrow = n, ncol = n)
  u <- rep(p, each = n)
  A_half <- diag(sqrt(u * (1 - u)))
  v <- A_half %*% R %*% A_half
  b <- v
  for (f in 2:n) {
    f1 <- f - 1
    gf <- v[1:f1, 1:f1]
    sf <- v[1:f1, f]
    bf <- solve(gf, sf)
    b[1:f1, f] <- bf
  }
  B <- cbind(B, b)
  return(B)
}

# simResponse function [sinew] ----
#' @importFrom stats rbinom
# simResponse function [sinew] ----
simResponse <- function(p, n, B) {
  u <- rep(p, each = n)
  y_out <- rep(-1, n)
  y_out[1] <- rbinom(1, 1, u[1])
  for (l in 2:n) {
    l1 <- l - 1
    res <- y_out[1:l1] - u[1:l1]
    cl <- u[l] + sum(res * B[1:l1, l])
    y_out[l] <- rbinom(1, 1, cl)
  }
  y <- c(y_out)
  return(y)
}

# iFrac function [sinew] ----
iFrac <- function(size_vec1, size_vec2, m_max, n_max, rho) {
  n_sq_vec1 <- size_vec1^2
  n_sq_vec2 <- size_vec2^2
  nif <- 2 * (1 + (n_max - 1) * rho) / (m_max * n_max)
  dif <- ((1 + ((sum(n_sq_vec1) / sum(size_vec1)) - 1) * rho) /
            sum(size_vec1) +
            (1 + ((sum(n_sq_vec2) / sum(size_vec2)) - 1) * rho) /
              sum(size_vec2))
  return(nif / dif)
}

# genContCRT roxy [sinew] ----
#' @importFrom stats var
# seContDiff function [sinew] ----
seContDiff <- function(x1, x2, size_vec1, size_vec2, rho, sigma_vec = NULL) {
  n_sq_vec1 <- size_vec1^2
  n_sq_vec2 <- size_vec2^2

  if (is.null(sigma_vec)) {
    pooled_num <- ((sum(size_vec1) - 1) * var(x1) +
                     (sum(size_vec2) - 1) * var(x2))
    pooled_denom <- (sum(size_vec1) - 1) + (sum(size_vec2) - 1)
    pooled <- pooled_num / pooled_denom
  } else {
    pooled_num <- ((sum(size_vec1) - 1) * sigma_vec[1]^2 +
                     (sum(size_vec2) - 1) * sigma_vec[2]^2)
    pooled_denom <- (sum(size_vec1) - 1) + (sum(size_vec2) - 1)
    pooled <- pooled_num / pooled_denom
  }
  se <- sqrt(pooled *
               (((1 + ((sum(n_sq_vec1) / sum(size_vec1)) - 1) * rho) /
                   sum(size_vec1)) +
                  ((1 + ((sum(n_sq_vec2) / sum(size_vec2)) - 1) * rho) /
                     sum(size_vec2))))
  return(se)
}

# seBinDiff function [sinew] ----
seBinDiff <- function(x1, x2, size_vec1, size_vec2, rho, p_vec = NULL) {
  n_sq_vec1 <- size_vec1^2
  n_sq_vec2 <- size_vec2^2

  if (is.null(p_vec)) {
    p <- (sum(x1) + sum(x2)) / (length(x1) + length(x2))
  } else {
    p <- mean(p_vec)
  }
  if ((p * (1 - p)) == 0) {
    return(1)
  } else {
    se <- sqrt(p * (1 - p) *
                 (((1 + ((sum(n_sq_vec1) / sum(size_vec1)) - 1) * rho) /
                     sum(size_vec1)) +
                    ((1 + ((sum(n_sq_vec2) / sum(size_vec2)) - 1) * rho) /
                       sum(size_vec2))))
    return(se)
  }
}
