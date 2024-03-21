# gsSimBinCRT roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial with continuous
#' outcomes
#'
#' @param t Number of analyses planned, including interim and final.
#' @param data Simulated trial data.
#' @param test_type \code{1=}one-sided test with early stopping for efficacy
#' \cr \code{2=}two-sided symmetric with early stopping for efficacy
#' \cr \code{3=}one-sided with early stopping for efficacy or futility
#' \cr \code{4=}two-sided with early stopping for efficacy or futility
#' @param recruit_type \code{1=}clusters \cr \code{2=}individuals per cluster
#' @param m Number of clusters.
#' @param n Mean size of each cluster.
#' @param ifrac_size Sample size at each interim look
#' @param alpha Type I error
#' @param beta Type II error
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
#' @return Results
#'
#' @importFrom lme4 glmer
#' @importFrom performance icc
#'
#' @export
#' @name gsSimBinCRT
# gsSimBinCRT function [sinew] ----
gsSimBinCRT <- function(t, data, test_type, recruit_type, icc_type, stat_type,
                        delta, p_vec, rho, alpha = 0.05, beta = 0.1, m, n,
                        ifrac_size_m = NULL, ifrac_size_n = NULL,
                        alpha_sf, alpha_sfpar = -4, beta_sf, beta_sfpar = -4,
                        tol = 0.000001, r = 42) {
  # Set interim look parameters
  ta_i <- 0
  tb_i <- 0

  info_vec <- c()
  info_vec_ret <- c()

  ifrac_prev <- 0
  ifrac_vec <- c()
  ifrac_vec_ret <- c()

  falsepos_vec <- c()
  falseneg_vec <- c()
  stop <- FALSE

  while (!stop && ta_i < t) {
    # Update number of looks
    ta_i <- ta_i + 1 # number of interim looks
    tb_i <- tb_i + 1 # number of crossing bounds (usually = ta_i)

    # "Recruit" clusters
    if (recruit_type == 1) {
      if (is.null(ifrac_size_m)) {
        ifrac_size_m <- ((1:t) / t) * m
      }

      if (ta_i == t) {
        m_i <- m
        n_i <- n
      } else {
        m_i <- max(ceiling(ifrac_size_m[ta_i]), 1)
        n_i <- n
      }
      data1_i <- data[(data$arm == 0 & data$cluster %in% 1:m_i), ]
      data2_i <- data[(data$arm == 1 & data$cluster %in% (m + 1):(m + m_i)), ]
      data_i <- rbind(data1_i, data2_i)
    } else if (recruit_type == 2) {
      if (is.null(ifrac_size_n)) {
        ifrac_size_n <- ((1:t) / t) * n
      }

      if (ta_i == t) {
        m_i <- m
        n_i <- n
      } else {
        m_i <- m
        n_i <- max(ceiling(ifrac_size_n[ta_i]), 1)
      }
      data1_i <- data[(data$arm == 0 & data$individual %in% 1:n_i), ]
      data2_i <- data[(data$arm == 1 & data$individual %in% 1:n_i), ]
      data_i <- rbind(data1_i, data2_i)
    } else if (recruit_type == 3) {
      if (is.null(ifrac_size_m)) {
        ifrac_size_m <- ((1:t) / t) * m
      }
      if (is.null(ifrac_size_n)) {
        ifrac_size_n <- ((1:t) / t) * n
      }

      if (ta_i == t) {
        m_i <- m
        n_i <- n
      } else {
        m_i <- max(ceiling(ifrac_size_m[ta_i]), 1)
        n_i <- max(ceiling(ifrac_size_n[ta_i]), 1)
      }
      data1_i <- data[(data$arm == 0 &
                         data$cluster %in% 1:m_i &
                         data$individual %in% 1:n_i), ]
      data2_i <- data[(data$arm == 1 &
                         data$cluster %in% (m + 1):(m + m_i) &
                         data$individual %in% 1:n_i), ]
      data_i <- rbind(data1_i, data2_i)
    }

    # Get responses
    x1_i <- data1_i$response
    x2_i <- data2_i$response

    # Compute corresponding proportions
    p_hat1_i <- mean(x1_i)
    p_hat2_i <- mean(x2_i)

    # Estimate ICC if needed
    rho_i <- rho
    if (icc_type == 2) {
      fit <- glmer(response ~ as.factor(arm) + (1 | cluster),
                   family = binomial(), data = data_i)
      icc_est_i <- icc(fit)
      if (length(icc_est_i) > 1) {
        rho_i <- icc_est_i$ICC_adjusted
      } else {
        rho_i <- 0
      }
    }

    # Compute test statistics
    if (stat_type == 1) {
      se <- seBinDiff(x1_i, x2_i, m_i, n_i, rho_i, p_vec)
    } else {
      se <- seBinDiff(x1_i, x2_i, m_i, n_i, rho_i, NULL)
    }
    z_i <- (p_hat2_i - p_hat1_i) / se

    # Conduct group sequential test
    info <- 1 / (se^2)
    info_vec_ret <- c(info_vec_ret, info)

    ifrac <- iFrac(m_i, n_i, m, n, rho_i)
    ifrac_vec_ret <- c(ifrac_vec_ret, ifrac)

    if (ifrac <= ifrac_prev) {
      tb_i <- tb_i - 1
    } else if (test_type == 1) {
      # Alpha spending
      spend_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
      spend_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
      falsepos_i <- c(falsepos_vec, spend_i$spend - spend_i_prev$spend)

      # Get stopping bounds depending on test
      bounds_i <- gsUpper1(0, I = c(info_vec, info), a = rep(-30, tb_i),
                           falsepos = falsepos_i, tol = tol, r = r)
      if (stat_type <= 2) {
        ub_i <- bounds_i$b[tb_i]
      } else {
        ub_i <- qt(p = pnorm(bounds_i$b[tb_i]), df = (2 * m_i * n_i) - 2)
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

    } else if (test_type == 2) {
      # Alpha spending
      spend_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
      spend_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
      falsepos_i <- c(falsepos_vec, spend_i$spend - spend_i_prev$spend)

      # Get stopping bounds depending on test
      bounds_i <- gsUpper2(0, I = c(info_vec, info), a = rep(0, tb_i),
                           falsepos = falsepos_i, tol = tol, r = r)
      if (stat_type <= 2) {
        ub_i <- bounds_i$b[tb_i]
      } else {
        df_i <- dim(data_i)[1] - 2
        ub_i <- qt(p = pnorm(bounds_i$b[tb_i]), df = df_i)
      }

      # Conduct group sequential test for interim look
      if (abs(z_i) >= ub_i) {
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

    } else if (test_type == 3) {
      # Alpha spending
      spend1_i <- alpha_sf(alpha, ifrac, alpha_sfpar)
      spend1_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sfpar)
      falsepos_i <- c(falsepos_vec, spend1_i$spend - spend1_i_prev$spend)

      # Beta spending
      spend2_i <- beta_sf(beta, ifrac, beta_sfpar)
      spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
      falseneg_i <- c(falseneg_vec, spend2_i$spend - spend2_i_prev$spend)

      # Get stopping bounds depending on test
      bounds_i <- gsBounds1(theta = delta, I = c(info_vec, info),
                            falseneg = falseneg_i, falsepos = falsepos_i,
                            tol = tol, r = r)
      if (stat_type <= 2) {
        lb_i <- bounds_i$a[tb_i]
        ub_i <- bounds_i$b[tb_i]
      } else {
        df_i <- dim(data_i)[1] - 2
        lb_i <- qt(p = pnorm(bounds_i$a[tb_i]), df = df_i)
        ub_i <- qt(p = pnorm(bounds_i$b[tb_i]), df = df_i)
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
    } else if (test_type == 4) {
      # Alpha spending
      spend1_i <- alpha_sf(alpha / 2, ifrac, alpha_sfpar)
      spend1_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sfpar)
      falsepos_i <- c(falsepos_vec, spend1_i$spend - spend1_i_prev$spend)

      # Beta spending
      spend2_i <- beta_sf(beta, ifrac, beta_sfpar)
      spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sfpar)
      falseneg_i <- c(falseneg_vec, spend2_i$spend - spend2_i_prev$spend)

      # Get stopping bounds depending on test
      bounds_i <- gsBounds2(theta = delta, I = c(info_vec, info),
                            falseneg = falseneg_i, falsepos = falsepos_i,
                            tol = tol, r = r)
      if (stat_type <= 2) {
        lb_i <- bounds_i$a[tb_i]
        ub_i <- bounds_i$b[tb_i]
      } else {
        df_i <- dim(data_i)[1] - 2
        lb_i <- qt(p = pnorm(bounds_i$a[tb_i]), df = df_i)
        ub_i <- qt(p = pnorm(bounds_i$b[tb_i]), df = df_i)
      }

      # Conduct group sequential test at interim look
      if (abs(z_i) <= lb_i) {
        reject <- FALSE
        stop <- TRUE
      } else if (abs(z_i) >= ub_i) {
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

  # Return results
  out <- list("reject" = reject,
              "t_i" = ta_i,
              "m_i" = m_i,
              "n_i" = n_i,
              "i_frac" = ifrac_prev,
              "rho" = rho_i,
              "ifrac_vec" = c(ifrac_vec_ret, rep(0, t - length(ifrac_vec_ret))))
  return(out)
}

# gsBinCRT roxy [sinew] ----
#' @title Simulate cluster-randomized trial data with continuous outcomes
#'
#' @param m Number of clusters.
#' @param n Mean size of each cluster.
#' @param mu_vec Vector of means for groups 1 and 2, respectively.
#' @param sigma_vec Vector of standard deviations for groups 1 and 2,
#' respectively.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param var_b Between-class variance.
#' @param var_w Within-class variance.
#' @param c Constant to convert probit to expit
#'
#' @return Results
#'
#' @export
#' @name gsBinCRT
# gsBinCRT function [sinew] ----
genBinCRT <- function(m = 1, n = 1, p_vec = c(0.5, 0.5),
                      rho = NULL, var_b = NULL, var_w = NULL) {
  # Specify between and within-cluster variances
  sigma_b_1 <- sqrt((rho / (1 - rho)) * (pi^2 / 3))
  sigma_b_2 <- sqrt((rho / (1 - rho)) * (pi^2 / 3))

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
    o1 <- log(p_vec[1] / (1 - p_vec[1])) + re1
    p1 <- exp(o1) / (1 + exp(o1))
    r1[(mi - 1) * n + 1:n] <- rbinom(n, 1, p1)

    re2 <- rnorm(1, mean = 0, sd = sigma_b_2)
    o2 <- log(p_vec[2] / (1 - p_vec[2])) + re2
    p2 <- exp(o2) / (1 + exp(o2))
    r2[(mi - 1) * n + 1:n] <- rbinom(n, 1, p2)
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

# iFrac function [sinew] ----
iFrac <- function(m_i, n_i, m_max, n_max, rho) {
  nif <- m_i * n_i * (1 + (n_max - 1) * rho)
  dif <- m_max * n_max * (1 + (n_i - 1) * rho)
  return(nif / dif)
}

# seBinDiff function [sinew] ----
seBinDiff <- function(x1_i, x2_i, m, n, rho, p_vec) {
  if (is.null(p_vec)) {
    p <- (sum(x1_i) + sum(x2_i)) / (length(x1_i) + length(x2_i))
  } else {
    p <- mean(p_vec)
  }

  if ((p * (1 - p)) == 0) {
    return(1)
  } else {
    se <- sqrt((2 * p * (1 - p)) * (1 + (n - 1) * rho) / (m * n))
    return(se)
  }
}