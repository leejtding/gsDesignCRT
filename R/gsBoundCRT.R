# gsUpperCRT roxy [sinew] ----
#' @title Boundary derivation for efficacy stopping only.
#'
#' @description  \code{gsUpperCRT} is used to calculate the stopping boundaries
#' for a group sequential trial with 1 or 2-sided efficacy-only stopping. Code
#' adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#' boundaries are calculated. At any given interim analysis this should include
#' the information levels at the current and previous analyses.
#' @param a Lower futility boundary values for the interim analyses at the
#' information levels specified by I.
#' @param falsepos Type I error spent for the interim analyses at the
#' information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{As input.}
#'  \item{b}{Computed upper efficacy boundaries at the specified interim
#'  analyses.} \item{r}{As input.} \item{error}{Error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gsupper1
#' @useDynLib gsDesignCRT gsupper2
#'
#' @rdname gsUpperCRT
# gsUpperCRT function [sinew] ----
gsUpperCRT <- function(theta = 0, I, a = NULL, falsepos, sides = 1,
                       tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  k <- as.integer(length(I))
  if (is.null(a)) {
    if (sides == 1) {
      a <- rep(-10, k)
    } else {
      a <- rep(0, k)
    }
  }

  checkScalar(sides, "integer", c(1, 2))
  if (sides == 2) {
    checkScalar(theta, "numeric", c(0, Inf))
    checkVector(a, "numeric", c(0, Inf))
  } else {
    checkScalar(theta, "numeric")
    checkVector(a, "numeric")
  }
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(falsepos, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(a, falsepos, I)

  # Set up numerical calculations
  if (falsepos[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)

  storage.mode(theta) <- "double"
  storage.mode(I) <- "double"
  storage.mode(a) <- "double"
  storage.mode(falsepos) <- "double"
  storage.mode(tol) <- "double"

  problo <- a
  b <- a
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    xx <- .C("gsupper1",
             k, theta, I, a, b, problo, falsepos, tol, r, retval, printerr)
  } else {
    xx <- .C("gsupper2",
             k, theta, I, a, b, problo, falsepos, tol, r, retval, printerr)
  }

  y <- list(k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
            falseneg = xx[[6]], falsepos = xx[[7]], tol = xx[[8]], r = xx[[9]],
            error = xx[[10]])

  # Set last interim analysis boundaries equal if needed
  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
  }
  y
}

# gsLowerCRT roxy [sinew] ----
#' @title Boundary derivation for binding or non-binding futility stopping only.
#'
#' @description  \code{gsLowerCRT} is used to calculate the stopping boundaries
#'  for a group sequential trial with 1 or 2-sided futility-only stopping. Code
#'  adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#'  boundaries are calculated. At any given interim analysis this should include
#'  the information levels at the current and previous analyses.
#' @param falseneg Type II error spent for the interim analyses at the
#'  information levels specified by I.
#' @param b Upper efficacy boundary values for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param binding \code{TRUE=} binding (default) \cr \code{FALSE=} non-binding
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{Computed
#'  lower futility boundaries at the specified interim analyses.} \item{b}{As
#'  input.} \item{r}{As input.} \item{error}{error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gslower1
#' @useDynLib gsDesignCRT gslower2
#'
#' @rdname gsLowerCRT
# gsLowerCRT function [sinew] ----
gsLowerCRT <- function(theta = 0, I, falseneg, b = NULL, sides = 1,
                       binding = TRUE, tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  k <- as.integer(length(I))
  if (is.null(b)) {
    b <- rep(10, k)
  }
  checkScalar(sides, "integer", c(1, 2))
  if (sides == 2) {
    checkScalar(theta, "numeric", c(0, Inf))
    checkVector(b, "numeric", c(0, Inf))
  } else {
    checkScalar(theta, "numeric")
    checkVector(b, "numeric")
  }
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(falseneg, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(as.numeric(binding), "integer", c(0, 1))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(b, falseneg, I)

  # Set up numerical calculations
  if (falseneg[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)

  storage.mode(theta) <- "double"
  storage.mode(I) <- "double"
  storage.mode(b) <- "double"
  storage.mode(falseneg) <- "double"
  storage.mode(tol) <- "double"

  probhi <- b
  a <- b
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    xx <- .C("gslower1",
             k, theta, I, a, b, falseneg, probhi, tol, r, retval, printerr)
  } else {
    xx <- .C("gslower2",
             k, theta, I, a, b, falseneg, probhi, tol, r, retval, printerr)
  }

  y <- list(k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
            falseneg = xx[[6]], falsepos = xx[[7]], tol = xx[[8]], r = xx[[9]],
            error = xx[[10]])

  # Set last interim analysis boundaries equal if needed
  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
  }
  y
}

# gsBoundsCRT roxy [sinew] ----
#' @title Boundary derivation for efficacy and binding or non-binding futility
#'  stopping.
#'
#' @description \code{gsBoundsCRT} is used to calculate the stopping boundaries
#'  for a group sequential trial with 1 or 2-sided efficacy and binding or
#'  non-binding futility stopping. Code adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#'  boundaries are calculated. At any given interim analysis this should include
#'  the information levels at the current and previous analyses.
#' @param falseneg Type II error spent for the interim analyses at the
#'  information levels specified by I.
#' @param falsepos Type I error spent for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param binding \code{TRUE=} binding (default) \cr \code{FALSE=} non-binding
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{Computed
#'  lower futility boundaries at the specified interim analyses.}
#'  \item{b}{Computed upper efficacy boundaries at the specified interim
#'  analyses.} \item{r}{As input.} \item{error}{error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gsbounds1
#' @useDynLib gsDesignCRT gsboundsnb1
#' @useDynLib gsDesignCRT gsbounds2
#' @useDynLib gsDesignCRT gsboundsnb2
#'
#' @rdname gsBoundsCRT
# gsBoundsCRT function [sinew] ----
gsBoundsCRT <- function(theta = 0, I, falseneg, falsepos, sides = 1,
                        binding = TRUE, tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  checkScalar(theta, "numeric", c(0, Inf))
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(falseneg, "numeric", c(0, 1), c(FALSE, FALSE))
  checkVector(falsepos, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(sides, "integer", c(1, 2))
  checkScalar(as.numeric(binding), "integer", c(0, 1))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(falseneg, falsepos, I)

  # Set up numerical calculations
  k <- as.integer(length(I))
  if (falseneg[k] < 0.) stop("Final futility spend must be >= 0")
  if (falsepos[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)
  storage.mode(I) <- "double"
  storage.mode(falseneg) <- "double"
  storage.mode(falsepos) <- "double"
  storage.mode(tol) <- "double"
  a <- falseneg
  b <- falsepos
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    if (binding) {
      xx <- .C("gsbounds1",
               k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
    } else {
      xx <- .C("gsboundsnb1",
               k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
    }
  } else {
    if (binding) {
      xx <- .C("gsbounds2",
               k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
    } else {
      xx <- .C("gsboundsnb2",
               k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
    }
  }

  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
    falseneg = xx[[7]], falsepos = xx[[6]], tol = xx[[8]], r = xx[[9]],
    error = xx[[10]], diff = xx[[5]] - xx[[4]]
  )
  y
}

# gsProbabilityCRT roxy [sinew] ----
#' @title Compute stopping boundary crossing probabilities.
#'
#' @param theta Effect size. Default value is 0.
#' @param I Information levels at interim analyses for which boundary crossing
#'  probabilities are computed.
#' @param a Lower futility boundaries for the interim analyses at the
#'  information levels specified by I.
#' @param b Upper efficacy boundaries for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{lower}{List
#'  containing the lower futility boundaries (\code{bound}) and the
#'  corresponding probabilities of crossing the lower boundaries given
#'  \code{theta} (\code{prob}).} \item{upper}{List containing the upper efficacy
#'  boundaries (\code{bound}) and the corresponding probabilities of crossing
#'  the upper boundaries given theta (\code{prob}).} \item{eI}{Expected
#'  information given \code{theta}} \item{r}{As input.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @rdname gsProbabilityCRT
#' @useDynLib gsDesignCRT probrej1
#' @useDynLib gsDesignCRT probrej2
#'
# gsProbabilityCRT function [sinew] ----
gsProbabilityCRT <- function(theta = 0, I = 1, a = 0, b = 1,
                             sides = 1, r = 18) {
  # Check input arguments
  checkVector(theta, "numeric")

  # Check remaining input arguments
  checkVector(I - c(0, I[1:length(I) - 1]),
              "numeric", c(0, Inf), c(FALSE, FALSE))
  checkScalar(sides, "integer", c(1, 2))
  checkScalar(r, "integer", c(1, 80))
  checkLengths(I, a, b)

  # Cast integer scalars
  k <- as.integer(length(I))
  ntheta <- as.integer(length(theta))
  r <- as.integer(r)

  phi <- as.double(c(1:(k * ntheta)))
  plo <- as.double(c(1:(k * ntheta)))

  # Compute crossing probabilities
  if (sides == 1) {
    xx <- .C("probrej1",
             k, ntheta, as.double(theta), as.double(I),
             as.double(a), as.double(b), plo, phi, r)
  } else if (sides == 2) {
    xx <- .C("probrej2",
             k, ntheta, as.double(theta), as.double(I),
             as.double(a), as.double(b), plo, phi, r)
  } else {
    stop("invalid test type")
  }
  plo <- matrix(xx[[7]], k, ntheta)
  phi <- matrix(xx[[8]], k, ntheta)
  powr <- as.vector(rep(1, k) %*% phi)
  futile <- rep(1, k) %*% plo

  # Compute expected information under null and alternative hypotheses
  if (k == 1) {
    eI <- (as.vector(I * (plo + phi)) +
             as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile)))
  } else {
    eI <- (as.vector(I %*% (plo + phi)) +
             as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile)))
  }

  x <- list(k = xx[[1]], theta = xx[[3]], I = xx[[4]],
            lower = list(bound = xx[[5]], prob = plo),
            upper = list(bound = xx[[6]], prob = phi), eI = eI, r = r)

  x
}

###
# Hidden Functions
###

# gsbetadiffCRT function [sinew] ----
gsbetadiffCRT <- function(Imax, theta, beta, time, a, b, sides, r = 18) {
  # Compute difference between actual and desired Type II error
  I <- time * Imax
  x <- gsprobCRT(theta = theta, I = I, a = a, b = b, sides = sides, r = r)
  (sum(x$probhi) + sum(x$problo)) - (1 - beta)
}

# gsalphadiffCRT function [sinew] ----
gsalphadiffCRT <- function(Imax, theta, alpha, falseneg, time, b, sides,
                           binding, tol = 0.000001, r = 18) {
  # Compute difference between actual and desired Type I error
  I <- time * Imax
  lower <- gsLowerCRT(theta = theta, I = I, falseneg = falseneg, b = b,
                      sides = sides, binding = binding, tol = tol, r = r)
  x <- gsprobCRT(theta = 0, I = I, a = lower$a, b = b, sides = sides, r = r)
  (sum(x$probhi) + sum(x$problo)) - alpha
}

# gsbounddiffCRT function [sinew] ----
gsbounddiffCRT <- function(Imax, theta, falseneg, falsepos, time, sides,
                           binding, tol = 0.000001, r = 18) {
  # Compute difference between last upper and lower stopping boundaries
  I <- time * Imax
  bounds <- gsBoundsCRT(theta = theta, I = I, falseneg = falseneg,
                        falsepos = falsepos, sides = sides, binding = binding,
                        tol = tol, r = r)
  bounds$diff[length(time)]
}

# gsprobCRT function [sinew] ----
gsprobCRT <- function(theta, I, a, b, sides = 1, r = 18) {
  nanal <- as.integer(length(I))
  ntheta <- as.integer(length(theta))
  phi <- as.double(c(1:(nanal * ntheta)))
  plo <- as.double(c(1:(nanal * ntheta)))
  if (sides == 1) {
    xx <- .C(
      "probrej1", nanal, ntheta, as.double(theta), as.double(I),
      as.double(a), as.double(b), plo, phi, as.integer(r)
    )
  } else if (sides == 2) {
    xx <- .C(
      "probrej2", nanal, ntheta, as.double(theta), as.double(I),
      as.double(a), as.double(b), plo, phi, as.integer(r)
    )
  } else {
    stop("invalid number of sides for test")
  }
  plo <- matrix(xx[[7]], nanal, ntheta)
  phi <- matrix(xx[[8]], nanal, ntheta)
  powr <- rep(1, nanal) %*% phi
  futile <- rep(1, nanal) %*% plo
  if (nanal == 1) {
    eI <- as.vector(I * (plo + phi)) + as.vector(I[nanal] * (t(rep(1, ntheta))
                                                             - powr - futile))
  } else {
    eI <- as.vector(I %*% (plo + phi)) + as.vector(I[nanal] * (t(rep(1, ntheta))
                                                               - powr - futile))
  }
  list(
    k = xx[[1]], theta = xx[[3]], I = xx[[4]], a = xx[[5]], b = xx[[6]],
    problo = plo, probhi = phi, powr = powr, eI = eI, r = r
  )
}