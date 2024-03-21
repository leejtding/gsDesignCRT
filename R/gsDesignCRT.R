# gsUpper1 roxy [sinew] ----
#' @title Boundary derivation - upper boundary
#' @export
#' @useDynLib gsDesignCRT gsupper1
#' @rdname gsUpper1
# gsUpper1 function [sinew] ----
gsUpper1 <- function(theta, I, a, falsepos, tol = 0.000001, r = 42, printerr = 0) {
  # check input arguments
  checkScalar(theta, "numeric")
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(a, "numeric")
  checkVector(falsepos, "numeric", c(0, 1), c(TRUE, FALSE))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(a, falsepos, I)

  # coerce type
  k <- as.integer(length(I))
  if (falsepos[k] <= 0.) stop("Final spend must be > 0")
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

  xx <- .C("gsupper1", k, theta, I, a, b, problo, falsepos, tol, r, retval, printerr)

  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
    problo = xx[[6]], falsepos = xx[[7]], tol = xx[[8]], r = xx[[9]], error = xx[[10]]
  )

  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
    z <- gsprob(theta = theta, I = I, a = a, b = y$b, r = r)
    y$falsepos <- z$falsepos
    y$problo <- z$problo
  }
  y
}

# gsUpper2 roxy [sinew] ----
#' @title Boundary derivation - upper boundary
#' @export
#' @useDynLib gsDesignCRT gsupper2
#' @rdname gsUpper2
# gsUpper2 function [sinew] ----
gsUpper2 <- function(theta, I, a, falsepos, tol = 0.000001, r = 42, printerr = 0) {
  # check input arguments
  checkScalar(theta, "numeric")
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(a, "numeric")
  checkVector(falsepos, "numeric", c(0, 1), c(TRUE, FALSE))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(a, falsepos, I)

  # coerce type
  k <- as.integer(length(I))
  if (falsepos[k] <= 0.) stop("Final spend must be > 0")
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

  xx <- .C("gsupper2", k, theta, I, a, b, problo, falsepos, tol, r, retval, printerr)

  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
    problo = xx[[6]], falsepos = xx[[7]], tol = xx[[8]], r = xx[[9]], error = xx[[10]]
  )

  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
    z <- gsprob(theta = theta, I = I, a = a, b = y$b, r = r)
    y$falsepos <- z$falsepos
    y$problo <- z$problo
  }
  y
}

# gsBounds1 roxy [sinew] ----
#' @title Boundaries for 1-sided test with efficacy and futility stopping
#' @export
#' @useDynLib gsDesignCRT gsbounds1
#' @useDynLib gsDesignCRT gsboundsnb1
#' @rdname gsBounds1
# gsBounds1 function [sinew] ----
gsBounds1 <- function(theta, I, falseneg, falsepos, binding = TRUE, tol = 0.000001, r = 42, printerr = 0) {
  # check input arguments
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(falseneg, "numeric", c(0, 1), c(TRUE, FALSE))
  checkVector(falsepos, "numeric", c(0, 1), c(TRUE, FALSE))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(falseneg, falsepos, I)

  k <- as.integer(length(I))
  if (falseneg[k] <= 0.) stop("Final futility spend must be > 0")
  if (falsepos[k] <= 0.) stop("Final efficacy spend must be > 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)
  storage.mode(I) <- "double"
  storage.mode(falseneg) <- "double"
  storage.mode(falsepos) <- "double"
  storage.mode(tol) <- "double"
  a <- falseneg
  b <- falsepos
  retval <- as.integer(0)
  if (binding) {
    xx <- .C("gsbounds1", k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
  } else {
    xx <- .C("gsboundsnb1", k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
  }
  rates <- list(falsepos = xx[[6]], falseneg = xx[[7]])
  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]], rates = rates, tol = xx[[8]],
    r = xx[[9]], error = xx[[10]], diff = xx[[5]] - xx[[4]]
  )

  y
}

# gsBounds2 roxy [sinew] ----
#' @title Boundaries for 2-sided test with efficacy and futility stopping
#' @export
#' @useDynLib gsDesignCRT gsbounds2
#' @useDynLib gsDesignCRT gsboundsnb2
#' @rdname gsBounds2
# gsBounds2 function [sinew] ----
gsBounds2 <- function(theta, I, falseneg, falsepos, binding = TRUE, tol = 0.000001, r = 42, printerr = 0) {
  # check input arguments
  checkVector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkVector(falseneg, "numeric", c(0, 1), c(TRUE, FALSE))
  checkVector(falsepos, "numeric", c(0, 1), c(TRUE, FALSE))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkScalar(r, "integer", c(1, 80))
  checkScalar(printerr, "integer")
  checkLengths(falseneg, falsepos, I)

  k <- as.integer(length(I))
  if (falseneg[k] <= 0.) stop("Final futility spend must be > 0")
  if (falsepos[k] <= 0.) stop("Final efficacy spend must be > 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)
  storage.mode(I) <- "double"
  storage.mode(falseneg) <- "double"
  storage.mode(falsepos) <- "double"
  storage.mode(tol) <- "double"
  a <- falseneg
  b <- falsepos
  retval <- as.integer(0)
  if (binding) {
    xx <- .C("gsbounds2", k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
  } else {
    xx <- .C("gsboundsnb2", k, theta, I, a, b, falseneg, falsepos, tol, r, retval, printerr)
  }
  rates <- list(falsepos = xx[[6]], falseneg = xx[[7]])
  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]], rates = rates, tol = xx[[8]],
    r = xx[[9]], error = xx[[10]], diff = xx[[5]] - xx[[4]]
  )

  y
}

# gsProbability1 roxy [sinew] ----
#' @title Boundary Crossing Probabilities for 1-sided test
#' @export
#' @rdname gsProbability1
#' @useDynLib gsDesignCRT probrej1
# gsProbability function [sinew] ----
gsProbability1 <- function(k = 0, theta, I, a, b, r = 42, d = NULL, overrun = 0) {
  # compute boundary crossing probabilities and return in a gsProbability structure

  # check input arguments
  checkScalar(k, "integer", c(0, 30))
  checkVector(theta, "numeric")

  if (k == 0) {
    if (!inherits(d, "gsDesignCRT")) {
      stop("d should be an object of class gsDesignCRT")
    }
    return(gsDProb(theta = theta, d = d))
  }

  # check remaining input arguments
  checkVector(x = overrun, isType = "numeric", interval = c(0, Inf),
              inclusion = c(TRUE, FALSE))
  if (length(overrun) != 1 && length(overrun) != k - 1) {
    stop(paste("overrun length should be 1 or ", as.character(k - 1)))
  }
  checkScalar(r, "integer", c(1, 80))
  checkLengths(I, a, b)
  if (k != length(a)) {
    stop("Lengths of I, a, and b must all equal k")
  }

  # cast integer scalars
  ntheta <- as.integer(length(theta))
  k <- as.integer(k)
  r <- as.integer(r)

  phi <- as.double(c(1:(k * ntheta)))
  plo <- as.double(c(1:(k * ntheta)))
  xx <- .C(
    "probrej1", k, ntheta, as.double(theta), as.double(I),
    as.double(a), as.double(b), plo, phi, r
  )
  plo <- matrix(xx[[7]], k, ntheta)
  phi <- matrix(xx[[8]], k, ntheta)
  powr <- as.vector(rep(1, k) %*% phi)
  futile <- rep(1, k) %*% plo
  if (k == 1) {
    IOver <- I[k]
  } else {
    IOver <- c(I[1:(k - 1)] + overrun, I[k])
  }
  IOver[IOver > I[k]] <- I[k]
  if (k == 1) {
    eI <- as.vector(IOver * (plo + phi)) + as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile))
  } else {
    eI <- as.vector(IOver %*% (plo + phi)) + as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile))
  }
  x <- list(
    k = xx[[1]], theta = xx[[3]], I = xx[[4]], lower = list(bound = xx[[5]], prob = plo),
    upper = list(bound = xx[[6]], prob = phi), eI = eI, r = r, overrun = overrun
  )

  class(x) <- "gsProbability"

  x
}

# gsProbability2 roxy [sinew] ----
#' @title Boundary Crossing Probabilities for 2-sided test
#' @export
#' @useDynLib gsDesignCRT probrej2
# gsProbability function [sinew] ----
gsProbability2 <- function(k = 0, theta, I, a, b, r = 42, d = NULL, overrun = 0) {
  # compute boundary crossing probabilities and return in a gsProbability structure

  # check input arguments
  checkScalar(k, "integer", c(0, 30))
  checkVector(theta, "numeric")

  if (k == 0) {
    if (!inherits(d, "gsDesignCRT")) {
      stop("d should be an object of class gsDesignCRT")
    }
    return(gsDProb(theta = theta, d = d))
  }

  # check remaining input arguments
  checkVector(x = overrun, isType = "numeric", interval = c(0, Inf),
              inclusion = c(TRUE, FALSE))
  if (length(overrun) != 1 && length(overrun) != k - 1) {
    stop(paste("overrun length should be 1 or ", as.character(k - 1)))
  }
  checkScalar(r, "integer", c(1, 80))
  checkLengths(I, a, b)
  if (k != length(a)) {
    stop("Lengths of I, a, and b must all equal k")
  }

  # cast integer scalars
  ntheta <- as.integer(length(theta))
  k <- as.integer(k)
  r <- as.integer(r)

  phi <- as.double(c(1:(k * ntheta)))
  plo <- as.double(c(1:(k * ntheta)))
  xx <- .C(
    "probrej2", k, ntheta, as.double(theta), as.double(I),
    as.double(a), as.double(b), plo, phi, r
  )
  plo <- matrix(xx[[7]], k, ntheta)
  phi <- matrix(xx[[8]], k, ntheta)
  powr <- as.vector(rep(1, k) %*% phi)
  futile <- rep(1, k) %*% plo
  if (k == 1) {
    IOver <- I[k]
  } else {
    IOver <- c(I[1:(k - 1)] + overrun, I[k])
  }
  IOver[IOver > I[k]] <- I[k]
  if (k == 1) {
    eI <- as.vector(IOver * (plo + phi)) + as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile))
  } else {
    eI <- as.vector(IOver %*% (plo + phi)) + as.vector(I[k] * (t(rep(1, ntheta)) - powr - futile))
  }
  x <- list(
    k = xx[[1]], theta = xx[[3]], I = xx[[4]], lower = list(bound = xx[[5]], prob = plo),
    upper = list(bound = xx[[6]], prob = phi), eI = eI, r = r, overrun = overrun
  )

  class(x) <- "gsProbability"

  x
}

###
# Hidden Functions
###

# gsbetadiff function [sinew] ----
gsbetadiff <- function(Imax, theta, beta, time, a, b, tol = 0.000001, r = 42) {
  # compute difference between actual and desired Type II error
  I <- time * Imax
  x <- gsprob1(theta, I, a, b, r)
  sum(x$probhi) - (1 - beta)
}

# gsbetadiff2 function [sinew] ----
gsbetadiff2 <- function(Imax, theta, beta, time, a, b, tol = 0.000001, r = 42) {
  # compute difference between actual and desired Type II error
  I <- time * Imax
  x <- gsprob2(theta, I, a, b, r)
  sum(x$probhi) - (1 - beta)
}

# gsbounddiff function [sinew] ----
gsbounddiff <- function(Imax, theta, falseneg, falsepos, time, binding, tol = 0.000001, r = 42) {
  I <- time * Imax
  bounds <- gsBounds1(theta, I, falseneg, falsepos, binding, tol, r)
  bounds$diff[length(time)]
}

# gsbounddiff2 function [sinew] ----
gsbounddiff2 <- function(Imax, theta, falseneg, falsepos, time, binding, tol = 0.000001, r = 42) {
  I <- time * Imax
  bounds <- gsBounds2(theta, I, falseneg, falsepos, binding, tol, r)
  return(bounds$diff[length(time)])
}

# gsprob1 function [sinew] ----
gsprob1 <- function(theta, I, a, b, r = 42) {
  nanal <- as.integer(length(I))
  ntheta <- as.integer(length(theta))
  phi <- as.double(c(1:(nanal * ntheta)))
  plo <- as.double(c(1:(nanal * ntheta)))
  xx <- .C(
    "probrej1", nanal, ntheta, as.double(theta), as.double(I),
    as.double(a), as.double(b), plo, phi, as.integer(r)
  )
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

# gsprob2 function [sinew] ----
gsprob2 <- function(theta, I, a, b, r = 42) {
  nanal <- as.integer(length(I))
  ntheta <- as.integer(length(theta))
  phi <- as.double(c(1:(nanal * ntheta)))
  plo <- as.double(c(1:(nanal * ntheta)))
  xx <- .C(
    "probrej2", nanal, ntheta, as.double(theta), as.double(I),
    as.double(a), as.double(b), plo, phi, as.integer(r)
  )
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
