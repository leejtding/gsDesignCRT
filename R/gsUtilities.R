# checkLengths roxy [sinew] ----
#' @rdname checkScalar
#' @export
# checkLengths function [sinew] ----
checkLengths <- function(..., allowSingle = FALSE) {
  lens <- unlist(lapply(list(...), length))

  if (allowSingle) {
    lens <- lens[lens > 1]
  }

  if (length(lens) > 0 && length(unique(lens)) != 1) {
    parent <- as.character(sys.call(-1)[[1]])
    stop(
      if (length(parent) > 0) paste("In function", parent, ":") else "",
      "lengths of inputs are not all equal"
    )
  }

  invisible(NULL)
}


# checkRange roxy [sinew] ----
#' @rdname checkScalar
#' @export
# checkRange function [sinew] ----
checkRange <- function(x, interval = 0:1, inclusion = c(TRUE, TRUE),
                       varname = deparse(substitute(x)), tol = 0) {
  # check inputs
  checkVector(interval, "numeric")
  if (length(interval) != 2) {
    stop("Interval input must contain two elements")
  }

  interval <- sort(interval)
  checkVector(inclusion, "logical")
  inclusion <- if (length(inclusion) == 1) rep(inclusion, 2) else inclusion[1:2]

  xrange <- range(x)
  left <- ifelse(inclusion[1],
                 xrange[1] >= interval[1] - tol,
                 xrange[1] > interval[1] - tol)
  right <- ifelse(inclusion[2],
                  xrange[2] <= interval[2] + tol,
                  xrange[2] < interval[2] + tol)

  if (!(left && right)) {
    stop(paste(varname, " not on interval ",
               if (inclusion[1]) "[" else "(", interval[1], ", ",
               interval[2], if (inclusion[2]) "]" else ")",
               sep = ""), call. = TRUE)
  }

  invisible(NULL)
}


# checkScalar roxy [sinew] ----
#' @title Utility functions to verify variable properties
#'
#' @description Utility functions to verify an objects's properties including
#' whether it is a scalar or vector, the class, the length, and (if numeric)
#' whether the range of values is on a specified interval. Additionally, the
#' \code{checkLengths} function can be used to ensure that all the supplied
#' inputs have equal lengths.
#'
#' \code{isInteger} is similar to \code{\link{is.integer}} except that
#' \code{isInteger(1)} returns \code{TRUE} whereas \code{is.integer(1)} returns
#' \code{FALSE}.
#'
#' \code{checkScalar} is used to verify that the input object is a scalar as
#' well as the other properties specified above.
#'
#' \code{checkVector} is used to verify that the input object is an atomic
#' vector as well as the other properties as defined above.
#'
#' \code{checkRange} is used to check whether the numeric input object's values
#' reside on the specified interval.  If any of the values are outside the
#' specified interval, a \code{FALSE} is returned.
#'
#' \code{checkLength} is used to check whether all of the supplied inputs have
#' equal lengths.
#'
#' @param x any object.
#' @param isType character string defining the class that the input object is
#' expected to be.
#' @param length integer specifying the expected length of the object in the
#' case it is a vector. If \code{length=NULL}, the default, then no length
#' check is performed.
#' @param interval two-element numeric vector defining the interval over which
#' the input object is expected to be contained.  Use the \code{inclusion}
#' argument to define the boundary behavior.
#' @param inclusion two-element logical vector defining the boundary behavior
#' of the specified interval. A \code{TRUE} value denotes inclusion of the
#' corresponding boundary. For example, if \code{interval=c(3,6)} and
#' \code{inclusion=c(FALSE,TRUE)}, then all the values of the input object are
#' verified to be on the interval (3,6].
#' @param varname character string defining the name of the input variable as
#' sent into the function by the caller.  This is used primarily as a mechanism
#' to specify the name of the variable being tested when \code{checkRange} is
#' being called within a function.
#' @param tol numeric scalar defining the tolerance to use in testing the
#' intervals of the
#'
#' \code{\link{checkRange}} function.
#' @param \dots For the \code{\link{checkScalar}} and \code{\link{checkVector}}
#' functions, this input represents additional arguments sent directly to the
#' \code{\link{checkRange}} function. For the
#'
#' \code{\link{checkLengths}} function, this input represents the arguments to
#' check for equal lengths.
#' @param allowSingle logical flag. If \code{TRUE}, arguments that are vectors
#' comprised of a single element are not included in the comparative length
#' test in the \code{\link{checkLengths}} function. Partial matching on the
#' name of this argument is not performed so you must specify 'allowSingle' in
#' its entirety in the call.
#' @return
#' \code{isInteger}: Boolean value as checking result
#' Other functions have no return value, called for side effects
#' @examples
#'
#' # check whether input is an integer
#' isInteger(1)
#' isInteger(1:5)
#' try(isInteger("abc")) # expect error
#'
#' # check whether input is an integer scalar
#' checkScalar(3, "integer")
#'
#' # check whether input is an integer scalar that resides
#' # on the interval on [3, 6]. Then test for interval (3, 6].
#' checkScalar(3, "integer", c(3, 6))
#' try(checkScalar(3, "integer", c(3, 6), c(FALSE, TRUE))) # expect error
#'
#' # check whether the input is an atomic vector of class numeric,
#' # of length 3, and whose value all reside on the interval [1, 10)
#' x <- c(3, pi, exp(1))
#' checkVector(x, "numeric", c(1, 10), c(TRUE, FALSE), length = 3)
#'
#' # do the same but change the expected length; expect error
#' try(checkVector(x, "numeric", c(1, 10), c(TRUE, FALSE), length = 2))
#'
#' # create faux function to check input variable
#' foo <- function(moo) checkVector(moo, "character")
#' foo(letters)
#' try(foo(1:5)) # expect error with function and argument name in message
#'
#' # check for equal lengths of various inputs
#' checkLengths(1:2, 2:3, 3:4)
#' try(checkLengths(1, 2, 3, 4:5)) # expect error
#'
#' # check for equal length inputs but ignore single element vectors
#' checkLengths(1, 2, 3, 4:5, 7:8, allowSingle = TRUE)
#'
#'
#' @aliases checkLengths isInteger
#' @keywords programming
#' @rdname checkScalar
#' @export
#' @importFrom methods is
# checkScalar function [sinew] ----
checkScalar <- function(x, isType = "numeric", ...) {
  # check inputs
  if (!is.character(isType)) {
    stop("isType must be an object of class character")
  }

  # check scalar type
  if (isType == "integer") {
    bad <- (!isInteger(x) || length(x) > 1)
  } else {
    bad <- (!methods::is(c(x), isType) || length(x) > 1)
  }
  if (bad) {
    # create error message
    parent <- as.character(sys.call(-1)[[1]])
    varstr <- paste(if (length(parent) > 0) paste("In function", parent, ": variable") else "", deparse(substitute(x)))
    stop(varstr, " must be scalar of class ", isType)
  }

  # check if input is on specified interval
  if (length(list(...)) > 0) {
    checkRange(x, ..., varname = deparse(substitute(x)))
  }

  invisible(NULL)
}


# checkVector roxy [sinew] ----
#' @rdname checkScalar
#' @export
#' @importFrom methods is
# checkVector function [sinew] ----
checkVector <- function(x, isType = "numeric", ..., length = NULL) {
  # check inputs
  checkScalar(isType, "character")
  if (!is.null(length)) {
    checkScalar(length, "integer")
  }

  # define local functions
  "isVectorAtomic" <- function(x)
    return(is.atomic(x) & any(c(NROW(x), NCOL(x)) == 1))

  # check vector type
  bad <- if (isType == "integer") {
    !isVectorAtomic(x) || !isInteger(x)
  } else {
    # wrap "x" in c() to strip dimension(s)
    !isVectorAtomic(x) || !methods::is(c(x), isType)
  }
  if (bad) {
    # create error message
    parent <- as.character(sys.call(-1)[[1]])
    varstr <- paste(if (length(parent) > 0) paste("In function", parent, ": variable") else "", deparse(substitute(x)))
    stop(paste(varstr, " must be vector of class ", isType))
  }
  # check vector length
  if (!is.null(length) && (length(x) != length)) {
    stop(paste(varstr, " is a vector of length ", length(x), " but should be of length", length))
  }

  # check if input is on specified interval
  if (length(list(...)) > 0) {
    checkRange(x, ..., varname = deparse(substitute(x)))
  }

  invisible(NULL)
}


# isInteger roxy [sinew] ----
#' @rdname checkScalar
#' @export
# isInteger function [sinew] ----
isInteger <- function(x) all(is.numeric(x)) && all(round(x, 0) == x)


# checkMatrix function [sinew] ----
checkMatrix <- function(x, isType = "numeric", ..., nrows = NULL, ncols = NULL) {
  # check inputs
  checkScalar(isType, "character")
  if (!is.null(nrows)) {
    checkScalar(nrows, "integer")
  }
  if (!is.null(ncols)) {
    checkScalar(ncols, "integer")
  }

  # define local functions
  "isMatrixAtomic" <- function(x)
    return(is.atomic(x) & all(c(NROW(x), NCOL(x)) > 0))

  # check matrix type
  bad <- if (isType == "integer") {
    !isMatrixAtomic(x) || !isInteger(x)
  } else {
    # wrap "x" in c() to strip dimension(s)
    !isMatrixAtomic(x) || !methods::is(c(x), isType)
  }
  if (bad) {
    # create error message
    parent <- as.character(sys.call(-1)[[1]])
    varstr <- paste(if (length(parent) > 0) paste("In function", parent, ": variable") else "", deparse(substitute(x)))
    stop(paste(varstr, " must be matrix of class ", isType))
  }
  # check matrix dimensions
  if (!is.null(nrows) && (NROW(x) != nrows)) {
    stop(paste(varstr, "is a matrix with", NROW(x), "rows, but should have", nrows, "rows"))
  }
  if (!is.null(ncols) && (NCOL(x) != ncols)) {
    stop(paste(varstr, "is a matrix with", NCOL(x), "columns, but should have", ncols, "columns"))
  }

  # check if input is on specified interval
  if (length(list(...)) > 0) {
    checkRange(x, ..., varname = deparse(substitute(x)))
  }

  invisible(NULL)
}

# expandTwo function [sinew] ----
expandTwo <- function(x) {
  if (length(x) == 1) {
    rep(x, 2L)
  } else if (length(x) == 2) {
    x
  } else {
    stop("Must be length 1 or 2.")
  }
}

# checkAlloc function [sinew] ----
checkAlloc <- function(m_alloc) {
  m_alloc <- as.numeric(m_alloc)
  if (length(m_alloc) != 2) {
    stop("m_alloc must be a vector of length 2.")
  }
  if (any(!is.finite(m_alloc)) || any(m_alloc < 0)) {
    stop("m_alloc must be finite and non-negative.")
  }
  if (sum(m_alloc) != 1) {
    stop("m_alloc must sum up to 1.")
  }
  invisible(NULL)
}

# checkSF function [sinew] ----
checkSF <- function(sf) {
  if (is.character(sf) &&
        !is.element(sf, c("OF", "Pocock", "WT"))) {
    stop("Character specification of upper spending may only be WT, OF or
         Pocock")
  } else if (!is.function(sf)) {
    stop("Upper spending function mis-specified")
  }
}