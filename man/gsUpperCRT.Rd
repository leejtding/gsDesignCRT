% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsBoundCRT.R
\name{gsUpperCRT}
\alias{gsUpperCRT}
\title{Boundary derivation for efficacy stopping only.}
\usage{
gsUpperCRT(
  theta = 0,
  I,
  a = NULL,
  falsepos,
  sides = 1,
  tol = 1e-06,
  r = 18,
  printerr = 0
)
}
\arguments{
\item{theta}{Effect size under null hypothesis. Default value is 0.}

\item{I}{Information levels at interim analyses for which stopping
boundaries are calculated. At any given interim analysis this should include
the information levels at the current and previous analyses.}

\item{a}{Lower futility boundary values for the interim analyses at the
information levels specified by I.}

\item{falsepos}{Type I error spent for the interim analyses at the
information levels specified by I.}

\item{sides}{\code{1=} 1-sided test (default) \cr \code{2=} 2-sided test}

\item{tol}{Tolerance for error (default is 0.000001). Normally this will not
be changed by the user.  This does not translate directly to number of
digits of accuracy, so use extra decimal places.}

\item{r}{Integer value controlling grid for numerical integration as in
Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
values provide larger number of grid points and greater accuracy.  Normally
\code{r} will not be changed by the user.}

\item{printerr}{Print output for debugging.}
}
\value{
Object containing the following elements: \item{k}{Number of interim
 analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{As input.}
 \item{b}{Computed upper efficacy boundaries at the specified interim
 analyses.} \item{r}{As input.} \item{error}{Error flag returned; 0 if
 convergence; 1 indicates error.}
}
\description{
\code{gsUpperCRT} is used to calculate the stopping boundaries
for a group sequential trial with 1 or 2-sided efficacy-only stopping. Code
adapted from gsDesign package.
}
\references{
Jennison C and Turnbull BW (2000), \emph{Group Sequential
Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
}
\author{
Lee Ding \email{lee_ding@g.harvard.edu}
}
