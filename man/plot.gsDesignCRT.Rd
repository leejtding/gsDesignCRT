% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsqplot.R
\name{plot.gsDesignCRT}
\alias{plot.gsDesignCRT}
\alias{plot.gsProbability}
\title{Plots for group sequential designs}
\usage{
\method{plot}{gsDesignCRT}(x, plottype = 1, base = FALSE, ...)

\method{plot}{gsProbability}(x, plottype = 2, base = FALSE, ...)
}
\arguments{
\item{x}{Object of class \code{gsDesignCRT} for \code{plot.gsDesignCRT()} or
\code{gsProbability} for

\code{plot.gsProbability()}.}

\item{plottype}{1=boundary plot (default for \code{gsDesignCRT}),

2=power plot (default for \code{gsProbability}),

3=approximate treatment effect at boundaries,

4=conditional power at boundaries,

5=spending function plot (only available if \code{class(x)=="gsDesignCRT"}),

6=expected sample size plot, and

7=B-values at boundaries.

Character values for \code{plottype} may also be entered: \code{"Z"} for
plot type 1, \code{"power"} for plot type 2, \code{"thetahat"} for plot type
3, \code{"CP"} for plot type 4, \code{"sf"} for plot type 5, \code{"ASN"},
\code{"N"} or \code{"n"} for plot type 6, and \code{"B"}, \code{"B-val"} or
\code{"B-value"} for plot type 7.}

\item{base}{Default is FALSE, which means ggplot2 graphics are used. If
true, base graphics are used for plotting.}

\item{...}{This allows many optional arguments that are standard when
calling \code{plot}.

Other arguments include:

\code{theta} which is used for \code{plottype=2}, \code{4}, \code{6};
normally defaults will be adequate; see details.

\code{ses=TRUE} which applies only when \code{plottype=3} and

\code{class(x)=="gsDesignCRT"}; indicates that approximate standardized effect
size at the boundary is to be plotted rather than the approximate natural parameter.

\code{xval="Default"} which is only effective when \code{plottype=2} or
\code{6}. Appropriately scaled (reparameterized) values for x-axis for power
and expected sample size graphs; see details.}
}
\value{
An object of \code{class(x)}; in many cases this is the input value
of \code{x}, while in others \code{x$theta} is replaced and corresponding
characteristics computed; see details.
}
\description{
The \code{plot()} function has been extended to work with objects returned
by \code{gsDesignCRT()} and \code{gsProbability()}.  For objects of type
\code{gsDesignCRT}, seven types of plots are provided: z-values at boundaries
(default), power, approximate treatment effects at boundaries, conditional
power at boundaries, spending functions, expected sample size, and B-values
at boundaries. For objects of type \code{gsProbability} plots are available
for z-values at boundaries, power (default), approximate treatment effects at
boundaries, conditional power, expected sample size and B-values at
boundaries.

The intent is that many standard \code{plot()} parameters will function as
expected; exceptions to this rule exist. In particular, \code{main, xlab,
ylab, lty, col, lwd, type, pch, cex} have been tested and work for most
values of \code{plottype}; one exception is that \code{type="l"} cannot be
overridden when \code{plottype=2}. Default values for labels depend on
\code{plottype} and the class of \code{x}.

Note that there is some special behavior for values plotted and returned for
power and expected sample size (ASN) plots for a \code{gsDesignCRT} object. A
call to \code{x<-gsDesignCRT()} produces power and expected sample size for
only two \code{theta} values: 0 and \code{x$delta}.  The call \code{plot(x,
plottype="Power")} (or \code{plot(x,plottype="ASN"}) for a \code{gsDesignCRT}
object produces power (expected sample size) curves and returns a
\code{gsDesignCRT} object with \code{theta} values determined as follows.  If
\code{theta} is non-null on input, the input value(s) are used. Otherwise,
for a \code{gsProbability} object, the \code{theta} values from that object
are used. For a \code{gsDesignCRT} object where \code{theta} is input as
\code{NULL} (the default), \code{theta=seq(0,2,.05)*x$delta}) is used.  For
a \code{gsDesignCRT} object, the x-axis values are rescaled to
\code{theta/x$delta} and the label for the x-axis \eqn{theta / delta}. For a
\code{gsProbability} object, the values of \code{theta} are plotted and are
labeled as \eqn{theta}. See examples below.

Approximate treatment effects at boundaries are computed dividing the Z-values
at the boundaries by the square root of \code{n.I} at that analysis.

Spending functions are plotted for a continuous set of values from 0 to 1.
This option should not be used if a boundary is used or a pointwise spending
function is used (\code{sfu} or \code{sfl="WT", "OF", "Pocock"} or
\code{sfPoints}).

Conditional power is computed using the function \code{gsBoundCP()}.  The
default input for this routine is \code{theta="thetahat"} which will compute
the conditional power at each bound using the approximate treatment effect at
that bound.  Otherwise, if the input is \code{gsDesignCRT} object conditional
power is computed assuming \code{theta=x$delta}, the original effect size
for which the trial was planned.

Average sample number/expected sample size is computed using \code{n.I} at
each analysis times the probability of crossing a boundary at that analysis.
If no boundary is crossed at any analysis, this is counted as stopping at
the final analysis.

B-values are Z-values multiplied by \code{sqrt(t)=sqrt(x$n.I/x$n.I[x$k])}.
Thus, the expected value of a B-value at an analysis is the true value of
\eqn{theta} multiplied by the proportion of total planned observations at
that time. See Proschan, Lan and Wittes (2006).
}
\note{
The gsDesignCRT technical manual is available at
  \url{https://keaven.github.io/gsd-tech-manual/}.
}
\examples{
library(ggplot2)
#  symmetric, 2-sided design with O'Brien-Fleming-like boundaries
#  lower bound is non-binding (ignored in Type I error computation)
#  sample size is computed based on a fixed design requiring n=100
x <- gsDesignCRT(k = 5, test.type = 2, n.fix = 100)
x

# the following translate to calls to plot.gsDesignCRT since x was
# returned by gsDesignCRT; run these commands one at a time
plot(x)
plot(x, plottype = 2)
plot(x, plottype = 3)
plot(x, plottype = 4)
plot(x, plottype = 5)
plot(x, plottype = 6)
plot(x, plottype = 7)

#  choose different parameter values for power plot
#  start with design in x from above
y <- gsProbability(
  k = 5, theta = seq(0, .5, .025), x$n.I,
  x$lower$bound, x$upper$bound
)

# the following translates to a call to plot.gsProbability since
# y has that type
plot(y)
}
\references{
Jennison C and Turnbull BW (2000), \emph{Group Sequential
Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.

Proschan, MA, Lan, KKG, Wittes, JT (2006), \emph{Statistical Monitoring of
Clinical Trials. A Unified Approach}.  New York: Springer.
}
\seealso{
\code{\link{gsDesignCRT}}, \code{\link{gsProbability}}
}
\author{
Keaven Anderson \email{keaven_anderson@merck.com}
}
\keyword{design}
