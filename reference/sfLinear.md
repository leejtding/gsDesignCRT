# Piecewise Linear and Step Function Spending Functions

The function `sfLinear()` allows specification of a piecewise linear
spending function. The function `sfStep()` specifies a step function
spending function. Both functions provide complete flexibility in
setting spending at desired timepoints in a group sequential design.
Normally these function will be passed to
[`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
in the parameter `sfu` for the upper bound or `sfl` for the lower bound
to specify a spending function family for a design. When passed to
[`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md),
the value of `param` would be passed to `sfLinear()` or `sfStep()`
through the
[`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
arguments `sfupar` for the upper bound and `sflpar` for the lower bound.

Note that `sfStep()` allows setting a particular level of spending when
the timing is not strictly known; an example shows how this can inflate
Type I error when timing of analyses are changed based on knowing the
treatment effect at an interim.

## Usage

``` r
sfLinear(alpha, t, param)

sfStep(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size or information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size or information for which the
  spending function will be computed.

- param:

  A vector with a positive, even length. Values must range from 0 to 1,
  inclusive. Letting `m <- length(param/2)`, the first `m` points in
  param specify increasing values strictly between 0 and 1 corresponding
  to interim timing (proportion of final total statistical information).
  The last `m` points in `param` specify non-decreasing values from 0 to
  1, inclusive, with the cumulative proportion of spending at the
  specified timepoints.

## Value

An object of type `spendfn`. The cumulative spending returned in
`sfLinear$spend` is 0 for `t <= 0` and `alpha` for `t>=1`. For `t`
between specified points, linear interpolation is used to determine
`sfLinear$spend`.

The cumulative spending returned in `sfStep$spend` is 0 for `t<param[1]`
and `alpha` for `t>=1`. Letting `m <- length(param/2)`, for
`i=1,2,...m-1` and ` param[i]<= t < param[i+1]`, the cumulative spending
is set at `alpha * param[i+m]` (also for `param[m]<=t<1`).

Note that if `param[2m]` is 1, then the first time an analysis is
performed after the last proportion of final planned information
(`param[m]`) will be the final analysis, using any remaining error that
was not previously spent.

See
[`vignette("SpendingFunctionOverview")`](https://leejtding.github.io/gsDesignCRT/articles/SpendingFunctionOverview.md)
for further details.

## Note

The gsDesign technical manual is available at
<https://keaven.github.io/gsd-tech-manual/>.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## See also

[`vignette("SpendingFunctionOverview")`](https://leejtding.github.io/gsDesignCRT/articles/SpendingFunctionOverview.md),
[`gsDesignCRT`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md),
`vignette("gsDesignCRTPackageOverview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
