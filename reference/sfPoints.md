# Pointwise Spending Function

The function `sfPoints` implements a spending function with values
specified for an arbitrary set of specified points. It is now
recommended to use sfLinear rather than sfPoints. If using `sfPoints()`
in a design, it is recommended to specify how to interpolate between the
specified points (e.g,, linear interpolation); also consider fitting
smooth spending functions; see
[`vignette("SpendingFunctionOverview")`](https://leejtding.github.io/gsDesignCRT/articles/SpendingFunctionOverview.md).

## Usage

``` r
sfPoints(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size/information.

- t:

  A vector of points with increasing values from \>0 and \<=1. Values of
  the proportion of sample size/information for which the spending
  function will be computed.

- param:

  A vector of the same length as `t` specifying the cumulative
  proportion of spending to corresponding to each point in `t`; must be
  \>=0 and \<=1.

## Value

An object of type `spendfn`. See
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
`vignette("gsDesignCRTPackageOverview")`,
[sfLogistic](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)

## Author

Keaven Anderson <keaven_anderson@merck.com>
