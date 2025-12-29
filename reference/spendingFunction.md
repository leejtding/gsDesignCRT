# Spending Function

Spending Function

## Usage

``` r
spendingFunction(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Defaults in calls to
  [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
  are `alpha=0.025` for one-sided Type I error specification and
  `alpha=0.1` for Type II error specification. However, this could be
  set to 1 if, for descriptive purposes, you wish to see the proportion
  of spending as a function of the proportion of sample
  size/information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size/information for which the
  spending function will be computed.

- param:

  A single real value or a vector of real values specifying the spending
  function parameter(s); this must be appropriately matched to the
  spending function specified.

## Value

`spendingFunction` and spending functions in general produce an object
of type `spendfn`.

- name:

  A character string with the name of the spending function.

- param:

  any parameters used for the spending function.

- parname:

  a character string or strings with the name(s) of the parameter(s) in
  `param`.

- sf:

  the spending function specified.

- spend:

  a vector of cumulative spending values corresponding to the input
  values in `t`.

- bound:

  this is null when returned from the spending function, but is set in
  [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
  if the spending function is called from there. Contains z-values for
  bounds of a design.

- prob:

  this is null when returned from the spending function, but is set in
  [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
  if the spending function is called from there. Contains probabilities
  of boundary crossing at `i`-th analysis for `j`-th theta value input
  to
  [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
  in `prob[i,j]`.

## Note

The gsDesign technical manual is available at
<https://keaven.github.io/gsd-tech-manual/>.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## See also

[`gsDesignCRT`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md),
[`sfHSD`](https://leejtding.github.io/gsDesignCRT/reference/sfHSD.md),
[`sfPower`](https://leejtding.github.io/gsDesignCRT/reference/sfPower.md),
[`sfLogistic`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md),
[`sfExponential`](https://leejtding.github.io/gsDesignCRT/reference/sfExponential.md),
[`sfTruncated`](https://leejtding.github.io/gsDesignCRT/reference/sfSpecial.md),
`vignette("gsDesignCRTPackageOverview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
