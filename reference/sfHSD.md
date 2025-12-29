# Hwang-Shih-DeCani Spending Function

The function `sfHSD` implements a Hwang-Shih-DeCani spending function.

A Hwang-Shih-DeCani spending function takes the form \$\$f(t;\alpha,
\gamma)=\alpha(1-e^{-\gamma t})/(1-e^{-\gamma})\$\$ where \\\gamma\\ is
the value passed in `param`. A value of \\\gamma=-4\\ is used to
approximate an O'Brien-Fleming design (see
[`sfExponential`](https://leejtding.github.io/gsDesignCRT/reference/sfExponential.md)
for a better fit), while a value of \\\gamma=1\\ approximates a Pocock
design well.

## Usage

``` r
sfHSD(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size/information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size/information for which the
  spending function will be computed.

- param:

  A single real value specifying the gamma parameter for which
  Hwang-Shih-DeCani spending is to be computed; allowable range is
  \[-40, 40\]

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
`vignette("gsDesignCRTPackageOverview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
