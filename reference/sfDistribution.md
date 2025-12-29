# Two-parameter Spending Function Families

The functions `sfLogistic()`, `sfNormal()`, `sfExtremeValue()`,
`sfExtremeValue2()`, `sfCauchy()`, and `sfBetaDist()` are all
2-parameter spending function families. These provide increased
flexibility in some situations where the flexibility of a one-parameter
spending function family is not sufficient. These functions all allow
fitting of two points on a cumulative spending function curve; in this
case, four parameters are specified indicating an x and a y coordinate
for each of 2 points.

`sfBetaDist(alpha,t,param)` is simply `alpha` times the incomplete beta
cumulative distribution function with parameters \\a\\ and \\b\\ passed
in `param` evaluated at values passed in `t`.

The other spending functions take the form \$\$f(t;\alpha,a,b)=\alpha
F(a+bF^{-1}(t))\$\$ where \\F()\\ is a cumulative distribution function
with values \\\> 0\\ on the real line (logistic for `sfLogistic()`,
normal for `sfNormal()`, extreme value for `sfExtremeValue()` and Cauchy
for `sfCauchy()`) and \\F^{-1}()\\ is its inverse.

For the logistic spending function this simplifies to
\$\$f(t;\alpha,a,b)=\alpha (1-(1+e^a(t/(1-t))^b)^{-1}).\$\$

For the extreme value distribution with \$\$F(x)=\exp(-\exp(-x))\$\$
this simplifies to \$\$f(t;\alpha,a,b)=\alpha \exp(-e^a (-\ln t)^b).\$\$
Since the extreme value distribution is not symmetric, there is also a
version where the standard distribution is flipped about 0. This is
reflected in `sfExtremeValue2()` where \$\$F(x)=1-\exp(-\exp(x)).\$\$

## Usage

``` r
sfLogistic(alpha, t, param)

sfBetaDist(alpha, t, param)

sfCauchy(alpha, t, param)

sfExtremeValue(alpha, t, param)

sfExtremeValue2(alpha, t, param)

sfNormal(alpha, t, param)
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

  In the two-parameter specification, `sfBetaDist()` requires 2 positive
  values, while `sfLogistic()`, `sfNormal()`, `sfExtremeValue()`,

  `sfExtremeValue2()` and `sfCauchy()` require the first parameter to be
  any real value and the second to be a positive value. The four
  parameter specification is `c(t1,t2,u1,u2)` where the objective is
  that `sf(t1)=alpha*u1` and `sf(t2)=alpha*u2`. In this
  parameterization, all four values must be between 0 and 1 and
  `t1 < t2`, `u1 < u2`.

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
