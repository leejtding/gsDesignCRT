# Lan-DeMets Spending function overview

Lan and DeMets (1983) first published the method of using spending
functions to set boundaries for group sequential trials. In this
publication they proposed two specific spending functions: one to
approximate an O'Brien-Fleming design and the other to approximate a
Pocock design. The spending function to approximate O'Brien-Fleming has
been generalized as proposed by Liu, et al (2012)

With `param=1=rho`, the Lan-DeMets (1983) spending function to
approximate an O'Brien-Fleming bound is implemented in the function
(`sfLDOF()`): \$\$f(t; \alpha)=2-2\Phi\left(\Phi^{-1}(1-\alpha/2)/
t^{\rho/2}\right).\$\$ For `rho` otherwise in `[.005,2]`, this is the
generalized version of Liu et al (2012). For `param` outside of
`[.005,2]`, `rho` is set to 1. The Lan-DeMets (1983) spending function
to approximate a Pocock design is implemented in the function
`sfLDPocock()`: \$\$f(t;\alpha)=\alpha ln(1+(e-1)t).\$\$ As shown in
examples below, other spending functions can be used to ge t as good or
better approximations to Pocock and O'Brien-Fleming bounds. In
particular, O'Brien-Fleming bounds can be closely approximated using
[`sfExponential`](https://leejtding.github.io/gsDesignCRT/reference/sfExponential.md).

## Usage

``` r
sfLDOF(alpha, t, param = NULL)

sfLDPocock(alpha, t, param)
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

  This parameter is not used for `sfLDPocock`, not required for `sfLDOF`
  and need not be specified. For `sfLDPocock` it is here so that the
  calling sequence conforms to the standard for spending functions used
  with
  [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md).
  For `sfLDOF` it will default to 1 (Lan-DeMets function to approximate
  O'Brien-Fleming) if `NULL` or if outside of the range `[.005,2]`.
  otherwise, it will be use to set rho from Liu et al (2012).

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

Lan, KKG and DeMets, DL (1983), Discrete sequential boundaries for
clinical trials. *Biometrika*;70: 659-663.

Liu, Q, Lim, P, Nuamah, I, and Li, Y (2012), On adaptive error spending
approach for group sequential trials with random information levels.
*Journal of biopharmaceutical statistics*; 22(4), 687-699.

## See also

[`vignette("SpendingFunctionOverview")`](https://leejtding.github.io/gsDesignCRT/articles/SpendingFunctionOverview.md),
[`gsDesignCRT`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md),
`vignette("gsDesignPackageOverview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
