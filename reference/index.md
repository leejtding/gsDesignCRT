# Package index

## Group Sequential CRT Sample Size Computation

For an overview of how to use the gsDesignCRT package, see
[`vignette("example")`](https://leejtding.github.io/gsDesignCRT/articles/example.md).

- [`gsDesignCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsDesignCRT.md)
  : Compute stopping boundaries, maximum sample size, and expected
  sample sizes for a group sequential cluster randomized trial.
- [`print(`*`<gsDesignCRT>`*`)`](https://leejtding.github.io/gsDesignCRT/reference/print.gsDesignCRT.md)
  : Print group sequential CRT design
- [`plot(`*`<gsDesignCRT>`*`)`](https://leejtding.github.io/gsDesignCRT/reference/plot.gsDesignCRT.md)
  : Plot group sequential CRT design

## Group Sequential CRT Boundary Computation

- [`gsUpperCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsUpperCRT.md)
  : Boundary derivation for efficacy stopping only.
- [`gsLowerCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsLowerCRT.md)
  : Boundary derivation for binding or non-binding futility stopping
  only.
- [`gsBoundsCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsBoundsCRT.md)
  : Boundary derivation for efficacy and binding or non-binding futility
  stopping.

## Group Sequential CRT Stopping Probability

- [`gsProbabilityCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsProbabilityCRT.md)
  : Compute stopping boundary crossing probabilities.

## Simulate Group Sequential CRT

- [`gsSimCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsSimCRT.md)
  : Simulate group sequential cluster-randomized trial.
- [`gsSimContCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsSimContCRT.md)
  : Simulate group sequential cluster-randomized trial with continuous
  outcomes.
- [`gsSimBinCRT()`](https://leejtding.github.io/gsDesignCRT/reference/gsSimBinCRT.md)
  : Simulate group sequential cluster-randomized trial with binary
  outcomes
- [`genContCRT()`](https://leejtding.github.io/gsDesignCRT/reference/genContCRT.md)
  : Simulate cluster-randomized trial data with continuous outcomes
- [`genBinCRT()`](https://leejtding.github.io/gsDesignCRT/reference/genBinCRT.md)
  : Simulate cluster-randomized trial data with binary outcomes

## Spending Functions

For an overview of spending functions, see
[`vignette("SpendingFunctionOverview")`](https://leejtding.github.io/gsDesignCRT/articles/SpendingFunctionOverview.md).

- [`spendingFunction()`](https://leejtding.github.io/gsDesignCRT/reference/spendingFunction.md)
  : Spending Function
- [`sfLDOF()`](https://leejtding.github.io/gsDesignCRT/reference/sfLDOF.md)
  [`sfLDPocock()`](https://leejtding.github.io/gsDesignCRT/reference/sfLDOF.md)
  : Lan-DeMets Spending function overview
- [`sfHSD()`](https://leejtding.github.io/gsDesignCRT/reference/sfHSD.md)
  : Hwang-Shih-DeCani Spending Function
- [`sfPower()`](https://leejtding.github.io/gsDesignCRT/reference/sfPower.md)
  : Kim-DeMets (power) Spending Function
- [`sfExponential()`](https://leejtding.github.io/gsDesignCRT/reference/sfExponential.md)
  : Exponential Spending Function
- [`sfLogistic()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  [`sfBetaDist()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  [`sfCauchy()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  [`sfExtremeValue()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  [`sfExtremeValue2()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  [`sfNormal()`](https://leejtding.github.io/gsDesignCRT/reference/sfDistribution.md)
  : Two-parameter Spending Function Families
- [`sfTDist()`](https://leejtding.github.io/gsDesignCRT/reference/sfTDist.md)
  : t-distribution Spending Function
- [`sfLinear()`](https://leejtding.github.io/gsDesignCRT/reference/sfLinear.md)
  [`sfStep()`](https://leejtding.github.io/gsDesignCRT/reference/sfLinear.md)
  : Piecewise Linear and Step Function Spending Functions
- [`sfPoints()`](https://leejtding.github.io/gsDesignCRT/reference/sfPoints.md)
  : Pointwise Spending Function
- [`sfTruncated()`](https://leejtding.github.io/gsDesignCRT/reference/sfSpecial.md)
  [`sfTrimmed()`](https://leejtding.github.io/gsDesignCRT/reference/sfSpecial.md)
  [`sfGapped()`](https://leejtding.github.io/gsDesignCRT/reference/sfSpecial.md)
  : Truncated, trimmed and gapped spending functions

## Utility Functions

- [`checkLengths()`](https://leejtding.github.io/gsDesignCRT/reference/checkScalar.md)
  [`checkRange()`](https://leejtding.github.io/gsDesignCRT/reference/checkScalar.md)
  [`checkScalar()`](https://leejtding.github.io/gsDesignCRT/reference/checkScalar.md)
  [`checkVector()`](https://leejtding.github.io/gsDesignCRT/reference/checkScalar.md)
  [`isInteger()`](https://leejtding.github.io/gsDesignCRT/reference/checkScalar.md)
  : Utility functions to verify variable properties
