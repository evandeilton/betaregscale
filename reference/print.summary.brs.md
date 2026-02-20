# Print a model summary (betareg style)

Print a model summary (betareg style)

## Usage

``` r
# S3 method for class 'summary.brs'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A `"summary.betaregscale"` object.

- digits:

  Number of digits.

- ...:

  Passed to `printCoefmat`.

## Value

Invisibly returns the input object `x`. The function is called for its
side effect of printing a comprehensive summary to the console,
including the model call, quantile residuals, coefficient tables for
mean and precision submodels with significance stars, goodness-of-fit
statistics (log-likelihood, pseudo R-squared), optimization details, and
censoring information.
