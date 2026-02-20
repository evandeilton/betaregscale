# Print a fitted model (brief betareg style)

Print a fitted model (brief betareg style)

## Usage

``` r
# S3 method for class 'brs'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A fitted `"betaregscale"` object.

- digits:

  Number of significant digits.

- ...:

  Included for consistency with generic methods. Currently passed to
  internal methods where applicable.

## Value

Invisibly returns the input object `x`. The function is called for its
side effect of printing a formatted summary of the fitted model to the
console, including the model call, mean coefficients (with link
function), and precision coefficients (with link function).
