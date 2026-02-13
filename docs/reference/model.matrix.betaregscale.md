# Extract design matrix

Extract design matrix

## Usage

``` r
# S3 method for class 'betaregscale'
model.matrix(object, model = c("mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: `"mean"` (default) or `"precision"`.

- ...:

  Ignored.

## Value

The design matrix for the specified submodel.
