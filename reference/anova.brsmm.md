# Model comparison by analysis of deviance (LR test) for \`brsmm\`

Model comparison by analysis of deviance (LR test) for \`brsmm\`

## Usage

``` r
# S3 method for class 'brsmm'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  A fitted \`brsmm\` model.

- ...:

  Additional fitted \`brsmm\` and/or \`brs\` models.

- test:

  Character; \`"Chisq"\` (default) or \`"none"\`.

## Value

An object of class \`anova\` and \`data.frame\` with model-wise
log-likelihood, information criteria, and (optionally) LR test columns.
