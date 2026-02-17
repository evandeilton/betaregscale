# Model comparison by analysis of deviance (LR test) for \`brs\`

Model comparison by analysis of deviance (LR test) for \`brs\`

## Usage

``` r
# S3 method for class 'brs'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  A fitted \`brs\` model.

- ...:

  Additional fitted \`brs\` and/or \`brsmm\` models.

- test:

  Character; \`"Chisq"\` (default) or \`"none"\`.

## Value

An object of class \`anova\` and \`data.frame\` with model-wise
log-likelihood, information criteria, and (optionally) LR test columns.
