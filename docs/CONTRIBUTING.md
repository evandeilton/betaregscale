# Contributing to betaregscale

Thank you for your interest in contributing to betaregscale! This
document provides guidelines for contributing to the project.

## How to contribute

### Reporting issues

If you find a problem, please open an issue on
[GitHub](https://github.com/evandeilton/betaregscale/issues) with:

- A minimal reproducible example
- The output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
- A clear description of the expected vs. actual behavior

### Suggesting features

Feature requests are welcome. Please open an issue describing the use
case and expected behavior.

### Pull requests

1.  Fork the repository and create a feature branch from `main`.
2.  Follow the existing code style (snake_case, roxygen2 documentation).
3.  Add tests for new functionality using `testthat`.
4.  Run
    [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
    and ensure 0 errors, 0 warnings, 0 notes.
5.  Update documentation with
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    if needed.
6.  Submit a pull request with a clear description of the changes.

## Development setup

``` r

# Install development dependencies
install.packages(c("devtools", "testthat", "roxygen2", "knitr", "rmarkdown"))

# Clone and install
# git clone https://github.com/evandeilton/betaregscale.git
devtools::install_dev_deps()
devtools::load_all()

# Run tests
devtools::test()

# Full check
devtools::check()
```

## Code style

- Use `snake_case` for function and variable names.
- Document all exported functions with roxygen2.
- Internal functions should use `@keywords internal`.
- Keep lines under 80 characters when possible.

## Package structure

- `R/` – R source files
- `src/` – C++ source (Rcpp/RcppArmadillo)
- `tests/testthat/` – unit tests
- `vignettes/` – package vignettes
- `man/` – auto-generated documentation (do not edit manually)

## License

By contributing, you agree that your contributions will be licensed
under the MIT License.
