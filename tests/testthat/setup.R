# Suppress plot output during tests to avoid Rplots.pdf creation
if (!interactive()) {
  pdf(nullfile())
}
