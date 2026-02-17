test_that("anova() compares brs and brsmm models in evolutionary workflow", {
  skip_on_cran()

  set.seed(411)
  g <- 10
  ni <- 10
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)
  b0 <- rnorm(g, sd = 0.45)
  b1 <- rnorm(g, sd = 0.20)

  eta <- 0.15 + 0.55 * x1 + b0[id] + b1[id] * x1
  mu <- plogis(eta)
  phi <- plogis(-0.25)
  shp <- brs_repar(mu = mu, phi = phi, repar = 2)
  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
  dat <- data.frame(y = y, x1 = x1, id = id)

  m0 <- brs(y ~ x1, data = dat, repar = 2)
  m1 <- brsmm(
    y ~ x1,
    random = ~ 1 | id,
    data = dat,
    repar = 2,
    int_method = "laplace",
    control = list(maxit = 700)
  )
  m2 <- brsmm(
    y ~ x1,
    random = ~ 1 + x1 | id,
    data = dat,
    repar = 2,
    int_method = "laplace",
    control = list(maxit = 900)
  )

  tab <- anova(m0, m1, m2, test = "Chisq")
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 3L)
  expect_true(all(c("Df", "logLik", "AIC", "BIC", "Chisq", "Chi Df", "Pr(>Chisq)") %in% colnames(tab)))
  expect_true(all(diff(tab$Df) >= 0))
  expect_true(is.finite(tab$Chisq[2]))
  expect_true(is.finite(tab$`Pr(>Chisq)`[2]))
  expect_true(is.finite(tab$Chisq[3]))
  expect_true(is.finite(tab$`Pr(>Chisq)`[3]))
})

test_that("anova() validates compatible observation counts", {
  skip_on_cran()

  dat <- data.frame(y = round(runif(80, 0, 100)), x1 = rnorm(80))
  m_a <- brs(y ~ x1, data = dat, repar = 2)
  m_b <- brs(y ~ x1, data = dat[1:70, , drop = FALSE], repar = 2)

  expect_error(
    anova(m_a, m_b),
    "same number of observations"
  )
})
