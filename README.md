
[![pkgdown](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# betaroti

<!-- badges: start -->
<!-- badges: end -->

Um pacote de funções para realizar ajuste de um mdoelo de regressão beta
para dados com resposta ordinal transformada intervalar

## Instalação

Você pode instalar o pacote com esse comando abaixo.

``` r
# install.packages("devtools")
devtools::install_github("evandeilton/betaroti")
```

## Exemplos

Esses são alguns exemplos de uso das funções do pacote.

``` r
library(betaroti)
```

### Simula dados do modelo beta ordinal com dispersão fixa

- Esta função simula dados de uma variável beta ordinal com dispersão
  fixa, aplicando diferentes funções de ligação.

``` r
# Criar um conjunto de dados de exemplo
set.seed(42)
dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

# Simular dados usando a função beta_ordinal_simula_dados com parâmetros personalizados
dados_simulados <- beta_ordinal_simula_dados(
  formula = ~ x1 + x2,
  dados = dados,
  betas = c(1,-0.3, 0.4),
  phi = 30,
  link = "probit",
  ncuts = 100,
  type = "m"
)
dados_simulados %>%
  head() %>%
  knitr::kable(digits = 4, caption = "Dados simulados")
```

|      y |  yr |  left | right |      x1 |      x2 |
|-------:|----:|------:|------:|--------:|--------:|
| 0.9623 |  96 | 0.955 | 0.965 |  1.3710 |  1.2010 |
| 0.9221 |  92 | 0.915 | 0.925 | -0.5647 |  1.0448 |
| 0.5581 |  56 | 0.555 | 0.565 |  0.3631 | -1.0032 |
| 0.9856 |  99 | 0.985 | 0.995 |  0.6329 |  1.8485 |
| 0.8270 |  83 | 0.825 | 0.835 |  0.4043 | -0.6668 |
| 0.9044 |  90 | 0.895 | 0.905 | -0.1061 |  0.1055 |

Dados simulados

### Simula dados do modelo beta ordinal com dispersão variável

Esta função simula dados de uma variável beta ordinal, aplicando
diferentes funções de ligação tanto para betas como para zetas de phi

``` r
# Criar um conjunto de dados de exemplo
set.seed(421)
n <- 50
fx <- ~ x1 + x2
fz <- ~ z1

dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = runif(n),
  z2 = runif(n)
)

dados_simulados <- beta_ordinal_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2,-0.5, 0.3),
  zetas = c(1, 1.2),
  link_x = "logit",
  link_z = "log",
  ncuts = 100,
  type = "m"
)

dados_simulados %>% 
  head() %>%
  knitr::kable(digits = 4, caption = "Dados simulados")
```

|      y |  yr |  left | right |      x1 |      x2 |     z1 |
|-------:|----:|------:|------:|--------:|--------:|-------:|
| 0.3492 |  35 | 0.345 | 0.355 |  0.8048 |  1.3698 | 0.4632 |
| 0.6844 |  68 | 0.675 | 0.685 |  0.5694 |  0.1462 | 0.3954 |
| 0.1605 |  16 | 0.155 | 0.165 |  1.0160 | -0.3881 | 0.3554 |
| 0.3401 |  34 | 0.335 | 0.345 |  1.2838 |  0.4173 | 0.1075 |
| 0.4315 |  43 | 0.425 | 0.435 |  0.0747 |  0.8039 | 0.5592 |
| 0.7773 |  78 | 0.775 | 0.785 | -0.4731 | -0.1909 | 0.5579 |

Dados simulados

### Log-verossimilhança do modelo beta ordinal com dispersão fixa

Esta função calcula a log-verossimilhança de um conjunto de dados para
uma variável beta ordinal, aplicando diferentes funções de ligação.

``` r
# Criar um conjunto de dados de exemplo
set.seed(421)
dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
# Calcular a log-verossimilhança usando a função log_vero_beta_ordinal
param <- c(0, 0.5,-0.2, 50)
phi <- 30
formula <- ~ x1 + x2
dados_simulados <-
  beta_ordinal_simula_dados(
    formula = formula,
    dados = dados,
    betas = c(0, 0.5,-0.2),
    phi = phi,
    link = "logit",
    ncuts = 10,
    type = "m"
  )
log_verossimilhanca <-
  beta_ordinal_log_vero(param, formula, dados_simulados)
print(log_verossimilhanca)
#> [1] -156.2426
```

### Log-verossimilhança do modelo beta ordinal com dispersão variável

Esta função calcula a log-verossimilhança de um conjunto de dados para
uma variável beta ordinal, aplicando diferentes funções de ligação tanto
no betas de mu como no zetas de phi.

``` r
# Criar um conjunto de dados de exemplo
n <- 50
fx <- ~ x1 + x2
fz <- ~ z1

dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = runif(n),
  z2 = runif(n)
)

dados_simulados <- beta_ordinal_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2,-0.5, 0.3),
  zetas = c(1, 1.2),
  link_x = "logit",
  link_z = "log",
  ncuts = 100
)
# Calcular a log-verossimilhança usando a função log_vero_beta_ordinal
log_verossimilhanca <- beta_ordinal_log_vero_z(
  param = c(c(0.2,-0.5, 0.3), c(1, 1.2)),
  formula_x = fx,
  formula_z = fz,
  dados = dados_simulados,
  link_x = "logit",
  link_z = "log",
  acumulada = TRUE
)
print(log_verossimilhanca)
#> [1] -212.9946
```

## Função para ajustar um modelo beta ordinal com dispersão fixa

A função beta_ordinal_fit ajusta um modelo beta ordinal usando a função
optim do pacote stats, retornando uma lista com as principais
informações do ajuste.

``` r
# Criar um conjunto de dados de exemplo
set.seed(42)
n <- 100
betas <- c(0.2, 0.3,-0.4, 0.1)
formula <- ~ x1 + x2 + x3
phi <- 50

dados <- data.frame(
  x1 = rnorm(n, mean = 1, sd = 0.5),
  x2 = rbinom(n, size = 1, prob = 0.5),
  x3 = rnorm(n, mean = 2, sd = 1)
)

dados_simulados <- beta_ordinal_simula_dados(
  formula = formula,
  dados = dados,
  betas = betas,
  phi = phi,
  link = "logit",
  ncuts = 100
)

fit <- beta_ordinal_fit(
  formula = formula,
  dados = dados_simulados,
  link = "probit",
  num_hessiana = TRUE
)
coe <- beta_ordinal_coef(fit)

coe$est %>% 
  knitr::kable(digits = 4, caption = "Estimativas")
```

| variable    | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept) |   0.2156 |   0.0353 |   0.3958 | 0.0920 |  2.3437 |  0.0212 |
| x1          |   0.3067 |   0.1965 |   0.4168 | 0.0562 |  5.4551 |  0.0000 |
| x2          |  -0.4558 |  -0.5734 |  -0.3381 | 0.0600 | -7.5940 |  0.0000 |
| x3          |   0.1003 |   0.0376 |   0.1630 | 0.0320 |  3.1357 |  0.0023 |
| (phi)       |  51.1296 |  36.7673 |  65.4919 | 7.3278 |  6.9774 |  0.0000 |

Estimativas

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "Bondade do ajuste")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|       330.5828 | -649.1657 | -633.5347 |

Bondade do ajuste

## Função para ajustar um modelo beta ordinal com dispersão variável

A função beta_ordinal_fit_z ajusta um modelo beta ordinal com dispersão
variável, isto é, com um preditor para os betas e outro para o phi, onde
covariáveis são aplicadas para explicar a dispersão. A função usa optim
do pacote stats e retorna uma lista com as principais informações do
ajuste.

``` r
n <- 50
fx <- ~ x1 + x2
fz <- ~ z1

dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = runif(n),
  z2 = runif(n)
)

dados_simulados <- beta_ordinal_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2,-0.5, 0.3),
  zetas = c(1, 1.2),
  link_x = "logit",
  link_z = "log",
  ncuts = 100
)
fit_z <- beta_ordinal_fit_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados_simulados,
  link_x = "logit",
  link_z = "log",
  num_hessiana = TRUE
)
coe <- beta_ordinal_coef(fit_z)

coe$est %>% 
  knitr::kable(digits = 4, caption = "Estimativas")
```

| variable           | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:-------------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept)        |   0.1225 |  -0.1384 |   0.3833 | 0.1331 |  0.9203 |  0.3623 |
| x1                 |  -0.3995 |  -0.6869 |  -0.1121 | 0.1466 | -2.7247 |  0.0091 |
| x2                 |   0.3626 |   0.0953 |   0.6299 | 0.1364 |  2.6591 |  0.0108 |
| (phi)\_(Intercept) |   0.5855 |  -0.2018 |   1.3728 | 0.4017 |  1.4576 |  0.1519 |
| (phi)\_z1          |   1.3097 |  -0.0996 |   2.7189 | 0.7190 |  1.8214 |  0.0752 |

Estimativas

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "Bondade do ajuste")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|       222.5784 | -433.1567 | -421.6846 |

Bondade do ajuste

### Coleta estatística do ajuste

Coleta as estimativas e suas estatísticas para um objeto de ajuste do
modelo beta ordinal com censura intervalar. Coleta também as
estatísticas de bondade do ajuste como a log-verossimilhança, o AIC e o
BIC do modelo tanto para o modelo de dispersão fixa como aquele com
dispersão variável.

``` r
n <- 50
fx <- ~ x1 + x2
fz <- ~ z1

dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = runif(n),
  z2 = runif(n)
)

dados_simulados <- beta_ordinal_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2,-0.5, 0.3),
  zetas = c(1, 1.2),
  link_x = "logit",
  link_z = "log",
  ncuts = 100
)
fit_z <- beta_ordinal_fit_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados_simulados,
  link_x = "logit",
  link_z = "log",
  num_hessiana = TRUE
)
coe <- beta_ordinal_coef(fit_z)

coe$est %>% 
  knitr::kable(digits = 4, caption = "Estimativas")
```

| variable           | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:-------------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept)        |   0.2752 |   0.0765 |   0.4738 | 0.1013 |  2.7154 |  0.0094 |
| x1                 |  -0.5384 |  -0.7388 |  -0.3380 | 0.1023 | -5.2649 |  0.0000 |
| x2                 |   0.4141 |   0.2219 |   0.6064 | 0.0981 |  4.2227 |  0.0001 |
| (phi)\_(Intercept) |   1.6446 |   0.8457 |   2.4436 | 0.4076 |  4.0346 |  0.0002 |
| (phi)\_z1          |   0.7760 |  -0.5975 |   2.1495 | 0.7008 |  1.1074 |  0.2740 |

Estimativas

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "Bondade do ajuste")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|        205.461 | -398.9221 | -387.4499 |

Bondade do ajuste

## Updating!
