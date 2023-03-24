
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betaroti

<!-- badges: start -->

[![pkgdown](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

Este pacote foi desenvolvido para fornecer uma biblioteca de funções em
R, especialmente projetadas para ajustar modelos de regressão beta em
dados com resposta ordinal transformada intervalar, tanto com dispersão
fixa quanto variável. Além disso, permite realizar simulações para
avaliar o desempenho dos modelos no processo de estimação. Acesse o
[repositório oficial](https://github.com/evandeilton/betaroti) no GitHub
para visualizar o código-fonte e contribuir com o projeto. Nesta página,
você encontrará informações detalhadas sobre a instalação, uso e
funcionamento interno das funções oferecidas pelo pacote “betaroti”.

O pacote “betaroti” foi criado para facilitar a modelagem e análise de
dados em situações onde a variável resposta é ordinal numérica, mas pode
ser transformada em um intervalo contínuo, que forma um tipo de dado
limitado, ex. $y = (y_s;y_i)$. Nesses casos, pode haver censura à
esquerda, censura à direita ou intervalar.

Aplicações desse tipo de modelo encontram espaço em dados de pesquisas
de opinião, avaliações de produtos, escala de dor, avaliação de reação
de compostos quimicos, etc. Ao utilizar a distribuição beta, o pacote
permite acomodar características específicas dos dados numa estrutura
que permite associar variáveis explicativas com a variável resposta
intervalar, através de uma estrutura de regressão covenciente e que
permite o uso de preditores lineares tanto pros coeficientes
relacionados com a média como para aqueles relacionados com a dispersão,
fornecendo estimativas robustas e confiáveis dos parâmetros do modelo.

## Principais funcionalidades

O pacote “betaroti” oferece uma série de funções úteis para lidar com
modelos de regressão beta e dados com resposta ordinal transformada
intervalar, abrangendo cenários com dispersão fixa e variável. As
principais funcionalidades incluem:

- Ajuste de modelos de regressão beta com dispersão fixa e variável.
- Funções para simulação de dados, permitindo a avaliação do desempenho
  dos modelos em diferentes cenários.
- Estatística de bondade do ajuste como AIC e BIC.
- Funções para ajuste e comparação de modelos com diferentes combinações
  de variáveis explicativas tanto para $\mu$ como $\phi$.

Acesse a documentação detalhada de cada função e exemplos de uso neste
site para obter informações sobre como utilizar o pacote “betaroti” em
suas análises.

## Instalação

Você pode instalar o pacote com esse comando abaixo.

``` r
if(!require(betaroti)){
  devtools::install_github("evandeilton/betaroti")  
}
require(betaroti)
```

## Exemplos

Esses são alguns exemplos de uso das funções do pacote.

### Simula dados do modelo beta ordinal com dispersão fixa

Esta função gera amostras de uma variável beta ordinal com dispersão
fixa, empregando diversas funções de ligação.

Neste bloco de código R, apresentamos um exemplo de como utilizar a
função beta_ordinal_simula_dados para simular dados de uma variável beta
ordinal com dispersão fixa. Segue uma descrição detalhada do processo:

- Criamos um conjunto de dados de exemplo com 100 observações e duas
  variáveis explicativas independentes (x1 e x2), geradas a partir de
  uma distribuição normal:

``` r

# Criar um conjunto de dados de exemplo
set.seed(42)
dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
```

- Em seguida, utilizamos a função `beta_ordinal_simula_dados` para
  simular dados com base nos parâmetros personalizados fornecidos. A
  função recebe os seguintes argumentos:
- `formula`: especifica a relação entre a variável resposta e as
  variáveis explicativas.
- `dados`: fornece o conjunto de dados de entrada.
- `betas`: vetor de coeficientes de regressão.
- `phi`: valor do parâmetro de dispersão.
- `link`: função de ligação a ser utilizada (neste caso, “probit”).
- `ncuts`: número de pontos de corte para a discretização da variável
  resposta.
- `type`: tipo de tratamento do intervalo. `m` centraliza `y` ao meio.
  Ex. Se foi coletado o valor $y = 6$, transforma-se $y_t = 6/10 = 0.6$.
  Assim, para tratar a incerteza do instrumento, sugere-se que a medida
  anotada pode estar limitada a $y_{left} = 5.5$ e $y_{right} = 6.6$.

``` r
dados_simulados <- beta_ordinal_simula_dados(
  formula = ~ x1 + x2,
  dados = dados,
  betas = c(1, -0.3, 0.4),
  phi = 30,
  link = "probit",
  ncuts = 100,
  type = "m"
)
dados_simulados %>%
  head() %>%
  knitr::kable(digits = 4, caption = "")
```

|      y |  yr |  left | right |      x1 |      x2 |
|-------:|----:|------:|------:|--------:|--------:|
| 0.9623 |  96 | 0.955 | 0.965 |  1.3710 |  1.2010 |
| 0.9221 |  92 | 0.915 | 0.925 | -0.5647 |  1.0448 |
| 0.5581 |  56 | 0.555 | 0.565 |  0.3631 | -1.0032 |
| 0.9856 |  99 | 0.985 | 0.995 |  0.6329 |  1.8485 |
| 0.8270 |  83 | 0.825 | 0.835 |  0.4043 | -0.6668 |
| 0.9044 |  90 | 0.895 | 0.905 | -0.1061 |  0.1055 |

### Simula dados provenientes de um modelo beta ordinal com dispersão variável.

Neste bloco de código R, é criado um conjunto de dados simulados de um
modelo beta ordinal com dispersão variável utilizando a função
`beta_ordinal_simula_dados_z.` O processo é resumido abaixo:

- Definir semente e tamanho da amostra, além das fórmulas para as
  variáveis explicativas x e z.

- Criar um conjunto de dados de exemplo com 50 observações e quatro
  variáveis independentes (x1, x2, z1 e z2), geradas a partir de
  distribuições normal e uniforme.

- Utilizar a função `beta_ordinal_simula_dados_z` para gerar dados
  simulados com base nos parâmetros fornecidos, como fórmulas,
  coeficientes de regressão, funções de ligação e número de pontos de
  corte.

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
  knitr::kable(digits = 4, caption = "")
```

|      y |  yr |  left | right |      x1 |      x2 |     z1 |
|-------:|----:|------:|------:|--------:|--------:|-------:|
| 0.3492 |  35 | 0.345 | 0.355 |  0.8048 |  1.3698 | 0.4632 |
| 0.6844 |  68 | 0.675 | 0.685 |  0.5694 |  0.1462 | 0.3954 |
| 0.1605 |  16 | 0.155 | 0.165 |  1.0160 | -0.3881 | 0.3554 |
| 0.3401 |  34 | 0.335 | 0.345 |  1.2838 |  0.4173 | 0.1075 |
| 0.4315 |  43 | 0.425 | 0.435 |  0.0747 |  0.8039 | 0.5592 |
| 0.7773 |  78 | 0.775 | 0.785 | -0.4731 | -0.1909 | 0.5579 |

### Log-verossimilhança do modelo beta ordinal com dispersão fixa

Esta função calcula a log-verossimilhança de um conjunto de dados para
uma variável beta ordinal, aplicando diferentes funções de ligação.

Neste bloco de código R, calculamos a log-verossimilhança de um modelo
beta ordinal com dispersão constante utilizando a função
`beta_ordinal_log_vero`. O processo é resumido abaixo:

- Definir semente e criar um conjunto de dados de exemplo com 100
  observações e duas variáveis explicativas independentes (x1 e x2),
  geradas a partir de uma distribuição normal.

- Definir os parâmetros do modelo, incluindo coeficientes de regressão,
  parâmetro de dispersão e a fórmula da relação entre as variáveis.

- Utilizar a função beta_ordinal_simula_dados para gerar dados simulados
  com base nos parâmetros fornecidos.

- Calcular a log-verossimilhança do modelo ajustado com os dados
  simulados usando a função `beta_ordinal_log_vero`.

Como resultado, obtemos a log-verossimilhança do modelo ajustado aos
dados simulados, que é uma medida de quão bem o modelo se ajusta aos
dados observados.

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

Neste bloco de código R, calculamos a log-verossimilhança de um modelo
beta ordinal com dispersão variável utilizando a função
`beta_ordinal_log_vero_z`. O processo é resumido abaixo:

- Definir o tamanho da amostra e as fórmulas para as variáveis
  explicativas x e z.

- Criar um conjunto de dados de exemplo com 50 observações e quatro
  variáveis independentes (x1, x2, z1 e z2), geradas a partir de
  distribuições normal e uniforme.

- Utilizar a função `beta_ordinal_simula_dados_z` para gerar dados
  simulados com base nos parâmetros fornecidos, como fórmulas,
  coeficientes de regressão e funções de ligação.

- Calcular a log-verossimilhança do modelo ajustado com os dados
  simulados usando a função `beta_ordinal_log_vero_z`.

Como resultado, obtemos a log-verossimilhança do modelo ajustado aos
dados simulados, que é uma medida de quão bem o modelo se ajusta aos
dados observados no caso de um modelo beta ordinal com dispersão
variável.

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

A função `beta_ordinal_fit` realiza o ajuste de um modelo beta ordinal
por meio da função `optim` do pacote stats, retornando uma lista
contendo as informações essenciais do ajuste realizado.

Neste bloco de código R, ajustamos um modelo beta ordinal a partir de um
conjunto de dados simulados, utilizando as funções
`beta_ordinal_simula_dados` e `beta_ordinal_fit`. O processo é resumido
abaixo:

- Definir semente, tamanho da amostra, coeficientes de regressão,
  fórmula da relação entre as variáveis e parâmetro de dispersão.

- Criar um conjunto de dados de exemplo com 100 observações e três
  variáveis independentes (x1, x2 e x3), geradas a partir de
  distribuições normal e binomial.

- Utilizar a função beta_ordinal_simula_dados para gerar dados simulados
  com base nos parâmetros fornecidos, como fórmula, coeficientes de
  regressão e função de ligação.

- Ajustar o modelo beta ordinal aos dados simulados utilizando a função
  beta_ordinal_fit, especificando a fórmula, os dados, a função de
  ligação e a opção de cálculo da matriz hessiana numérica.

- Extrair os coeficientes estimados do ajuste do modelo usando a função
  beta_ordinal_coef.

Como resultado, obtemos os coeficientes estimados do modelo beta ordinal
ajustado aos dados simulados, que podem ser utilizados para análise e
interpretação das relações entre a variável resposta ordinal e as
variáveis explicativas.

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
  knitr::kable(digits = 4, caption = "")
```

| variable    | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept) |   0.2156 |   0.0353 |   0.3958 | 0.0920 |  2.3437 |  0.0212 |
| x1          |   0.3067 |   0.1965 |   0.4168 | 0.0562 |  5.4551 |  0.0000 |
| x2          |  -0.4558 |  -0.5734 |  -0.3381 | 0.0600 | -7.5940 |  0.0000 |
| x3          |   0.1003 |   0.0376 |   0.1630 | 0.0320 |  3.1357 |  0.0023 |
| (phi)       |  51.1296 |  36.7673 |  65.4919 | 7.3278 |  6.9774 |  0.0000 |

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|       330.5828 | -649.1657 | -633.5347 |

## Função para ajustar um modelo beta ordinal com dispersão variável

A função `beta_ordinal_fit_z` ajusta um modelo beta ordinal com
dispersão variável, isto é, com um preditor para os betas e outro para o
phi, onde covariáveis são aplicadas para explicar a dispersão. A função
usa optim do pacote stats e retorna uma lista com as principais
informações do ajuste.

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
  knitr::kable(digits = 4, caption = "")
```

| variable           | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:-------------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept)        |   0.1225 |  -0.1384 |   0.3833 | 0.1331 |  0.9203 |  0.3623 |
| x1                 |  -0.3995 |  -0.6869 |  -0.1121 | 0.1466 | -2.7247 |  0.0091 |
| x2                 |   0.3626 |   0.0953 |   0.6299 | 0.1364 |  2.6591 |  0.0108 |
| (phi)\_(Intercept) |   0.5855 |  -0.2018 |   1.3728 | 0.4017 |  1.4576 |  0.1519 |
| (phi)\_z1          |   1.3097 |  -0.0996 |   2.7189 | 0.7190 |  1.8214 |  0.0752 |

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|       222.5784 | -433.1567 | -421.6846 |

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
  knitr::kable(digits = 4, caption = "")
```

| variable           | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:-------------------|---------:|---------:|---------:|-------:|--------:|--------:|
| (Intercept)        |   0.2752 |   0.0765 |   0.4738 | 0.1013 |  2.7154 |  0.0094 |
| x1                 |  -0.5384 |  -0.7388 |  -0.3380 | 0.1023 | -5.2649 |  0.0000 |
| x2                 |   0.4141 |   0.2219 |   0.6064 | 0.0981 |  4.2227 |  0.0001 |
| (phi)\_(Intercept) |   1.6446 |   0.8457 |   2.4436 | 0.4076 |  4.0346 |  0.0002 |
| (phi)\_z1          |   0.7760 |  -0.5975 |   2.1495 | 0.7008 |  1.1074 |  0.2740 |

``` r

coe$gof %>% 
  knitr::kable(digits = 4, caption = "")
```

| log_likelihood |       AIC |       BIC |
|---------------:|----------:|----------:|
|        205.461 | -398.9221 | -387.4499 |

## Updating!
