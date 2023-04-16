
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betaregesc

<!-- badges: start -->

[![pkgdown](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/evandeilton/betaroti/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

O pacote `betaroti` oferece uma biblioteca de funções em R para ajuste
de modelos de regressão beta em dados ordinais transformados
intervalares, com dispersão fixa ou variável. Possibilita simulações
para avaliação do desempenho dos modelos no processo de estimação. O
código-fonte e contribuições podem ser acessados no repositório oficial
do GitHub. Informações detalhadas sobre instalação e uso estão
disponíveis na documentação do pacote.

O “betaroti” é voltado para modelagem de dados com variável resposta
ordinal numérica transformável em intervalo contínuo (e.g.,
$y = (y_s;y_i)$), abrangendo censura à esquerda, direita ou intervalar.
Aplica-se em pesquisas de opinião, avaliações de produtos, escalas de
dor, e avaliações de compostos químicos, entre outros. Utilizando a
distribuição beta, acomoda características dos dados em estrutura de
regressão, associando variáveis explicativas à variável resposta
intervalar e permitindo preditores lineares para coeficientes
relacionados à média e dispersão, fornecendo estimativas robustas e
confiáveis dos parâmetros do modelo.

## Principais funcionalidades

O pacote `betaroti` oferece uma série de funções úteis para lidar com
modelos de regressão beta e dados com resposta ordinal transformada
intervalar, abrangendo cenários com dispersão fixa e variável. As
principais funcionalidades incluem:

- Ajuste de modelos de regressão beta com dispersão fixa e variável.
- Funções para simulação de dados, permitindo a avaliação do desempenho
  dos modelos em diferentes cenários.
- Estatística de bondade do ajuste como AIC e BIC, por exemplo em
  `gof()`.
- Funções genéricas como `coef`, `vcov`, `fitted`, `residuals`,
  `summary` e `print` foram implementadas para a classe `betaroti` para
  facilitar o acesso às medidas do ajuste.
- Funções para ajuste e comparação de modelos com diferentes combinações
  de variáveis explicativas tanto para $\mu$ como $\phi$.

Acesse a documentação detalhada de cada função e exemplos de uso neste
site para obter informações sobre como utilizar o pacote “betaroti” em
suas análises.

## Instalação

Você pode instalar o pacote com esse comando abaixo.

``` r
if(!require(betaregesc)){
  devtools::install_github("evandeilton/betaregesc")  
}
require(betaregesc, quietly = TRUE)
```

## Exemplos

Esses são alguns exemplos de uso das funções do pacote.

### Simula dados do modelo beta ordinal com dispersão fixa

Esta função gera amostras de variável beta ordinal com dispersão fixa
usando várias funções de ligação.

No exemplo a seguir em código R, demonstramos como usar a função
beta_ordinal_simula_dados para simular dados de variável beta ordinal
com dispersão fixa:

- Criamos um conjunto de dados com 100 observações e duas variáveis
  explicativas independentes (x1 e x2) a partir de uma distribuição
  normal.
- Utilizamos a função beta_ordinal_simula_dados para simular dados com
  parâmetros personalizados fornecidos.

> OBS.: `type` é o tipo de tratamento do intervalo `m` centraliza `y` ao
> meio. Ex. Se foi coletado o valor $y = 6$, transforma-se
> $y_t = 6/10 = 0.6$. Assim, para tratar a incerteza do instrumento,
> sugere-se que a medida anotada pode estar limitada a $y_{left} = 5.5$
> e $y_{right} = 6.6$.

``` r
# Criar um conjunto de dados de exemplo
set.seed(4255)
dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

dados_simulados <- betaregesc_simula_dados(
  formula = ~ x1 + x2,
  dados = dados,
  betas = c(0.3, -0.6, 0.4),
  phi = 30,
  link = "probit",
  ncuts = 100,
  type = "m"
)
dados_simulados %>%
  head() %>%
  knitr::kable(digits = 4, caption = "")
```

|  left | right |   yt |   y |      x1 |      x2 |
|------:|------:|-----:|----:|--------:|--------:|
| 0.215 | 0.225 | 0.22 |  22 |  1.9510 |  0.6403 |
| 0.385 | 0.395 | 0.39 |  39 |  0.7725 | -0.2677 |
| 0.495 | 0.505 | 0.50 |  50 |  0.7264 |  0.0222 |
| 0.595 | 0.605 | 0.60 |  60 |  0.0487 |  0.0113 |
| 0.815 | 0.825 | 0.82 |  82 | -0.5445 |  0.3109 |
| 0.595 | 0.605 | 0.60 |  60 |  0.3600 |  0.4195 |

### Ajuste de modelos com dispersão fixa

- Exemplo do ajuste com optim direto para uma lista de links

``` r
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_fixo <- purrr::map(links, .f = function(link){
  betaregesc(
    formula = y ~ x1 + x2,
    dados = dados_simulados,
    link = link,
    link_phi = "sqrt",
    num_hessiana = TRUE)
})
```

- Resumo das estimativas e bondade

- Estimativas do ajuste e Bondade

``` r
resumo <- purrr::map(fit_fixo, function(fit){
  summary(fit)
})
```

``` r
purrr::map_df(resumo, function(res){
  res$est
  }, .id = "link") %>% 
  knitr::kable(digits = 4, caption = "")  
```

| link    | variable           | estimate | ci_lower | ci_upper |     se |  t_value | p_value |
|:--------|:-------------------|---------:|---------:|---------:|-------:|---------:|--------:|
| logit   | (Intercept)        |   0.4695 |   0.3906 |   0.5484 | 0.0403 |  11.6636 |       0 |
| logit   | x1                 |  -0.9485 |  -1.0447 |  -0.8522 | 0.0491 | -19.3154 |       0 |
| logit   | x2                 |   0.6897 |   0.6101 |   0.7692 | 0.0406 |  16.9906 |       0 |
| logit   | (phi)\_(Intercept) |   5.6981 |   4.9004 |   6.4958 | 0.4070 |  14.0000 |       0 |
| probit  | (Intercept)        |   0.2812 |   0.2341 |   0.3283 | 0.0240 |  11.7041 |       0 |
| probit  | x1                 |  -0.5695 |  -0.6256 |  -0.5134 | 0.0286 | -19.8981 |       0 |
| probit  | x2                 |   0.4116 |   0.3656 |   0.4577 | 0.0235 |  17.5212 |       0 |
| probit  | (phi)\_(Intercept) |   5.6762 |   4.8761 |   6.4764 | 0.4083 |  13.9033 |       0 |
| cloglog | (Intercept)        |  -0.1180 |  -0.1699 |  -0.0661 | 0.0265 |  -4.4557 |       0 |
| cloglog | x1                 |  -0.5929 |  -0.6545 |  -0.5314 | 0.0314 | -18.8795 |       0 |
| cloglog | x2                 |   0.4243 |   0.3736 |   0.4750 | 0.0259 |  16.4066 |       0 |
| cloglog | (phi)\_(Intercept) |   5.4515 |   4.6605 |   6.2424 | 0.4036 |  13.5084 |       0 |

``` r
purrr::map_df(resumo, function(res){
  res$gof
  }, .id = "link") %>% 
  knitr::kable(digits = 4, caption = "")
```

| link    |   logLik |       AIC |       BIC |
|:--------|---------:|----------:|----------:|
| logit   | 341.1732 | -672.3465 | -659.3206 |
| probit  | 341.9841 | -673.9683 | -660.9424 |
| cloglog | 347.7286 | -685.4573 | -672.4314 |

- Exemplo do ajuste com `bbmle` direto para uma lista de links

``` r
require(bbmle, quietly = TRUE)
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_fixo_bbmle <- purrr::map(links, .f = function(link){
  betaregesc_bbmle(
    formula = y ~ x1 + x2,
    dados = dados_simulados,
    link = link,
    link_phi = "sqrt",
    num_hessiana = TRUE)
})
```

- Gráficos dos perfis de verossimilhança

``` r
fit_fixo_profiles <- purrr::map(fit_fixo_bbmle, profile)
purrr::walk(names(fit_fixo_profiles), function(p){
  cat("\n+", p, "\n")
  plot(fit_fixo_profiles[[p]])
})
```

- logit
  <img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
- probit
  <img src="man/figures/README-unnamed-chunk-10-2.png" width="100%" />
- cloglog
  <img src="man/figures/README-unnamed-chunk-10-3.png" width="100%" />

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
set.seed(2222)
n <- 50
fx <- ~ x1 + x2 + x3
fz <- ~ z1 + z2

dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rbinom(n, size = 1, prob = 1/2),
  z1 = runif(n),
  z2 = runif(n)
)

dados_simulados <- betaregesc_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2, -0.6, 0.2, 0.2),
  zetas = c(0.5, 1, 2),
  link = "logit",
  link_phi = "sqrt",
  ncuts = 100,
  type = "m"
)

dados_simulados %>% 
  head() %>%
  knitr::kable(digits = 4, caption = "")
```

|  left | right |   yt |   y |      x1 |      x2 |  x3 |     z1 |     z2 |
|------:|------:|-----:|----:|--------:|--------:|----:|-------:|-------:|
| 0.745 | 0.755 | 0.75 |  75 | -0.3381 | -0.5606 |   0 | 0.9154 | 0.6295 |
| 0.545 | 0.555 | 0.55 |  55 |  0.9392 | -0.4519 |   1 | 0.3437 | 0.3239 |
| 0.145 | 0.155 | 0.15 |  15 |  1.7377 |  0.5993 |   1 | 0.5909 | 0.1956 |
| 0.435 | 0.445 | 0.44 |  44 |  0.6963 | -0.4836 |   0 | 0.7916 | 0.7053 |
| 0.275 | 0.285 | 0.28 |  28 |  0.4623 | -0.7956 |   1 | 0.9538 | 0.5312 |
| 0.485 | 0.495 | 0.49 |  49 | -0.3151 | -0.9410 |   0 | 0.2544 | 0.8612 |

### Ajuste de modelos com dispersão variável

- Exemplo do ajuste com optim direto para uma lista de links

``` r
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_variavel <- purrr::map(links, .f = function(link){
  betaregesc(
    formula = y ~x1 + x2 + x3 | z1 + z2,
    dados = dados_simulados,
    link = link,
    link_phi = "log",
    num_hessiana = TRUE)
})
```

- Resumo das estimativas e bondade

- Estimativas do ajuste e Bondade

``` r
resumo <- purrr::map(fit_variavel, function(fit){
  summary(fit)
})
```

``` r
purrr::map_df(resumo, function(res){
  res$est
  }, .id = "link") %>% 
  knitr::kable(digits = 4, caption = "")  
```

| link    | variable           | estimate | ci_lower | ci_upper |     se | t_value | p_value |
|:--------|:-------------------|---------:|---------:|---------:|-------:|--------:|--------:|
| logit   | (Intercept)        |   0.4605 |   0.2144 |   0.7066 | 0.1256 |  3.6673 |  0.0007 |
| logit   | x1                 |  -0.5781 |  -0.8006 |  -0.3556 | 0.1135 | -5.0921 |  0.0000 |
| logit   | x2                 |   0.0381 |  -0.1771 |   0.2534 | 0.1098 |  0.3474 |  0.7300 |
| logit   | x3                 |  -0.1979 |  -0.5819 |   0.1861 | 0.1959 | -1.0102 |  0.3181 |
| logit   | (phi)\_(Intercept) |  -0.3645 |  -1.4231 |   0.6942 | 0.5402 | -0.6747 |  0.5034 |
| logit   | (phi)\_z1          |   2.4141 |   1.0334 |   3.7949 | 0.7045 |  3.4269 |  0.0014 |
| logit   | (phi)\_z2          |   1.6193 |   0.0688 |   3.1697 | 0.7911 |  2.0469 |  0.0468 |
| probit  | (Intercept)        |   0.2855 |   0.1353 |   0.4356 | 0.0766 |  3.7261 |  0.0006 |
| probit  | x1                 |  -0.3536 |  -0.4846 |  -0.2226 | 0.0668 | -5.2915 |  0.0000 |
| probit  | x2                 |   0.0221 |  -0.1094 |   0.1535 | 0.0671 |  0.3291 |  0.7437 |
| probit  | x3                 |  -0.1245 |  -0.3609 |   0.1118 | 0.1206 | -1.0325 |  0.3076 |
| probit  | (phi)\_(Intercept) |  -0.3716 |  -1.4281 |   0.6848 | 0.5390 | -0.6894 |  0.4942 |
| probit  | (phi)\_z1          |   2.4347 |   1.0525 |   3.8170 | 0.7052 |  3.4524 |  0.0013 |
| probit  | (phi)\_z2          |   1.6180 |   0.0698 |   3.1661 | 0.7899 |  2.0484 |  0.0467 |
| cloglog | (Intercept)        |  -0.0760 |  -0.2299 |   0.0779 | 0.0785 | -0.9680 |  0.3385 |
| cloglog | x1                 |  -0.3754 |  -0.5074 |  -0.2435 | 0.0673 | -5.5783 |  0.0000 |
| cloglog | x2                 |   0.0233 |  -0.1090 |   0.1557 | 0.0675 |  0.3452 |  0.7316 |
| cloglog | x3                 |  -0.1213 |  -0.3722 |   0.1296 | 0.1280 | -0.9476 |  0.3486 |
| cloglog | (phi)\_(Intercept) |  -0.4400 |  -1.4888 |   0.6088 | 0.5351 | -0.8223 |  0.4155 |
| cloglog | (phi)\_z1          |   2.5081 |   1.1222 |   3.8939 | 0.7071 |  3.5470 |  0.0010 |
| cloglog | (phi)\_z2          |   1.7133 |   0.1706 |   3.2559 | 0.7871 |  2.1768 |  0.0350 |

``` r
purrr::map_df(resumo, function(res){
  res$gof
  }, .id = "link") %>% 
  knitr::kable(digits = 4, caption = "")
```

| link    |   logLik |       AIC |       BIC |
|:--------|---------:|----------:|----------:|
| logit   | 205.0974 | -394.1949 | -378.8987 |
| probit  | 205.0196 | -394.0392 | -378.7431 |
| cloglog | 204.5744 | -393.1489 | -377.8527 |

- Exemplo do ajuste com `bbmle` direto para uma lista de links

``` r
require(bbmle, quietly = TRUE)
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_variavel_bbmle <- purrr::map(links, .f = function(link){
  betaregesc_bbmle(
    formula = y ~ x1 + x2 + x3 | z1,
    dados = dados_simulados,
    link = link,
    link_phi = "sqrt",
    num_hessiana = TRUE)
})
```

- Gráficos dos perfis de verossimilhança

``` r
fit_variavel_profiles <- purrr::map(fit_variavel_bbmle, profile)
purrr::walk(names(fit_variavel_profiles), function(p){
  cat("\n+", p, "\n")
  plot(fit_variavel_profiles[[p]])
})
```

- logit
  <img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />
- probit
  <img src="man/figures/README-unnamed-chunk-17-2.png" width="100%" />
- cloglog
  <img src="man/figures/README-unnamed-chunk-17-3.png" width="100%" />

### Outras funções genéricas

``` r
## Resumo das estimativas e bondades
summary(fit_fixo$logit)

## Coeficientes do modelo
coef(fit_fixo$logit)

## Matriz de covariâncias
vcov(fit_fixo$logit)

## Resíduo dos valores preditos em relação ao ponto médio do intervalo de y, 
## isto é (left + right) / 2
resid(fit_fixo$logit)

## Valores preditos
fitted(fit_fixo$logit)

## Print do modelo
print(fit_fixo$logit)
```

------------------------------------------------------------------------

## Work in progress …

![](https://media4.giphy.com/media/LHZyixOnHwDDy/giphy.gif)
