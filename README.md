
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betaroti

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
if(!require(betaroti)){
  devtools::install_github("evandeilton/betaroti")  
}
require(betaroti, quietly = TRUE)
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
set.seed(42)
dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

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

### Ajuste de modelos com dispersão fixa

- Exemplo do ajuste com optim direto para uma lista de links

``` r
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_fixo <- purrr::map(links, .f = function(link){
  betaroti(
    formula = ~x1 + x2,
    dados = dados_simulados,
    link = link,
    link_phi = "identity",
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

| link    | variable    | estimate | ci_lower | ci_upper |     se |  t_value | p_value |
|:--------|:------------|---------:|---------:|---------:|-------:|---------:|--------:|
| logit   | (Intercept) |   1.7168 |   1.6138 |   1.8199 | 0.0526 |  32.6556 |       0 |
| logit   | x1          |  -0.5225 |  -0.6149 |  -0.4302 | 0.0471 | -11.0865 |       0 |
| logit   | x2          |   0.7326 |   0.6262 |   0.8390 | 0.0543 |  13.4911 |       0 |
| logit   | (phi)       |  32.9500 |  23.5197 |  42.3804 | 4.8115 |   6.8482 |       0 |
| probit  | (Intercept) |   1.0091 |   0.9552 |   1.0630 | 0.0275 |  36.7154 |       0 |
| probit  | x1          |  -0.2893 |  -0.3395 |  -0.2391 | 0.0256 | -11.2963 |       0 |
| probit  | x2          |   0.4088 |   0.3512 |   0.4664 | 0.0294 |  13.9075 |       0 |
| probit  | (phi)       |  34.1445 |  24.3752 |  43.9138 | 4.9844 |   6.8502 |       0 |
| cloglog | (Intercept) |   0.5935 |   0.5500 |   0.6371 | 0.0222 |  26.7189 |       0 |
| cloglog | x1          |  -0.2367 |  -0.2786 |  -0.1947 | 0.0214 | -11.0477 |       0 |
| cloglog | x2          |   0.3384 |   0.2903 |   0.3865 | 0.0245 |  13.7909 |       0 |
| cloglog | (phi)       |  34.3315 |  24.3100 |  44.3529 | 5.1131 |   6.7144 |       0 |

``` r
purrr::map_df(resumo, function(res){
  res$gof
  }, .id = "link") %>% 
  knitr::kable(digits = 4, caption = "")
```

| link    |   logLik |       AIC |       BIC |
|:--------|---------:|----------:|----------:|
| logit   | 322.0392 | -634.0784 | -621.0526 |
| probit  | 321.3533 | -632.7067 | -619.6808 |
| cloglog | 322.9948 | -635.9895 | -622.9637 |

- Exemplo do ajuste com `bbmle` direto para uma lista de links

``` r
require(bbmle, quietly = TRUE)
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_fixo_bbmle <- purrr::map(links, .f = function(link){
  betaroti_bbmle(
    formula = ~x1 + x2,
    dados = dados_simulados,
    link = link,
    link_phi = "identity",
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

dados_simulados <- beta_ordinal_simula_dados_z(
  formula_x = fx,
  formula_z = fz,
  dados = dados,
  betas = c(0.2, -0.6, 0.2, 0.2),
  zetas = c(0.5, 1, 2),
  link_x = "logit",
  link_z = "sqrt",
  ncuts = 100,
  type = "m"
)

dados_simulados %>% 
  head() %>%
  knitr::kable(digits = 4, caption = "")
```

|      y |  yr |  left | right |      x1 |      x2 |  x3 |     z1 |     z2 |
|-------:|----:|------:|------:|--------:|--------:|----:|-------:|-------:|
| 0.7458 |  75 | 0.745 | 0.755 | -0.3381 | -0.5606 |   0 | 0.9154 | 0.6295 |
| 0.5467 |  55 | 0.545 | 0.555 |  0.9392 | -0.4519 |   1 | 0.3437 | 0.3239 |
| 0.1459 |  15 | 0.145 | 0.155 |  1.7377 |  0.5993 |   1 | 0.5909 | 0.1956 |
| 0.4424 |  44 | 0.435 | 0.445 |  0.6963 | -0.4836 |   0 | 0.7916 | 0.7053 |
| 0.2753 |  28 | 0.275 | 0.285 |  0.4623 | -0.7956 |   1 | 0.9538 | 0.5312 |
| 0.4916 |  49 | 0.485 | 0.495 | -0.3151 | -0.9410 |   0 | 0.2544 | 0.8612 |

### Ajuste de modelos com dispersão variável

- Exemplo do ajuste com optim direto para uma lista de links

``` r
links <- c("logit", "probit", "cloglog")
names(links) <- links

fit_variavel <- purrr::map(links, .f = function(link){
  betaroti(
    formula = ~x1 + x2 + x3 | z1 + z2,
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
| logit   | (Intercept)        |   0.4605 |   0.2144 |   0.7066 | 0.1256 |  3.6671 |  0.0007 |
| logit   | x1                 |  -0.5781 |  -0.8006 |  -0.3556 | 0.1135 | -5.0921 |  0.0000 |
| logit   | x2                 |   0.0381 |  -0.1771 |   0.2534 | 0.1098 |  0.3474 |  0.7300 |
| logit   | x3                 |  -0.1979 |  -0.5819 |   0.1861 | 0.1959 | -1.0100 |  0.3182 |
| logit   | (phi)\_(Intercept) |  -0.3644 |  -1.4231 |   0.6943 | 0.5401 | -0.6746 |  0.5035 |
| logit   | (phi)\_z1          |   2.4138 |   1.0331 |   3.7946 | 0.7045 |  3.4264 |  0.0014 |
| logit   | (phi)\_z2          |   1.6194 |   0.0689 |   3.1699 | 0.7911 |  2.0471 |  0.0468 |
| probit  | (Intercept)        |   0.2855 |   0.0368 |   0.5341 | 0.1269 |  2.2498 |  0.0296 |
| probit  | x1                 |  -0.3536 |  -0.5645 |  -0.1428 | 0.1076 | -3.2875 |  0.0020 |
| probit  | x2                 |   0.0221 |  -0.1824 |   0.2266 | 0.1043 |  0.2115 |  0.8335 |
| probit  | x3                 |  -0.1245 |  -0.5011 |   0.2520 | 0.1921 | -0.6482 |  0.5203 |
| probit  | (phi)\_(Intercept) |  -0.3717 |  -1.4069 |   0.6636 | 0.5282 | -0.7037 |  0.4854 |
| probit  | (phi)\_z1          |   2.4348 |   1.1570 |   3.7125 | 0.6519 |  3.7348 |  0.0005 |
| probit  | (phi)\_z2          |   1.6180 |   0.2073 |   3.0287 | 0.7198 |  2.2480 |  0.0298 |
| cloglog | (Intercept)        |  -0.0760 |  -0.4926 |   0.3406 | 0.2126 | -0.3575 |  0.7225 |
| cloglog | x1                 |  -0.3754 |  -0.5976 |  -0.1533 | 0.1134 | -3.3118 |  0.0019 |
| cloglog | x2                 |   0.0233 |  -0.1792 |   0.2258 | 0.1033 |  0.2256 |  0.8225 |
| cloglog | x3                 |  -0.1213 |  -0.5014 |   0.2588 | 0.1939 | -0.6255 |  0.5349 |
| cloglog | (phi)\_(Intercept) |  -0.4400 |  -1.3810 |   0.5010 | 0.4801 | -0.9164 |  0.3646 |
| cloglog | (phi)\_z1          |   2.5081 |   1.2985 |   3.7176 | 0.6171 |  4.0641 |  0.0002 |
| cloglog | (phi)\_z2          |   1.7133 |   0.4721 |   2.9544 | 0.6333 |  2.7054 |  0.0097 |

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
  betaroti_bbmle(
    formula = ~x1 + x2 + x3 | z1,
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
