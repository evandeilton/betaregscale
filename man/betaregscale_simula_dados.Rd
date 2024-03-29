% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaregscale.R
\name{betaregscale_simula_dados}
\alias{betaregscale_simula_dados}
\title{Simula dados do modelo beta ordinal com dispersão fixa}
\usage{
betaregscale_simula_dados(
  formula,
  dados,
  betas,
  phi = 1/5,
  link = "logit",
  link_phi = "logit",
  ncuts = 100,
  type = "m",
  lim = 0.5,
  repar = "2"
)
}
\arguments{
\item{formula}{Fórmula que indica a relação entre a variável dependente e as
variáveis independentes. O padrão é ~x1 + x2.}

\item{dados}{Um conjunto de dados que contém as variáveis independentes
especificadas na fórmula.}

\item{betas}{Vetor numérico de coeficientes de regressão. O padrão é c(0, 0.5, -0.2).}

\item{phi}{Parâmetro positivo que controla a dispersão da distribuição beta.
O padrão é 50.}

\item{link}{Nome da função de ligação a ser usada. Pode ser uma das
seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".}

\item{link_phi}{Nome da função de ligação a ser usada para transformar phi.}

\item{ncuts}{Número de cortes para a variável ordinal. O padrão é 100.}

\item{type}{Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.}

\item{lim}{Limite numérico a ser utilizado para ajustar os intervalos. O padrão é 0.5.}

\item{repar}{Tipo de reparametrização sendo 0, 1 ou 2. Padrão 2.}
}
\value{
Retorna um data.frame contendo os dados simulados da variável
 beta ordinal e as variáveis independentes.
}
\description{
Esta função simula dados de uma variável beta ordinal com dispersão fixa,
aplicando diferentes funções de ligação.
}
\examples{
# Criar um conjunto de dados de exemplo
set.seed(421)
n <- 100
dados <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
 )
formula <- ~ x1 + x2
dados_simulados <- betaregscale_simula_dados(
 formula = formula,
 dados = dados,
 betas = c(0.2,-0.5, 0.3),
 link = "logit",
 link_phi = "logit",
 ncuts = 100,
 phi = 1/2,
 repar = "2"
 )
}
