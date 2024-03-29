% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaregscale.R
\name{fn_switch_link}
\alias{fn_switch_link}
\title{Seleção das funções de ligação inversas}
\usage{
fn_switch_link(eta, link)
}
\arguments{
\item{eta}{Preditor linear dado por X*param.}

\item{link}{Função de ligação escolhida.}
}
\value{
Valor de mu, resultante da aplicação da transformação inversa em eta.
}
\description{
Função para escolher e aplicar as funções de ligação inversas aos preditores de X e Z
}
