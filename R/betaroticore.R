## -------------------------------------------------------------------------- ##
## titulo:  REGRESSÃO BETA PARA DADOS DE ESCALA
##          INTERVALAR
## Autor:   José Lopes / Wagner Bonatl
## Code:    Simulação da beta
## -------------------------------------------------------------------------- ##

#' Seleção das funções de ligação inversas
#'
#' Função para escolher e aplicar as funções de ligação inversas aos preditores de X e Z
#' @param eta Preditor linear dado por X*param.
#' @param link Função de ligação escolhida.
#' @return Valor de mu, resultante da aplicação da transformação inversa em eta.
fn_switch_link <- function(eta, link){
  switch(link,
         logit = make.link("logit")$linkinv(eta),
         probit = make.link("probit")$linkinv(eta),
         cauchit = make.link("cauchit")$linkinv(eta),
         cloglog = make.link("cloglog")$linkinv(eta),
         log = make.link("log")$linkinv(eta),
         sqrt = make.link("sqrt")$linkinv(eta),
         identity = eta
  )
}

#' Simula dados do modelo beta ordinal com dispersão fixa
#'
#' Esta função simula dados de uma variável beta ordinal com dispersão fixa,
#' aplicando diferentes funções de ligação.
#'
#' @param formula Fórmula que indica a relação entre a variável dependente e as
#'  variáveis independentes. O padrão é ~x1 + x2.
#' @param dados Um conjunto de dados que contém as variáveis independentes
#' especificadas na fórmula.
#' @param betas Vetor numérico de coeficientes de regressão. O padrão é c(0, 0.5, -0.2).
#' @param phi Parâmetro positivo que controla a dispersão da distribuição beta.
#'  O padrão é 50.
#' @param link Nome da função de ligação a ser usada. Pode ser uma das
#'  seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Limite numérico a ser utilizado para ajustar os intervalos. O padrão é 0.5.
#' @return Retorna um data.frame contendo os dados simulados da variável
#'  beta ordinal e as variáveis independentes.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(42)
#' dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#'
#' # Simular dados usando a função betaregesc_simula_dados com parâmetros personalizados
#' dados_simulados <- betaregesc_simula_dados(
#'  formula = ~x1 + x2, dados = dados,
#'  betas = c(1, -0.3, 0.4), phi = 30,
#'  link = "logit", ncuts = 100, type = "m")
#' head(dados_simulados)
#' @export
betaregesc_simula_dados <- function(formula, dados, betas, phi = 50, link = "logit", 
                                    ncuts = 100, type = "m", lim = 0.5){

  mfx <- model.frame(formula, data = dados)
  X   <- model.matrix(mfx, data = dados)
  n   <- nrow(dados)

  # Aplicação do preditor linear nos betas
  eta <- X%*%betas
  mu <- fn_switch_link(eta = eta, link = link)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu * phi)
  beta <- as.numeric((1 - mu)*phi)
  y <- rbeta(n = n, shape1 = mu*phi, shape2 = (1-mu)*phi)
  y_meio <- round(y*ncuts, 0)

  out_y <- fn_check_response(y = y_meio, type = type, ncuts = ncuts, lim = lim)

  out <- data.frame(out_y, X[,-1, drop = FALSE])
  return(out)
}

#' Simula dados do modelo beta ordinal com dispersão variável
#'
#' Esta função simula dados de uma variável beta ordinal, aplicando diferentes
#' funções de ligação tanto para betas como para zetas de phi
#'
#' @param formula_x Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas.
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param formula_z Fórmula para expressar a relação das preditoras Z1, Z2, Zn
#' relacionadas com os phi's no modelo se regressão beta com preditoas em phi.
#' Ela deve ser referenciada em Y. Ex. formula = ~Z1 + Z2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param link Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity".
#'  O padrão é "logit".
#' @param link_phi Nome da função de ligação a ser usada para as preditoras
#' Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param dados Um conjunto de dados que contém a variável dependente e as
#'  variáveis independentes especificadas nas fórmulas.
#' @param betas Vetor de betas associados aos preditores de X utilizados para simular da beta.
#' @param zetas Vetor de zetas associados aos preditores de Z utilizados para simular da beta com dispersão variável.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Limite numérico a ser utilizado para ajustar os intervalos. O padrão é 0.5.
#' @return Retorna um data.frame contendo os dados simulados da variável beta
#' ordinal e as variáveis independentes.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(421)
#' n <- 50
#' dados <- data.frame(
#' x1 = rnorm(n), x2 = rnorm(n),
#' z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- betaregesc_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link = "logit",
#' link_phi = "log",
#' ncuts = 100,
#' type = "m")
#' @export
betaregesc_simula_dados_z <- function(formula_x = ~x1 + x2,
                                      formula_z = ~z1 + z2,
                                      dados,
                                      betas = c(0, 0.5, -0.2),
                                      zetas = c(1, 0.5, 0.2),
                                      link = "logit",
                                      link_phi = "log",
                                      ncuts = 100,
                                      type = "m", 
                                      lim = 0.5){
  mfx <- model.frame(formula_x, data = dados)
  mfz <- model.frame(formula_z, data = dados)
  X   <- model.matrix(mfx, data = dados)
  Z   <- model.matrix(mfz, data = dados)
  n   <- nrow(X)

  link <- match.arg(link, c("logit","probit","cauchit","cloglog","identity"))
  link_phi <- match.arg(link_phi, c("log","sqrt","identity"))

  # Aplicação do preditor linear nos betas
  mu_x <- fn_switch_link(eta = X%*%betas, link = link)
  mu_p <- fn_switch_link(eta = Z%*%zetas, link = link_phi)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu_x * mu_p)
  beta <- as.numeric((1 - mu_x)*mu_p)
  y <- rbeta(n = n, shape1 = mu_x*mu_p, shape2 = (1-mu_x)*mu_p)
  y_meio <- round(y*ncuts, 0)
  
  out_y <- fn_check_response(y = y_meio, type = type, ncuts = ncuts, lim = lim)

  out <- data.frame(out_y, X[,-1, drop = FALSE], Z[,-1, drop = FALSE])
  return(out)
}

#' Log-verossimilhança do modelo beta ordinal com dispersão fixa
#'
#' Esta função calcula a log-verossimilhança de um conjunto de dados para uma
#' variável beta ordinal, aplicando diferentes funções de ligação.
#'
#' @param param Vetor numérico de parâmetros, incluindo coeficientes de
#' regressão e o parâmetro phi.
#' @param formula Fórmula que indica a relação entre a variável dependente
#' e as variáveis independentes.
#' @param dados Um conjunto de dados que contém a variável dependente e as
#' variáveis independentes especificadas na fórmula.
#' @param link Nome da função de ligação a ser usada. Pode ser uma das
#' seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity".
#' O padrão é "logit".
#' @param link_phi Função de ligação para phi. Uma dentre "log","sqrt","identity"
#' @param acumulada Um valor lógico indicando se a log-verossimilhança
#' acumulada deve ser calculada. O padrão é TRUE.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @return Retorna a soma da log-verossimilhança dos dados.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(421)
#' dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' # Calcular a log-verossimilhança usando a função log_vero_beta_ordinal
#' param <- c(0, 0.5, -0.2, 50)
#' phi <- 30
#' formula <- y ~ x1 + x2
#' dados_simulados <- betaregesc_simula_dados(formula = ~ x1 + x2, dados = dados,
#' betas = c(0, 0.5, -0.2), phi = phi, link = "logit", ncuts = 100, type = "m")
#' log_verossimilhanca <- betaregesc_log_vero(param, formula, dados_simulados)
#' print(log_verossimilhanca)
#' @importFrom stats as.formula delete.response model.response
#' @export
betaregesc_log_vero <- function(param, formula, dados, link = "logit", link_phi = "log",
                                acumulada = TRUE, ncuts = 100, type = "m", lim = 0.5
                                ){
  link <- match.arg(link, c("logit","probit","cauchit","cloglog"))
  link_phi <- match.arg(link_phi, c("log","sqrt","identity"))
  
  mfx <- model.frame(formula, data = dados)
  Y <- fn_check_response(model.response(mfx), ncuts = ncuts, type = type, lim = lim)
  X <- model.matrix(mfx, data = dados)
  
  betas  <- param[1:ncol(X)]

  # Aplicação do preditor linear nos betas
  eta <- X%*%betas
  mu <- fn_switch_link(eta = eta, link = link)
  phi <- fn_switch_link(eta = param[(length(betas)+1)], link = link_phi)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu * phi)
  beta <- as.numeric((1 - mu)*phi)

  ll <- if(acumulada){
    p1 <- pbeta(q = as.numeric(Y[,"left"]), shape1 = alpha, shape2 = beta)
    p2 <- pbeta(q = as.numeric(Y[,"right"]), shape1 = alpha, shape2 = beta)
    area <- p2 - p1
    area <- area + 0.00001
    log(area)
  } else {
    suppressWarnings(dbeta(dados[,"y"], shape1 = alpha, shape2 = beta, log = TRUE))
  }

  ll[is.infinite(ll)] <- NaN
  return(sum(ll, na.rm = TRUE))
}


#' Log-verossimilhança do modelo beta ordinal com dispersão variável
#'
#' Esta função calcula a log-verossimilhança de um conjunto de dados para uma
#' variável beta ordinal, aplicando diferentes funções de ligação tanto no betas
#' de mu como no zetas de phi.
#'
#' @param param Vetor numérico de parâmetros, incluindo coeficientes de
#' regressão e o parâmetro phi.
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas. Ex. formula = y ~ X1 + X2, y ~ X1 + X2|z1, y ~ X1 + X2 | z1 + z2.
#' @param link Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit" ou "cloglog".
#' O padrão é "logit".
#' @param link_phi Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn
#'  relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#' independentes especificadas na fórmula.
#' @param acumulada Um valor lógico indicando se a log-verossimilhança acumulada
#' deve ser calculada. O padrão é TRUE.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @return Retorna a soma da log-verossimilhança dos dados.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' n <- 50
#' dados <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
#'                     z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- betaregesc_simula_dados_z(
#'   formula_x = fx,
#'   formula_z = fz,
#'   dados = dados,
#'   betas = c(0.2, -0.5, 0.3),
#'   zetas = c(1, 1.2),
#'   link = "logit",
#'   link_phi = "log",
#'   ncuts = 100)
#' # Calcular a log-verossimilhança usando a função betaregesc_log_vero_z
#' log_verossimilhanca <- betaregesc_log_vero_z(
#'   param = c(c(0.2, -0.5, 0.3), c(1, 1.2)),
#'   formula = y ~ x1 + x2 | z1,
#'   dados = dados_simulados,
#'   link = "logit",
#'   link_phi = "log",
#'   acumulada = TRUE)
#'   print(log_verossimilhanca)
#' @export
betaregesc_log_vero_z <- function(param,
                                  formula = y ~x1 + x2 | z1, dados, 
                                  link = "logit", link_phi = "log",
                                  ncuts = 100, type = "m", lim = 0.5,
                                  acumulada = TRUE){
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  } else if(length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }
  
  mf <- model.frame(formula, data = dados)
  
  mtX <- terms(formula, data = dados, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = dados, rhs = 2L))
  Y <- fn_check_response(model.response(mf, "numeric"), ncuts = ncuts, type = type, lim = lim)
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  n   <- nrow(X)

  link_x <- match.arg(link, c("logit","probit","cauchit","cloglog"))
  link_z <- match.arg(link_phi, c("log","sqrt","identity"))

  betas  <- param[1:ncol(X)]
  zetas  <- param[(ncol(X)+1):length(param)]
  names(zetas) <- paste0("phi", names(zetas))

  # Aplicação do preditor linear nos betas
  mu_x <- fn_switch_link(eta = X%*%betas, link = link_x)
  mu_p <- fn_switch_link(eta = Z%*%zetas, link = link_z)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu_x * mu_p)
  beta <- as.numeric((1 - mu_x)*mu_p)

  ll <- if(acumulada){
    p1 <- pbeta(q = as.numeric(Y[,"left"]), shape1 = alpha, shape2 = beta)
    p2 <- pbeta(q = as.numeric(Y[,"right"]), shape1 = alpha, shape2 = beta)
    area <- p2 - p1
    area <- area + 0.00001
    log(area)
  } else {
    suppressWarnings(dbeta(Y[,"yt"], shape1 = alpha, shape2 = beta, log = TRUE))
  }

  ll[is.infinite(ll)] <- NaN
  return(sum(ll, na.rm = TRUE))
}

#' Função para ajustar um modelo beta 
#'
#' A função fit_beta_ ajusta um modelo beta usando a função optim
#' do pacote stats, retornando uma tabela com estatísticas sumarizadas do ajuste,
#' incluindo intervalos de confiança e cálculo do BIAS.
#'
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn.
#' Ela deve ser referenciada em Y. Ex. formula = y ~ X1 + X2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula.
#' @param link Nome da função de ligação a ser usada. Pode ser uma das
#' seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_phi Função de ligação para phi. Uma dentre "log","sqrt","identity"
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @param acumulada Se TRUE, retorna a verossimilhança pela pbeta, dbeta caso contrário.
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'  x1 = rnorm(n, mean = 1, sd = 0.5),
#'  x2 = rbinom(n, size = 1, prob = 0.5),
#'  x3 = rnorm(n, mean = 2, sd = 1))
#'  betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- y ~ x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'  formula = ~ x1 + x2 + x3,
#'  dados = dados,
#'  betas = betas,
#'  phi = phi,
#'  link = "logit",
#'  ncuts = 100)
#' fit <- betaregesc_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#'  fit$par
#' @importFrom betareg betareg
#' @importFrom numDeriv hessian
#' @importFrom stats AIC cor fitted hatvalues logLik qlogis runif terms
#' @export
betaregesc_fit <- function(formula, dados, link = "logit", link_phi = "identity",
                           ncuts = 100, type = "m", lim = 0.5, num_hessiana = TRUE, acumulada = TRUE) {

  ini <- coef(betareg::betareg(formula = as.formula(paste0("yt ~", as.character(formula)[3])),
                               data = dados, link = link, link.phi = link_phi))

  # Ajustando o modelo com a função optim
  opt_result <- optim(par = ini, 
                      fn = betaregesc_log_vero, formula = formula, 
                      dados = dados, link = link, link_phi = link_phi,
                      hessian = !num_hessiana, ncuts = ncuts, type = type, lim = lim,
                      method = "BFGS", acumulada = acumulada,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(betaregesc_log_vero, opt_result$par, 
                                  formula = formula, dados = dados, acumulada = acumulada,
                                  ncuts = ncuts, type = type, lim = lim,
                                  link = link, link_phi = link_phi)
    opt_result$hessian <- hessiana
  }

  # Mu chapéu
  est <- opt_result$par
  
  mf <- model.frame(formula, data = dados)
  mtX <- terms(formula, data = dados, rhs = 1L)
  Y <- fn_check_response(model.response(mf, "numeric"), ncuts = ncuts, type = type, lim = lim)
  y <- Y[,c("left","right")]
  X <- model.matrix(mtX, mf)
  hatmu <- fn_switch_link(eta = X%*%est[1:ncol(X)], link = link)
  hatphi <- est[ncol(X)+1]
  pseudor2 <- cor(X%*%est[1:ncol(X)], make.link(link)$linkfun(apply(y, 1, mean)))^2
  
  opt_result$dados <- data.frame(Y, X[,-1, drop = FALSE])
  opt_result$link <- link
  opt_result$link_phi <- link_phi
  opt_result$formula_x <- formula
  opt_result$formula_z <- ~ 1
  
  opt_result$residuals <- apply(y, 1, mean) - hatmu
  opt_result$hatmu  <- hatmu
  opt_result$hatphi <- hatphi
  opt_result$pseudo.r.squared = pseudor2
  
  class(opt_result) <- c("betaregesc","betaregescdv")
  
  return(invisible(opt_result))
}



#' Função para ajustar um modelo beta ordinal
#'
#' A função betaregesc ajusta um modelo beta para escala usando a função optim
#' do pacote stats, retornando uma tabela com estatísticas sumarizadas do ajuste,
#' incluindo intervalos de confiança e cálculo do BIAS.
#'
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas e Z1 + Z2 com phi. A formula pode ser escrita em
#' forma composta. Isto é, com três partes, como no exemplo y ~ x1 + x2 | z1 + z2.
#' Para mais detalhes veja \code{\link{Formula}}.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula.
#' @param link Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_phi Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @param acumulada Se TRUE, retorna a verossimilhança pela pbeta, dbeta caso contrário.
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @importFrom betareg betareg
#' @importFrom numDeriv hessian
#' @examples
#' n <- 50
#' dados <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n), z2 = runif(n))
#' 
#' formula <- y ~ x1 + x2 | z1
#' dados_simulados <- betaregesc_simula_dados_z(
#'  formula_x = ~x1 + x2,
#'  formula_z = ~z1,
#'  dados = dados,
#'  betas = c(0.2, -0.5, 0.3),
#'  zetas = c(1, 1.2),
#'  link = "logit",
#'  link_phi = "log",
#'  ncuts = 100)
#'  
#' fit_z <- betaregesc_fit_z(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "logit",
#'  link_phi = "log",
#'  num_hessiana = TRUE)
#' coe <- betaregesc_coef(fit_z)
#' coe$est
#' coe$gof
#' @export
betaregesc_fit_z <- function(formula,
                             dados,
                             link = "logit",
                             link_phi = "log",
                             num_hessiana = TRUE, 
                             acumulada = TRUE,
                             ncuts = 100,
                             type = "m",
                             lim = 0.5
                             ) {
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  } else if(length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }
  
  ini <- c(coef(betareg::betareg(formula = as.formula(paste0("yt ~", as.character(formula)[3])), data = dados)))

  # Ajustando o modelo com a função optim
  opt_result <- optim(par = ini,
                      fn = betaregesc_log_vero_z,
                      formula = formula, 
                      link = link,
                      link_phi = link_phi,
                      dados = dados,
                      hessian = !num_hessiana,
                      method = "BFGS",
                      acumulada = acumulada,
                      ncuts = ncuts, type = type, lim = lim,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(betaregesc_log_vero_z,
                                  opt_result$par,
                                  ncuts = ncuts, type = type, lim = lim,
                                  formula = formula, link = link, link_phi = link_phi, dados = dados)
    opt_result$hessian <- hessiana
  }
  # Mu chapéu
  est <- opt_result$par
  
  mf <- model.frame(formula, data = dados)
  mtX <- terms(formula, data = dados, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = dados, rhs = 2L))
  Y <- fn_check_response(model.response(mf, "numeric"), ncuts = ncuts, type = type, lim = lim)
  y <- Y[,c("left","right")]
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  p <- 1:ncol(X)
  q <- (ncol(X)+1):length(est)
  hatmu  <- fn_switch_link(eta = X%*%est[p], link = link)
  hatphi <- fn_switch_link(eta = Z%*%est[q], link = link_phi)
  pseudor2 <- cor(X%*%est[p], make.link(link)$linkfun(apply(y, 1, mean)))^2

  opt_result$dados <- data.frame(Y, X[,-1, drop = FALSE], Z[,-1, drop = FALSE])
  opt_result$link <- link
  opt_result$link_phi <- link_phi
  opt_result$formula <- formula
  opt_result$residuals <- apply(y, 1, mean) - hatmu
  opt_result$hatmu  <- hatmu
  opt_result$hatphi <- hatphi
  opt_result$pseudo.r.squared = pseudor2
  
  class(opt_result) <- c("betaregesc","betaregescdv")
  
  return(invisible(opt_result))
}


#' Coleta estatística do ajuste
#'
#' Coleta as estimativas e suas estatísticas para um objeto de ajuste
#' do modelo beta ordinal com censura intervalar. Coleta também as estatísticas
#' de bondade do ajuste como a log-verossimilhança, o AIC e o BIC do modelo
#'
#' @param fit Objeto do sjuste retornado das funções \link{betaregesc_fit} e
#' \link{betaregesc_fit_z}
#' @param alpha Nível de significância do alpha para os intervalos de confiança.
#' Padrão 0.05 bilateral.
#'
#' @examples
#' # Simulando dados
#' n <- 50
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   z1 = runif(n),
#'   z2 = runif(n)
#'  )
#'  dados_simulados <- betaregesc_simula_dados_z(
#'    formula_x = fx,
#'    formula_z = fz,
#'    dados = dados,
#'    betas = c(0.2,-0.5, 0.3),
#'    zetas = c(1, 1.2),
#'    link = "logit",
#'    link_phi = "log",
#'    ncuts = 100
#'    )
#'
#'  fit_z <- betaregesc_fit_z(
#'   formula = y ~ x1 + x2,
#'   dados = dados_simulados,
#'   link = "logit",
#'   link_phi = "log",
#'   num_hessiana = TRUE
#'   )
#'
#' coe <- betaregesc_coef(fit_z)
#' coe$est
#' coe$gof
#' @return Lista contendo as estimativas coeficientes, suas estatísticas e a bondade do ajuste.
#' @export
betaregesc_coef <- function(fit, alpha = 0.05){
  
  if(!inherits(fit, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(fit), collapse = ","), "' nao suportada.\n"))
  }
  
  dados <- fit$dados
  # Obtendo os parâmetros ajustados
  beta_hat <- fit$par
  # Criando uma tabela com as estatísticas sumarizadas do ajuste
  summary_table <- data.frame(
    #beta = beta_hat,
    logLik = -fit$value,
    AIC = 2 * (length(beta_hat) + 1) - 2 * (-fit$value),
    BIC = log(nrow(dados)) * (length(beta_hat) + 1) - 2 * (-fit$value)
  )

  # Obtendo os intervalos de confiança para as estimativas dos coeficientes
  cov_mat <- -solve(fit$hessian)
  se_beta <- sqrt(diag(cov_mat))
  z_alpha <- qnorm(1 - alpha/2)
  ci_lower <- beta_hat - z_alpha * se_beta
  ci_upper <- beta_hat + z_alpha * se_beta
  ci_table <- data.frame(
    variable = names(beta_hat),
    estimate = beta_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    se = se_beta,
    t_value = beta_hat/se_beta,
    p_value = 2 * pt(-abs(beta_hat/se_beta), df = nrow(dados) - length(beta_hat))
  )
  row.names(ci_table) <- NULL
  # Retornando o resultado
  list(
    gof = summary_table,
    est = ci_table
  )
}

#' @title Log-verossimilhança
#'
#' @description
#' Esta função calcula a log-verossimilhança negativa para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param ... Argumentos extras
#' 
#' @examples
#' \dontrun{
#' # Exemplo de uso da função logLik
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' # Em seguida, use a função logLik.betaregesc
#' log_likelihood <- logLik(fit)
#'}
#' @method logLik betaregesc
#' @export
logLik.betaregesc <- function(object, ...){
  
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  -object$value
}


#' @title Critério de Informação de Akaike (AIC)
#'
#' @description
#' Esta função calcula o Critério de Informação de Akaike (AIC) para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param ... Argumentos extras
#' @param k Penalidade por parâmetro a ser utilizado; o padrão k = 2 é o AIC clássico.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função AIC
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' # Em seguida, use a função AIC
#' aic_result <- AIC(fit)
#'}
#' @method AIC betaregesc
#' @export
AIC.betaregesc <- function(object, ..., k = 2){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  k * (length(object$par) + 1) - 2 * (-object$value)
}


#' @title Critério de Informação Bayesiano (BIC)
#'
#' @description
#' Esta função calcula o Critério de Informação Bayesiano (BIC) para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função BIC
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' # Em seguida, use a função BIC
#' bic_result <- BIC(fit)
#' }
#' @export
BIC <- function(object){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  log(nrow(object$dados)) * (length(object$par) + 1) - 2 * (-object$value)
}

#' @title Matriz hessiana
#'
#' @description
#' Esta função retorna a matriz hessiana para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' # Em seguida, use a função hessian
#' hessian_result <- hessian(fit)
#' }
#' @export
hessian <- function(object){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  object$hessian
}

#' Coeficientes
#' @description
#' Esta função retorna os coeficientes para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param ... outros argumentos.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' coef(fit)
#' }
#' @method coef betaregesc
#' @export
coef.betaregesc <- function(object, ...){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  object$par
}

#' @title Medidas de ajuste
#'
#' @description
#' Esta função retorna um conjunto de medidas de ajuste, incluindo log-verossimilhança, AIC e BIC, para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' gof(fit)
#' }
#' @return Retorna um data.frame contendo log-verossimilhança, AIC e BIC do objeto fornecido.
#' 
#' @export
gof <- function(object){
  data.frame(
    logLik = logLik(object),
    AIC = AIC(object),
    BIC = betaregesc::BIC(object)
  )
}

#' @title Estimativas e intervalos de confiança
#'
#' @description
#' Esta função retorna estimativas, erros padrão, intervalos de confiança e valores-p para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param alpha Nível de significância para os intervalos de confiança (padrão é 0,05).
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' est(fit)
#' }
#' @return
#' Retorna um data.frame contendo estimativas, erros padrão, intervalos de confiança, estatísticas t e valores-p do objeto fornecido.
#' @export
est <- function(object, alpha = 0.05){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  beta_hat <- object$par
  dfs <- nrow(object$dados) - length(beta_hat)
  cov_mat <- -solve(object$hessian)
  se_beta <- sqrt(diag(cov_mat))
  z_alpha <- qnorm(1 - alpha/2)
  ci_lower <- beta_hat - z_alpha * se_beta
  ci_upper <- beta_hat + z_alpha * se_beta
  
  ci_table <- data.frame(
    variable = names(beta_hat),
    estimate = beta_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    se = se_beta,
    t_value = beta_hat/se_beta,
    p_value = 2 * pt(-abs(beta_hat/se_beta), df = dfs)
  )
  row.names(ci_table) <- NULL
  return(ci_table)
}

#' @title Summary
#'
#' @description
#' Esta função retorna um resumo das estimativas, erros padrão, intervalos de confiança e valores-p para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param alpha Nível de significância para os intervalos de confiança (padrão é 0,05).
#' @param ... Outros argumentos
#'
#' @return
#' Retorna um resumo das estimativas, erros padrão, intervalos de confiança e valores-p do objeto fornecido.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' # Em seguida, use a função summary.betaregesc
#' summary_result <- summary.betaregesc(fit, alpha = 0.05)
#'}
#' @method summary betaregesc
#' @export
summary.betaregesc <- function(object, ..., alpha = 0.05){
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  betaregesc_coef(object, alpha = alpha)
}

#' @title Resíduos do modelo
#'
#' @description
#' Esta função retorna os resíduos para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param type Tipo de resíduo. Um entre "rqa","deviance", "pearson", "response", "weighted", "sweighted". O padrão é Resíduo Quantílico Aleatorizado (rqa)
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' residuals(fit)
#' }
#' @method residuals betaregesc
#' @export
residuals.betaregesc <- function(object, type = c("rqa","deviance", "pearson", "response", "weighted", "sweighted"), ...){
  type <- match.arg(type)
  
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  
  y  <- apply(as.matrix(object$dados[,c("left","right")]), 1, mean)
  mu <- object$hatmu
  phi <- object$hatphi
  res <- object$residuals
  wts <- 1
  
  if(type == "response"){
	return(res)
  }

  res <- switch(type,

    "pearson" = {
      sqrt(wts) * res / sqrt(mu * (1 - mu) / (1 + phi))
    },

    "deviance" = {
      ll <- function(mu, phi){
        (lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
        (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y))	  
	  }		
      sqrt(wts) * sign(res) * sqrt(2 * abs(ll(y, phi) - ll(mu, phi)))
    },
	
	"rqa" = {
      ll <- function(mu, phi){
	    p <- mu * phi
		q <- (1 - mu) * phi
		rqa <- qbeta(runif(length(y)), p, q)
		qnorm(rqa)
	  }
		ll(mu, phi)
    },

    "weighted" = {
      ystar <- qlogis(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(phi * v)
    },

    "sweighted" = {
      ystar <- qlogis(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(v)
    }
	)

	return(res)
}


#' @title Valores ajustados
#'
#' @description
#' Esta função retorna os valores ajustados (mu ou phi) para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param type Tipo de valor ajustado a ser retornado: "mu" ou "phi" (padrão é "mu").
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' fitted(fit)
#' }
#' @method fitted betaregesc
#' @export
fitted.betaregesc <- function(object, type = "mu", ...){
  type <- match.arg(type, c("mu","phi"))
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  if(type == "mu"){
    return(object$hatmu)
  } else {
    return(object$hatphi)
  }
}

#' @title Quantis preditos
#'
#' @description
#' Esta função retorna as previsões de quantis para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param type Tipo de previsão a ser retornado: "quantile" (padrão é "quantile").
#' @param at Valor numérico entre 0 e 1 para o qual o quantil é desejado (padrão é 0,5).
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' predict(fit)
#' }
#' @importFrom stats qbeta
#' @method predict betaregesc
#' @export
predict.betaregesc <- function(object, type = "quantile", at = 0.5, ...){
  type <- match.arg(type, c("quantile"))
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  
  qfun <- function(at, mu, phi) {
    rval <- sapply(at, function(p) stats::qbeta(p, mu * phi, (1 - mu) * phi))
    if(length(at) > 1L) {
      if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
                                          dimnames = list(unique(names(rval)), NULL))
      colnames(rval) <- paste("q_", at, sep = "")
    } else {
      rval <- drop(rval)
    }
    rval   
  }
  if(type == "quantile") {
    mu  <- fitted(object, type = "mu")
    phi <- fitted(object, type = "phi")
    return(qfun(at, mu, phi))
  }
}


#' @title Matriz de covariância
#'
#' @description
#' Esta função retorna a matriz de covariância para objetos das classes 'betaregesc' e 'betaregescdv'.
#'
#' @param object Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param ... Outros parametros
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' vcov(fit)
#' }
#' @method vcov betaregesc
#' @export
vcov.betaregesc <- function(object, ...){
  
  if(!inherits(object, c("betaregesc","betaregescdv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaregesc' ou 'betaregescdv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  return(solve(-as.matrix(object$hessian)))
}

#' @title Print
#'
#' @description
#' Esta função imprime um resumo do objeto das classes 'betaregesc' e 'betaregescdv' fornecido, incluindo estimativas, erros padrão, valores t e valores-p.
#'
#' @param x Um objeto das classes 'betaregesc' ou 'betaregescdv'.
#' @param digits Número de dígitos significativos para a impressão dos valores (padrão é max(3, getOption("digits") - 2)).
#' @param ... outros parametros
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaregesc' ou 'betaregescdv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- betaregesc_simula_dados(
#'   formula = formula,
#'   dados = dados,
#'   betas = betas,
#'   phi = phi,
#'   link = "logit",
#'   ncuts = 100)
#' fit <- beta_ordinal_fit(
#'  formula = formula,
#'  dados = dados_simulados,
#'  link = "probit",
#'  num_hessiana = TRUE)
#' print(fit)
#' }
#' @importFrom stats printCoefmat qbeta
#' @method print betaregesc
#' @export
print.betaregesc <- function(x, digits = max(3, getOption("digits") - 2), ...){
  cf <- betaregesc_coef(x)
  be <- cf$est[,c("estimate","se","t_value","p_value")]
  colnames(be) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(be) <- cf$est$variable
  stats::printCoefmat(be, digits = digits, P.values = TRUE)
  cat("\nLog-Likelihood:", formatC(logLik(x), digits = digits),
      "AIC:", formatC(AIC(x), digits = digits),
      "BIC:", formatC(betaregesc::BIC(x), digits = digits)
  )
}


#' Ajuste modelo beta com resposta ordinal transformada intervalar
#' 
#'
#' Função core que ajusta um modelo beta ordinal para resposta intervalar 
#' usando a função optim do pacote stats, retornando uma tabela com estatísticas
#' sumarizadas do ajuste, tanto para dispesão fixa como variável.
#'
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas e também Z1, Z2, ..., Zn para aquelas relacionadas com phi, se houver.
#' Ela deve ser referenciada em Y. Ex. formula = y ~ X1 + X2 ou formula = y ~ X1 + X2 | Z1
#' ou formula = y ~ X1 + X2 | Z1 + Z2, etc. Como a variável resposta é intervalar
#' ela deverá ser passada como left e right no objeto dados.
#' Veja os detalhes para mais informação.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula. Ele dever conter o limite inferior (left)
#'  e o limite superior (right) da variável resposta intervalar Y.
#' @param link Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_phi Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @param acumulada Se TRUE, retorna a verossimilhança pela pbeta, dbeta caso contrário.
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @importFrom Formula as.Formula
#' @examples
#' \dontrun{
#' n <- 100
#' dados <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
#'                     z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- betaregesc_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link = "logit",
#' link_phi = "log",
#' ncuts = 100)
#' fit1 <- betaregesc(y ~ x1 + x2,
#' dados = dados_simulados,
#' link = "logit",
#' link_phi = "log",
#' num_hessiana = TRUE)
#' print(fit1)
#' 
#' fit2 <- betaregesc(y ~ x1 + x2 | z1,
#' dados = dados_simulados,
#' link = "logit",
#' link_phi = "log",
#' num_hessiana = TRUE)
#' print(fit2)
#' }
#' @export
betaregesc <- function(formula, dados, link = "logit", link_phi = "identity", acumulada = TRUE,
                       ncuts = 100, type = "m", lim = 0.5, num_hessiana = TRUE){
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    #fx <- formula(terms(formula, rhs = 1L))
    out <- betaregesc_fit(formula = formula, dados = dados, link = link, link_phi = link_phi,
                          ncuts = ncuts, type = type, lim = lim, acumulada = acumulada,
                          num_hessiana = num_hessiana)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
    }
    #fx <- formula(terms(formula, rhs = 1L))
    #fz <- formula(terms(formula, rhs = 2L))
    out <- betaregesc_fit_z(formula = formula, dados = dados, link = link, link_phi = link_phi,
                            ncuts = ncuts, type = type, lim = lim, num_hessiana = num_hessiana, 
                            acumulada = acumulada)
  }
  class(out) <- c("betaregesc","betaregescdv")
  return(out)
}


#' Ajuste modelo beta com resposta ordinal transformada intervalar via bbmle
#' 
#'
#' Função core que ajusta um modelo beta ordinal para resposta intervalar 
#' usando a função optim do pacote stats, retornando uma tabela com estatísticas
#' sumarizadas do ajuste, tanto para dispesão fixa como variável.
#'
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas e também Z1, Z2, ..., Zn para aquelas relacionadas com phi, se houver.
#' Ela deve ser referenciada em Y. Ex. formula = y ~ X1 + X2 ou formula = y ~ X1 + X2 | Z1
#' ou formula = y ~X1 + X2 | Z1 + Z2, etc. Como a variável resposta é intervalar
#' ela deverá ser passada como left e right no objeto dados.
#' Veja os detalhes para mais informação.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula. Ele dever conter o limite inferior (left)
#'  e o limite superior (right) da variável resposta intervalar Y.
#' @param link Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_phi Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim.
#' 
#' @param optimizer Algoritmo de otimização. Pode ser "optim" (padrão) ou "nlminb". Veja
#' \code{\link{mle2}} para mais detalhes.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
#' @param lim Região de incerteza da medida. Padrão 0.5.
#' @param acumulada Se TRUE, retorna a verossimilhança pela pbeta, dbeta caso contrário.
#' 
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @importFrom Formula as.Formula
#' @examples
#' \dontrun{
#' require(bbmle)
#' n <- 100
#' dados <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
#'                     z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- betaregesc_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link = "logit",
#' link_phi = "log",
#' ncuts = 100)
#' fit <- betaregesc_bbmle(formula = y ~ x1 + x2|z1,
#' dados = dados_simulados, link = "logit", link_phi = "log")
#' p <- profile(fit)
#' plot(p)
#' }
#' @importFrom betareg betareg
#' @importFrom bbmle parnames
#' @importFrom stats as.formula delete.response model.response
#' @export
betaregesc_bbmle <- function(formula, dados, link = "logit", link_phi = "identity", 
                           num_hessiana = TRUE, acumulada = TRUE, optimizer = "optim",
                           ncuts = 100, type = "m", lim = 0.5
                           ){
  
  optimizer <- match.arg(optimizer, c("optim","nlminb"))
  formula <- Formula::as.Formula(formula)
  
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    fx <- formula(terms(formula, rhs = 1L))
    
    mf <- model.frame(fx, data = dados)
    mtX <- terms(formula, data = dados, rhs = 1L)
    Y <- fn_check_response(model.response(mf, "numeric"), ncuts = ncuts, type = type, lim = lim)
    
    fb <- betareg::betareg(formula = formula(paste0("yt ~ ", paste0(fx)[2])), 
                           data = data.frame(Y, X[,-1, drop = FALSE]), link = "logit")
    ini <- coef(fb)
    ll <- function(param, formula, dados, link, link_phi, acumulada){
      -betaregesc::betaregesc_log_vero(
        param = param, formula = formula, dados = dados,
        link = link, link_phi = link_phi, 
        acumulada = acumulada
      )
    }
    bbmle::parnames(ll) <- names(ini)
    out <- bbmle::mle2(minuslogl = ll, 
                       vecpar = TRUE,
                       optimizer = optimizer,
                       data = list(hessian=TRUE,
                                   dados = dados,
                                   formula = formula,
                                   link = link, 
                                   link_phi = link_phi, 
                                   ncuts = ncuts, type = type, lim = lim,
                                   acumulada = acumulada),
                       start = as.list(ini))
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
    }

    mf <- model.frame(formula, data = dados)
    mtX <- terms(formula, data = dados, rhs = 1L)
    mtZ <- delete.response(terms(formula, data = dados, rhs = 2L))
    Y <- fn_check_response(model.response(mf, "numeric"), ncuts = ncuts, type = type, lim = lim)
    y <- Y[,c("left","right")]
    X <- model.matrix(mtX, mf)
    Z <- model.matrix(mtZ, mf)

    fb <- betareg::betareg(formula = formula(paste0("yt ~ ", paste0(formula)[3])),
                           data = data.frame(Y, X[,-1, drop = FALSE], Z[,-1, drop = FALSE]),
                           link = "logit")
    ini <- coef(fb)
    ll2 <- function(param, dados, formula_x, formula_z, link_x, link_z, acumulada){
      -betaregesc::betaregesc_log_vero_z(
        param = param, formula = formula, 
        dados = dados, link = link, link_phi = link_phi, 
        acumulada = acumulada
      )
    }
    bbmle::parnames(ll2) <- names(ini)
    out <- bbmle::mle2(minuslogl = ll2, 
                       vecpar = TRUE,
                       optimizer = optimizer,
                       data = list(hessian=TRUE,
                                   dados = dados,
                                   formula_x = fx,
                                   formula_z = fz,
                                   link_x = link,
                                   link_z = link_phi,
                                   ncuts = ncuts, type = type, lim = lim,
                                   acumulada = acumulada),
                       start = as.list(ini))
  }
  return(out)
}



#' @title Verifica e ajusta os limites superior e inferior de uma variável resposta.
#'
#' @description Essa função verifica e ajusta os limites superior e inferior de uma variável resposta, com base no tipo de intervalo especificado ("m", "l" ou "r") e no número de quebras.
#'
#' @param y Variável resposta numérica.
#' @param type Caractere que indica o tipo de intervalo a ser utilizado. Deve ser um dos seguintes: "m" (meio), "l" (esquerda) ou "r" (direita). O padrão é "m".
#' @param ncuts Número inteiro que representa a quantidade de quebras a serem aplicadas. O padrão é 100L.
#' @param lim Limite numérico a ser utilizado para ajustar os intervalos. O padrão é 0.5.
#'
#' @return Retorna uma matriz com os limites ajustados para a variável resposta.
#' @examples
#' y <- c(0, 3, 5, 7, 9, 10)
#' result <- fn_check_response(y, type = "m", ncuts = 10)
#'
#' @export
fn_check_response <- function(y, type = "m", ncuts = 100L, lim = 0.5){
  if(all(y > 0 & y < 1)){
    message("Variavel recebida esta no intervalo unitario.")
    y <- y*ncuts
  }
  type  <- match.arg(type, c("m","l","r"))
  if(ncuts < max(y)){
    message("O numero de quebras maior que o valor maximo de da variavel resposta.\nO ideal que fosse no maximo igual.")
  }
  
  if(type == "m"){
    y_inf <- (y - lim)/ncuts
    y_sup <- (y + lim)/ncuts
  } else if(type == "l"){
    y_inf <- (y - lim*2)/ncuts
    y_sup <- (y)/ncuts
  } else if(type == "r"){
    y_inf <- (y)/ncuts
    y_sup <- (y + lim*2)/ncuts
  }
  y_inf[y_inf <= 0] <- 0.00001
  y_sup[y_sup <= 0] <- 0.00001
  y_inf[y_inf >= 1] <- 0.99999
  y_sup[y_sup >= 1] <- 0.99999
  
  Y <- cbind(left = y_inf, right = y_sup, yt=y/ncuts, y)
  return(Y)
}

