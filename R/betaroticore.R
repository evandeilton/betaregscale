## -------------------------------------------------------------------------- ##
## titulo:  REGRESSÃO BETA PARA DADOS COM RESPOSTA ORDINAL TRANSFORMADA
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
#' @return Retorna um data.frame contendo os dados simulados da variável
#'  beta ordinal e as variáveis independentes.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(42)
#' dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#'
#' # Simular dados usando a função beta_ordinal_simula_dados com parâmetros personalizados
#' dados_simulados <- beta_ordinal_simula_dados(
#'  formula = ~x1 + x2, dados = dados,
#'  betas = c(1, -0.3, 0.4), phi = 30,
#'  link = "probit", ncuts = 100, type = "m")
#' head(dados_simulados)
#' @export
beta_ordinal_simula_dados <- function(formula, dados, betas, phi = 50, link = "logit", ncuts = 100, type = "m"){
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

  if(type == "m"){
    y_inf <- (y_meio - 0.5)/ncuts
    y_sup <- (y_meio + 0.5)/ncuts
  } else if(type == "l"){
    y_inf <- (y_meio - 1.0)/ncuts
    y_sup <- (y_meio)/ncuts
  } else if(type == "r"){
    y_inf <- (y_meio)/ncuts
    y_sup <- (y_meio + 1.0)/ncuts
  }
  y_inf[y_inf <= 0] <- 0.00001
  y_sup[y_sup <= 0] <- 0.00001
  y_inf[y_inf >= 1] <- 0.99999
  y_sup[y_sup >= 1] <- 0.99999

  Y <- cbind(left = y_inf, right = y_sup)
  out <- data.frame(y, yr = y_meio, Y, X[,-1, drop = FALSE])
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
#' @param link_x Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity".
#'  O padrão é "logit".
#' @param link_z Nome da função de ligação a ser usada para as preditoras
#' Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param dados Um conjunto de dados que contém a variável dependente e as
#'  variáveis independentes especificadas nas fórmulas.
#' @param betas Vetor de betas associados aos preditores de X utilizados para simular da beta.
#' @param zetas Vetor de zetas associados aos preditores de Z utilizados para simular da beta com dispersão variável.
#' @param ncuts Número de cortes para a variável ordinal. O padrão é 100.
#' @param type Tipo de intervalo. "m" = meio; "l" = esquerda e "r" = direita.
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
#' dados_simulados <- beta_ordinal_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link_x = "logit",
#' link_z = "log",
#' ncuts = 100,
#' type = "m")
#' @export
beta_ordinal_simula_dados_z <- function(formula_x = ~x1 + x2,
                                        formula_z = ~z1 + z2,
                                        dados,
                                        betas = c(0, 0.5, -0.2),
                                        zetas = c(1, 0.5, 0.2),
                                        link_x = "logit",
                                        link_z = "log",
                                        ncuts = 100,
                                        type = "m"){
  mfx <- model.frame(formula_x, data = dados)
  mfz <- model.frame(formula_z, data = dados)
  X   <- model.matrix(mfx, data = dados)
  Z   <- model.matrix(mfz, data = dados)
  n   <- nrow(X)

  link_x <- match.arg(link_x, c("logit","probit","cauchit","cloglog","identity"))
  link_z <- match.arg(link_z, c("log","sqrt","identity"))

  # Aplicação do preditor linear nos betas
  mu_x <- fn_switch_link(eta = X%*%betas, link = link_x)
  mu_p <- fn_switch_link(eta = Z%*%zetas, link = link_z)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu_x * mu_p)
  beta <- as.numeric((1 - mu_x)*mu_p)
  y <- rbeta(n = n, shape1 = mu_x*mu_p, shape2 = (1-mu_x)*mu_p)
  y_meio <- round(y*ncuts, 0)

  if(type == "m"){
    y_inf <- (y_meio - 0.5)/ncuts
    y_sup <- (y_meio + 0.5)/ncuts
  } else if(type == "l"){
    y_inf <- (y_meio - 1.0)/ncuts
    y_sup <- (y_meio)/ncuts
  } else if(type == "r"){
    y_inf <- (y_meio)/ncuts
    y_sup <- (y_meio + 1.0)/ncuts
  }
  y_inf[y_inf <= 0] <- 0.00001
  y_sup[y_sup <= 0] <- 0.00001
  y_inf[y_inf >= 1] <- 0.99999
  y_sup[y_sup >= 1] <- 0.99999

  Y <- cbind(left = y_inf, right = y_sup)
  out <- data.frame(y, yr = y_meio, Y, X[,-1, drop = FALSE], Z[,-1, drop = FALSE])
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
#' @return Retorna a soma da log-verossimilhança dos dados.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(421)
#' dados <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' # Calcular a log-verossimilhança usando a função log_vero_beta_ordinal
#' param <- c(0, 0.5, -0.2, 50)
#' phi <- 30
#' formula <- ~x1 + x2
#' dados_simulados <- beta_ordinal_simula_dados(formula = formula, dados = dados,
#' betas = c(0, 0.5, -0.2), phi = phi, link = "logit", ncuts = 10, type = "m")
#' log_verossimilhanca <- beta_ordinal_log_vero(param, formula, dados_simulados)
#' print(log_verossimilhanca)
#' @export
beta_ordinal_log_vero <- function(param, formula, dados, link = "logit", link_phi = "log", acumulada = TRUE){
  link <- match.arg(link, c("logit","probit","cauchit","cloglog","identity"))
  link_phi <- match.arg(link_phi, c("log","sqrt","identity"))
  
  mfx <- model.frame(formula, data = dados)
  X   <- model.matrix(mfx, data = dados)
  betas  <- param[1:ncol(X)]

  # Aplicação do preditor linear nos betas
  eta <- X%*%betas
  mu <- fn_switch_link(eta = eta, link = link)
  phi <- fn_switch_link(eta = param[(length(betas)+1)], link = link_phi)

  # Prepara a log-verossimilhança com dados intervalares
  alpha <- as.numeric(mu * phi)
  beta <- as.numeric((1 - mu)*phi)

  ll <- if(acumulada){
    p1 <- pbeta(q = as.numeric(dados[,"left"]), shape1 = alpha, shape2 = beta)
    p2 <- pbeta(q = as.numeric(dados[,"right"]), shape1 = alpha, shape2 = beta)
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
#' @param formula_x Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas.
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param formula_z Fórmula para expressar a relação das preditoras Z1, Z2, Zn
#' relacionadas com os phi's no modelo se regressão beta com preditoas em phi.
#' Ela deve ser referenciada em Y. Ex. formula = ~Z1 + Z2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param link_x Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity".
#' O padrão é "logit".
#' @param link_z Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn
#'  relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#' independentes especificadas na fórmula.
#' @param acumulada Um valor lógico indicando se a log-verossimilhança acumulada
#' deve ser calculada. O padrão é TRUE.
#' @return Retorna a soma da log-verossimilhança dos dados.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' n <- 50
#' dados <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
#'                     z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- beta_ordinal_simula_dados_z(
#'   formula_x = fx,
#'   formula_z = fz,
#'   dados = dados,
#'   betas = c(0.2, -0.5, 0.3),
#'   zetas = c(1, 1.2),
#'   link_x = "logit",
#'   link_z = "log",
#'   ncuts = 100)
#' # Calcular a log-verossimilhança usando a função log_vero_beta_ordinal
#' log_verossimilhanca <- beta_ordinal_log_vero_z(
#'   param = c( c(0.2, -0.5, 0.3), c(1, 1.2)),
#'   formula_x = fx,
#'   formula_z = fz,
#'   dados = dados_simulados,
#'   link_x = "logit",
#'   link_z = "log",
#'   acumulada = TRUE)
#'   print(log_verossimilhanca)
#' @export
beta_ordinal_log_vero_z <- function(param,
                                    formula_x = ~x1 + x2,
                                    formula_z = ~z1 + z2,
                                    dados,
                                    link_x = "logit",
                                    link_z = "log",
                                    acumulada = TRUE){

  mfx <- model.frame(formula_x, data = dados)
  mfz <- model.frame(formula_z, data = dados)
  X   <- model.matrix(mfx, data = dados)
  Z   <- model.matrix(mfz, data = dados)
  n   <- nrow(X)

  link_x <- match.arg(link_x, c("logit","probit","cauchit","cloglog","identity"))
  link_z <- match.arg(link_z, c("log","sqrt","identity"))

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
    p1 <- pbeta(q = as.numeric(dados[,"left"]), shape1 = alpha, shape2 = beta)
    p2 <- pbeta(q = as.numeric(dados[,"right"]), shape1 = alpha, shape2 = beta)
    area <- p2 - p1
    area <- area + 0.00001
    log(area)
  } else {
    suppressWarnings(dbeta(dados[,"y"], shape1 = alpha, shape2 = beta, log = TRUE))
  }

  ll[is.infinite(ll)] <- NaN
  return(sum(ll, na.rm = TRUE))
}

#' Função para ajustar um modelo beta ordinal
#'
#' A função fit_beta_ordinal ajusta um modelo beta ordinal usando a função optim
#' do pacote stats, retornando uma tabela com estatísticas sumarizadas do ajuste,
#' incluindo intervalos de confiança e cálculo do BIAS.
#'
#' @param formula Fórmula para expressar a relação das preditoras X1, X2, Xn.
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula.
#' @param link Nome da função de ligação a ser usada. Pode ser uma das
#' seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_phi Função de ligação para phi. Uma dentre "log","sqrt","identity"
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @examples
#' # Criar um conjunto de dados de exemplo
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' coe <- beta_ordinal_coef(fit)
#' coe$est
#' @importFrom betareg betareg
#' @importFrom numDeriv hessian
#' @importFrom stats AIC cor fitted hatvalues logLik qlogis runif terms
#' @export
beta_ordinal_fit <- function(formula, dados, link = "logit", link_phi = "identity", num_hessiana = TRUE) {
  formula_ini <- paste0("y ", paste0(as.character(formula), collapse = " "))
  ini <- coef(betareg::betareg(formula = formula_ini, data = dados, 
                               link = link, link.phi = link_phi))

  # Ajustando o modelo com a função optim
  opt_result <- optim(par = ini, 
                      fn = beta_ordinal_log_vero, formula = formula, 
                      dados = dados, link = link, link_phi = link_phi,
                      hessian = !num_hessiana,
                      method = "BFGS", acumulada = TRUE,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(beta_ordinal_log_vero, opt_result$par, 
                                  formula = formula, dados = dados,
                                  link = link, link_phi = link_phi)
    opt_result$hessian <- hessiana
  }

  # Mu chapéu
  est <- opt_result$par
  y  <- as.matrix(dados[,c("left","right")])
  X <- model.matrix(formula, data = dados)
  hatmu <- fn_switch_link(eta = X%*%est[1:ncol(X)], link = link)
  hatphi <- est[ncol(X)+1]
  pseudor2 <- cor(X%*%est[1:ncol(X)], make.link(link)$linkfun(apply(y, 1, mean)))^2
  
  opt_result$dados <- dados
  opt_result$link <- link
  opt_result$link_phi <- link_phi
  opt_result$formula_x <- formula
  opt_result$formula_z <- ~ 1
  
  opt_result$residuals <- apply(y, 1, mean) - hatmu
  opt_result$hatmu  <- hatmu
  opt_result$hatphi <- hatphi
  opt_result$pseudo.r.squared = pseudor2
  
  class(opt_result) <- c("betaroti","betarotidv")
  
  return(invisible(opt_result))
}


#' Função para ajustar um modelo beta ordinal
#'
#' A função fit_beta_ordinal ajusta um modelo beta ordinal usando a função optim
#' do pacote stats, retornando uma tabela com estatísticas sumarizadas do ajuste,
#' incluindo intervalos de confiança e cálculo do BIAS.
#'
#' @param formula_x Fórmula para expressar a relação das preditoras X1, X2, Xn
#' relacionadas com os betas.
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param formula_z Fórmula para expressar a relação das preditoras Z1, Z2, Zn
#' relacionadas com os phi's no modelo se regressão beta com preditoas em phi.
#' Ela deve ser referenciada em Y. Ex. formula = ~Z1 + Z2. Isso porque a variável
#' resposta é intervalar. Veja os detalhes para mais informação.
#' @param dados Um conjunto de dados que contém a variável dependente e as variáveis
#'  independentes especificadas na fórmula.
#' @param link_x Nome da função de ligação a ser usada para as preditoras X1, X2, ..., Xn.
#' Pode ser uma das seguintes: "logit", "probit", "cauchit", "cloglog" ou "identity". O padrão é "logit".
#' @param link_z Nome da função de ligação a ser usada para as preditoras Z1, Z2, ..., Zn relacionadas com phi.
#' Pode ser uma das seguintes: "log", "sqrt e "identity". O padrão é "log".
#' @param num_hessiana Se TRUE, calcula a matriz Hessian numericamente com o
#' pacote numDeriv. Se FALSE, calcula com o padrão da optim
#' @return Retorna uma lista contendo o resultado da otimização e uma tabela com
#' estatísticas sumarizadas do ajuste.
#' @importFrom betareg betareg
#' @importFrom numDeriv hessian
#' @examples
#' n <- 50
#' dados <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
#'                     z1 = runif(n), z2 = runif(n))
#' fx <- ~ x1 + x2
#' fz <- ~ z1
#' dados_simulados <- beta_ordinal_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link_x = "logit",
#' link_z = "log",
#' ncuts = 100)
#' fit_z <- beta_ordinal_fit_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados_simulados,
#' link_x = "logit",
#' link_z = "log",
#' num_hessiana = TRUE)
#' beta_ordinal_coef(fit_z)$est
#' @export
beta_ordinal_fit_z <- function(formula_x,
                               formula_z,
                               dados,
                               link_x = "logit",
                               link_z = "log",
                               num_hessiana = TRUE) {

  formula_ini_composta <- paste0("y ", paste0(as.character(formula_x), collapse = " "), "|", as.character(formula_z)[2])
  ini <- c(coef(betareg::betareg(formula = formula_ini_composta, data = dados)))

  # Ajustando o modelo com a função optim
  opt_result <- optim(par = ini,
                      fn = beta_ordinal_log_vero_z,
                      formula_x = formula_x,
                      formula_z = formula_z,
                      link_x = link_x,
                      link_z = link_z,
                      dados = dados,
                      hessian = !num_hessiana,
                      method = "BFGS",
                      acumulada = TRUE,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(beta_ordinal_log_vero_z,
                                  opt_result$par,
                                  formula_x = formula_x, formula_z = formula_z, dados = dados)
    opt_result$hessian <- hessiana
  }
  # Mu chapéu
  est <- opt_result$par
  y  <- as.matrix(dados[,c("left","right")])
  X <- model.matrix(formula_x, data = dados)
  Z <- model.matrix(formula_z, data = dados)
  p <- 1:ncol(X)
  q <- (ncol(X)+1):length(est)
  hatmu  <- fn_switch_link(eta = X%*%est[p], link = link_x)
  hatphi <- fn_switch_link(eta = Z%*%est[q], link = link_z)
  pseudor2 <- cor(X%*%est[p], make.link(link_x)$linkfun(apply(y, 1, mean)))^2

  opt_result$dados <- dados
  opt_result$link <- link_x
  opt_result$link_phi <- link_z
  opt_result$formula_x <- formula_x
  opt_result$formula_z <- formula_z
  opt_result$residuals <- apply(y, 1, mean) - hatmu
  opt_result$hatmu  <- hatmu
  opt_result$hatphi <- hatphi
  opt_result$pseudo.r.squared = pseudor2
  
  class(opt_result) <- c("betaroti","betarotidv")
  
  return(invisible(opt_result))
}


#' Coleta estatística do ajuste
#'
#' Coleta as estimativas e suas estatísticas para um objeto de ajuste
#' do modelo beta ordinal com censura intervalar. Coleta também as estatísticas
#' de bondade do ajuste como a log-verossimilhança, o AIC e o BIC do modelo
#'
#' @param fit Objeto do sjuste retornado das funções \link{beta_ordinal_fit} e
#' \link{beta_ordinal_fit_z}
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
#'  dados_simulados <- beta_ordinal_simula_dados_z(
#'    formula_x = fx,
#'    formula_z = fz,
#'    dados = dados,
#'    betas = c(0.2,-0.5, 0.3),
#'    zetas = c(1, 1.2),
#'    link_x = "logit",
#'    link_z = "log",
#'    ncuts = 100
#'    )
#'
#'  fit_z <- beta_ordinal_fit_z(
#'   formula_x = fx,
#'   formula_z = fz,
#'   dados = dados_simulados,
#'   link_x = "logit",
#'   link_z = "log",
#'   num_hessiana = TRUE
#'   )
#'
#' coe <- beta_ordinal_coef(fit_z)
#' coe$est
#' coe$gof
#' @return Lista contendo as estimativas coeficientes, suas estatísticas e a bondade do ajuste.
#' @export
beta_ordinal_coef <- function(fit, alpha = 0.05){
  
  if(!inherits(fit, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(fit), collapse = ","), "' nao suportada.\n"))
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
#' Esta função calcula a log-verossimilhança negativa para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param ... Argumentos extras
#' 
#' @examples
#' \dontrun{
#' # Exemplo de uso da função logLik
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' # Em seguida, use a função logLik.betaroti
#' log_likelihood <- logLik(fit)
#'}
#' @method logLik betaroti
#' @export
logLik.betaroti <- function(object, ...){
  
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  -object$value
}


#' @title Critério de Informação de Akaike (AIC)
#'
#' @description
#' Esta função calcula o Critério de Informação de Akaike (AIC) para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param ... Argumentos extras
#' @param k Penalidade por parâmetro a ser utilizado; o padrão k = 2 é o AIC clássico.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função AIC
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method AIC betaroti
#' @export
AIC.betaroti <- function(object, ..., k = 2){
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  k * (length(object$par) + 1) - 2 * (-object$value)
}


#' @title Critério de Informação Bayesiano (BIC)
#'
#' @description
#' Esta função calcula o Critério de Informação Bayesiano (BIC) para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função BIC
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  log(nrow(object$dados)) * (length(object$par) + 1) - 2 * (-object$value)
}

#' @title Matriz hessiana
#'
#' @description
#' Esta função retorna a matriz hessiana para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  object$hessian
}

#' @title Coeficientes
#'
#' @description
#' Esta função retorna os coeficientes para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param ... outros argumentos.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method coef betaroti
#' @export
coef.betaroti <- function(object, ...){
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  object$par
}

#' @title Medidas de ajuste
#'
#' @description
#' Esta função retorna um conjunto de medidas de ajuste, incluindo log-verossimilhança, AIC e BIC, para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
    BIC = betaroti::BIC(object)
  )
}

#' @title Estimativas e intervalos de confiança
#'
#' @description
#' Esta função retorna estimativas, erros padrão, intervalos de confiança e valores-p para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param alpha Nível de significância para os intervalos de confiança (padrão é 0,05).
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
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
#' Esta função retorna um resumo das estimativas, erros padrão, intervalos de confiança e valores-p para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param alpha Nível de significância para os intervalos de confiança (padrão é 0,05).
#' @param ... Outros argumentos
#'
#' @return
#' Retorna um resumo das estimativas, erros padrão, intervalos de confiança e valores-p do objeto fornecido.
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' # Em seguida, use a função summary.betaroti
#' summary_result <- summary.betaroti(fit, alpha = 0.05)
#'}
#' @method summary betaroti
#' @export
summary.betaroti <- function(object, ..., alpha = 0.05){
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  beta_ordinal_coef(object, alpha = alpha)
}

#' @title Resíduos do modelo
#'
#' @description
#' Esta função retorna os resíduos para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param type Tipo de resíduo. Um entre "rqa","deviance", "pearson", "response", "weighted", "sweighted". O padrão é Resíduo Quantílico Aleatorizado (rqa)
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method residuals betaroti
#' @export
residuals.betaroti <- function(object, type = c("rqa","deviance", "pearson", "response", "weighted", "sweighted"), ...){
  type <- match.arg(type)
  
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
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
#' Esta função retorna os valores ajustados (mu ou phi) para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param type Tipo de valor ajustado a ser retornado: "mu" ou "phi" (padrão é "mu").
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method fitted betaroti
#' @export
fitted.betaroti <- function(object, type = "mu", ...){
  type <- match.arg(type, c("mu","phi"))
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
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
#' Esta função retorna as previsões de quantis para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param type Tipo de previsão a ser retornado: "quantile" (padrão é "quantile").
#' @param at Valor numérico entre 0 e 1 para o qual o quantil é desejado (padrão é 0,5).
#' @param ... Outros argumentos
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method predict betaroti
#' @export
predict.betaroti <- function(object, type = "quantile", at = 0.5, ...){
  type <- match.arg(type, c("quantile"))
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
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
#' Esta função retorna a matriz de covariância para objetos das classes 'betaroti' e 'betarotidv'.
#'
#' @param object Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param ... Outros parametros
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method vcov betaroti
#' @export
vcov.betaroti <- function(object, ...){
  
  if(!inherits(object, c("betaroti","betarotidv"))){
    stop(paste0("log: Preciso de um objeto da classe 'betaroti' ou 'betarotidv'. Classe '", paste0(class(object), collapse = ","), "' nao suportada.\n"))
  }
  return(solve(-as.matrix(object$hessian)))
}

#' @title Print
#'
#' @description
#' Esta função imprime um resumo do objeto das classes 'betaroti' e 'betarotidv' fornecido, incluindo estimativas, erros padrão, valores t e valores-p.
#'
#' @param x Um objeto das classes 'betaroti' ou 'betarotidv'.
#' @param digits Número de dígitos significativos para a impressão dos valores (padrão é max(3, getOption("digits") - 2)).
#' @param ... outros parametros
#' @examples
#' \dontrun{
#' # Exemplo de uso da função hessian
#' # Primeiro, gere um objeto de classe 'betaroti' ou 'betarotidv'
#' set.seed(42)
#' n <- 100
#' dados <- data.frame(
#'   x1 = rnorm(n, mean = 1, sd = 0.5),
#'   x2 = rbinom(n, size = 1, prob = 0.5),
#'   x3 = rnorm(n, mean = 2, sd = 1))
#' betas <- c(0.2, 0.3, -0.4, 0.1)
#' formula <- ~x1 + x2 + x3
#' phi <- 50
#' dados_simulados <- beta_ordinal_simula_dados(
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
#' @method print betaroti
#' @export
print.betaroti <- function(x, digits = max(3, getOption("digits") - 2), ...){
  cf <- beta_ordinal_coef(x)
  be <- cf$est[,c("estimate","se","t_value","p_value")]
  colnames(be) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(be) <- cf$est$variable
  stats::printCoefmat(be, digits = digits, P.values = TRUE)
  cat("\nLog-Likelihood:", formatC(logLik(x), digits = digits),
      "AIC:", formatC(AIC(x), digits = digits),
      "BIC:", formatC(betaroti::BIC(x), digits = digits)
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
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2 ou formula = ~X1 + X2 | Z1
#' ou formula = ~X1 + X2 | Z1 + Z2, etc. Como a variável resposta é intervalar
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
#' @details
#' O objeto dados precisa ter duas colunas nomeadas contendo os dados do 
#' limite inferior (\code{left}) e superior (\code{right}) da variável resposta.
#' Isto é, além de conter as preditoras que serão inserida no modelo via formula,
#' os dados devem conter as duas colunas de y em forma de intervalo.
#' As formulas não devem conter mensão a y. Como descrito no parâmtro formula, 
#' opções válidas formula = ~X1 + X2 ou formula = ~X1 + X2 | Z1 ou 
#' formula = ~X1 + X2 | Z1 + Z2 sem parte esquerda. Essa formula vai definir
#' apenas as preditoras que entram no modelo sejam para mu ou para phi.
#' 
#' 
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
#' dados_simulados <- beta_ordinal_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link_x = "logit",
#' link_z = "log",
#' ncuts = 100)
#' fit1 <- betaroti( ~ x1 + x2,
#' dados = dados_simulados,
#' link = "logit",
#' link_phi = "log",
#' num_hessiana = TRUE)
#' print(fit1)
#' 
#' fit2 <- betaroti( ~ x1 + x2 | z1,
#' dados = dados_simulados,
#' link = "logit",
#' link_phi = "log",
#' num_hessiana = TRUE)
#' print(fit2)
#' }
#' @export
betaroti <- function(formula, dados, link = "logit", link_phi = "identity", num_hessiana = TRUE){
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    fx <- formula(terms(formula, rhs = 1L))
    out <- beta_ordinal_fit(formula = fx, dados = dados, link = link, 
                            link_phi = link_phi, num_hessiana = num_hessiana)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
    }
    fx <- formula(terms(formula, rhs = 1L))
    fz <- formula(terms(formula, rhs = 2L))
    out <- beta_ordinal_fit_z(formula_x = fx, formula_z = fz, dados = dados, 
                              link_x = link, link_z = link_phi, num_hessiana = num_hessiana)
  }
  class(out) <- c("betaroti","betarotidv")
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
#' Ela deve ser referenciada em Y. Ex. formula = ~X1 + X2 ou formula = ~X1 + X2 | Z1
#' ou formula = ~X1 + X2 | Z1 + Z2, etc. Como a variável resposta é intervalar
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
#' 
#' @param acumulada Se acumulada ou distribuição.
#' 
#' @details
#' O objeto dados precisa ter duas colunas nomeadas contendo os dados do 
#' limite inferior (\code{left}) e superior (\code{right}) da variável resposta.
#' Isto é, além de conter as preditoras que serão inserida no modelo via formula,
#' os dados devem conter as duas colunas de y em forma de intervalo.
#' As formulas não devem conter mensão a y. Como descrito no parâmtro formula, 
#' opções válidas formula = ~X1 + X2 ou formula = ~X1 + X2 | Z1 ou 
#' formula = ~X1 + X2 | Z1 + Z2 sem parte esquerda. Essa formula vai definir
#' apenas as preditoras que entram no modelo sejam para mu ou para phi.
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
#' dados_simulados <- beta_ordinal_simula_dados_z(
#' formula_x = fx,
#' formula_z = fz,
#' dados = dados,
#' betas = c(0.2, -0.5, 0.3),
#' zetas = c(1, 1.2),
#' link_x = "logit",
#' link_z = "log",
#' ncuts = 100)
#' fit <- betaroti_bbmle(formula =  ~ x1 + x2|z1,
#' dados = dados_simulados, link = "logit", link_phi = "log")
#' p <- profile(fit)
#' plot(p)
#' }
#' @importFrom betareg betareg
#' @importFrom bbmle parnames
#' @export
betaroti_bbmle <- function(formula, dados, link = "logit", link_phi = "identity", 
                           num_hessiana = TRUE, acumulada = TRUE, optimizer = "optim"){
  
  optimizer <- match.arg(optimizer, c("optim","nlminb"))
  
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    fx <- formula(terms(formula, rhs = 1L))
    fb <- betareg::betareg(formula = formula(paste0("y ~ ", paste0(fx)[2])),
                           data = dados, link = "logit")
    ini <- coef(fb)
    ll <- function(param, formula, dados, link, link_phi, acumulada){
      -betaroti::beta_ordinal_log_vero(
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
                                   formula = fx,
                                   link = link, 
                                   link_phi = link_phi, 
                                   acumulada = acumulada),
                       start = as.list(ini))
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
    }
    fx <- formula(terms(formula, rhs = 1L))
    fz <- formula(terms(formula, rhs = 2L))
    fb <- betareg::betareg(formula = formula(paste0("y ~ ", paste0(formula)[2])),
                           data = dados, link = "logit")
    ini <- coef(fb)
    ll2 <- function(param, dados, formula_x, formula_z, link_x, link_z, acumulada){
      -betaroti::beta_ordinal_log_vero_z(
        param = param, formula_x = formula_x, formula_z = formula_z,
        dados = dados, link_x = link_x, link_z = link_z, 
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
                                   acumulada = acumulada),
                       start = as.list(ini))
  }
  return(out)
}
