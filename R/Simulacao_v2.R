## -------------------------------------------------------------------------- ##
## titulo:  REGRESSÃO BETA PARA DADOS COM RESPOSTA ORDINAL TRANSFORMADA
##          INTERVALAR
## Autor:   José Lopes / Wagner Bonatl
## Code:    Simulação da beta
## -------------------------------------------------------------------------- ##

#' Escolha das funções de ligação inversas
#'
#' Chave para escolha das funções de ligação inversas para os preditores de X e Z
#' @param eta Preditor linear dado por X*param.
#' @param link Função de ligação preferida.
#' @return Transformação inversa em eta, ou seja, mu.
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
beta_ordinal_log_vero <- function(param, formula, dados, link = "logit", acumulada = TRUE){
  mfx <- model.frame(formula, data = dados)
  X   <- model.matrix(mfx, data = dados)
  betas  <- param[1:ncol(X)]

  # Aplicação do preditor linear nos betas
  eta <- X%*%betas
  mu <- fn_switch_link(eta = eta, link = link)
  phi <- param[(length(betas)+1)]

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
#' @export
beta_ordinal_fit <- function(formula, dados, link = "logit", num_hessiana = TRUE) {
  formula_ini <- paste0("y ", paste0(as.character(formula), collapse = " "))
  ini <- coef(betareg::betareg(formula = formula_ini, data = dados))

  # Ajustando o modelo com a função optim
  opt_result <- optim(par = ini, fn = beta_ordinal_log_vero, formula = formula, dados = dados,
                      hessian = !num_hessiana,
                      method = "BFGS", acumulada = TRUE,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(beta_ordinal_log_vero, opt_result$par, formula = formula, dados = dados)
    opt_result$hessian <- hessiana
  }

  # Mu chapéu
  est <- opt_result$par
  X <- model.matrix(formula, data = dados)
  hatmu <- fn_switch_link(eta = X%*%est[1:ncol(X)], link = link)
  hatphi <- est[ncol(X)+1]
  opt_result$dados <- cbind(dados, hatmu, hatphi)
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
                      method = "BFGS", acumulada = TRUE,
                      control=list(fnscale = -1))

  if(num_hessiana){
    hessiana <- numDeriv::hessian(beta_ordinal_log_vero_z,
                                  opt_result$par,
                                  formula_x = formula_x, formula_z = formula_z, dados = dados)
    opt_result$hessian <- hessiana
  }
  # Mu chapéu
  est <- opt_result$par
  X <- model.matrix(formula_x, data = dados)
  Z <- model.matrix(formula_z, data = dados)
  p <- 1:ncol(X)
  q <- (ncol(X)+1):length(est)
  hatmu  <- fn_switch_link(eta = X%*%est[p], link = link_x)
  hatphi <- fn_switch_link(eta = Z%*%est[q], link = link_z)
  opt_result$dados <- cbind(dados, hatmu, hatphi)
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
beta_ordinal_coef <- function(fit){
  dados <- fit$dados
  # Obtendo os parâmetros ajustados
  beta_hat <- fit$par
  # Criando uma tabela com as estatísticas sumarizadas do ajuste
  summary_table <- data.frame(
    #beta = beta_hat,
    log_likelihood = -fit$value,
    AIC = 2 * (length(beta_hat) + 1) - 2 * (-fit$value),
    BIC = log(nrow(dados)) * (length(beta_hat) + 1) - 2 * (-fit$value)
  )

  # Obtendo os intervalos de confiança para as estimativas dos coeficientes
  cov_mat <- -solve(fit$hessian)
  se_beta <- sqrt(diag(cov_mat))
  z_alpha <- qnorm(1 - 0.05/2)
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
