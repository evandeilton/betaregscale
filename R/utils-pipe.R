#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats coef dbeta make.link model.frame model.matrix optim pbeta pt qnorm rbeta
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Classe betaregscale
#' 
#' @name cls betaregscale
#' @keywords internal
#' @importFrom methods setOldClass
methods::setOldClass(Classes = c("betaregscale","betaregscaledv"))  

#' Variaveis globais
#' @name Variaveis
#' @keywords internal
#@export
utils::globalVariables(c(".","::",":::","ll","fit","est","fz"))


