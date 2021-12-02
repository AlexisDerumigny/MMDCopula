
#' @title Convert between parameter and Kendall's tau for Marshall-Olkin copulas
#'
#' @param par the parameter of the Marshall-Olkin copula
#'
#' @return Either the Kendall's tau or the parameter of the Marshall-Olkin copula.
#'
#' @examples
#' BiCopPar2Tau.MO(par = 0.5)
#' BiCopTau2Par.MO(tau = 1/3)
#'
#' @references
#' Nelsen, R. B. (2007). An introduction to copulas. Springer Science & Business Media.
#' (Example 5.5)
#'
#' @export
#'
BiCopPar2Tau.MO <- function(par)
{
  return (par / (2 - par))
}

#' @param tau the Kendall's tau of the Marshall-Olkin copula
#'
#' @rdname BiCopPar2Tau.MO
#' @export
BiCopTau2Par.MO <- function(tau)
{
  return ( 2*tau / (1+tau) )
}

