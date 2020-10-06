
#' Compute the distance between 2 parametric copulas
#'
#' This function uses the numerical integration procedure
#' \code{cubature::\link[cubature]{hcubature}()} to numerical integrate the distance between
#' the distribution or between the densities of two bivariate copulas.
#'
#' @param family family of the first copula.
#'
#' @param par first parameter of the first copula.
#'
#' @param par_p first parameter of the second copula.
#'
#' @param par2 second parameter of the first copula
#' (only useful for two-parameter families of copulas).
#'
#' @param par2_p second parameter of the first copula
#' (only useful for two-parameter families of copulas).
#'
#' @param family_p family of the second copula.
#'
#' @param p determines the \eqn{L_p} distance that is used.
#'
#' @param type type of the functions considered.
#' Can be \code{cdf} for the distance between the two cumulative distribution functions
#' or \code{pdf} for the distance between the two probability density functions.
#'
#' @param maxEval maximum number of evaluation of the function be integrated.
#' If 0, then no maximum limit is given.
#'
#' @return a list of four items
#'   \itemize{
#'     \item \code{distance} the value of the distance
#'     \item \code{integral} the value of the integral, which is
#'      the \eqn{p}-th power of the distance.
#'     \item \code{error} the estimated relative error of the integral
#'     \item \code{returnCode} the integer return code of the C routine
#'      called by \code{cubature::\link[cubature]{hcubature}()}.
#'      This should be 0 if there is no error.
#'   }
#'
#' @examples
#' # Distance between the densities of a Gaussian copula with correlation 0.5
#' # and a Gaussian copula with correlation 0.2
#' BiCopParamDistLp(family = 1, par = 0.5, par_p = 0.2, p = 2, type = "cdf", maxEval = 10)
#'
#' # Distance between the cdf of a Student copula
#' # with correlation 0.5 and 4 degrees of freedom
#' # and a Student copula with the same correlation but 20 degrees of freedom
#' BiCopParamDistLp(family = 2, par = 0.5, par_p = 0.5,
#' par2 = 5, par2_p = 20, p = 2, type = "pdf", maxEval = 10)
#'
#' # Distance between the densities of a Gaussian copula with correlation 0.5
#' # and of a Student copula with correlation 0.5 and 15 degrees of freedom
#' BiCopParamDistLp(family = 1, par = 0.5, par_p = 0.5, par2_p = 15,
#' family_p = 2, p = 2, type = "pdf", maxEval = 10)
#'
#' @export
#'
BiCopParamDistLp <- function(family, par, par_p, par2 = par, par2_p = par_p, family_p = family,
                             p, type, maxEval = 0)
{
  switch(
    type,
    "cdf" = {
      toIntegrate <- function(u) {
        return ( matrix(( abs(
          VineCopula::BiCopCDF(u1 = u[1,], u2 = u[2,],
                               family = family, par = par, par2 = par2) -
            VineCopula::BiCopCDF(u1 = u[1,], u2 = u[2,],
                                 family = family_p, par = par_p, par2 = par2_p)
        ))^p, nrow = 1, ncol = length(u[1,])) )
      }
    },

    "pdf" = {
      toIntegrate <- function(u) {
        return ( matrix(( abs(
          VineCopula::BiCopPDF(u1 = u[1,], u2 = u[2,],
                               family = family, par = par, par2 = par2) -
            VineCopula::BiCopPDF(u1 = u[1,], u2 = u[2,],
                                 family = family_p, par = par_p, par2 = par2_p)
        ))^p, nrow = 1, ncol = length(u[1,])) )
      }
    })

  result = cubature::hcubature(f = toIntegrate,
                               lower = c(0,0), upper = c(1,1),
                               vectorInterface = TRUE, maxEval = maxEval)
  result[["distance"]] = result$integral^(1/p)

  return (result[c("distance", "integral", "error", "returnCode")])
}


