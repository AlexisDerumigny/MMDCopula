
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
#' Can be \code{"cdf"} for the distance between the two cumulative distribution functions
#' or \code{"pdf"} for the distance between the two probability density functions.
#'
#' @param maxEval maximum number of evaluation of the function to be integrated.
#' If 0, then no maximum limit is given. (Only used if \code{p < Inf}).
#'
#' @param truncVal the distance is computed using the supremum or the integral
#' of the function on \eqn{[truncVal, 1 - truncVal]^2}.
#'
#' @return If \code{p < Inf}, it returns a list of four items
#'   \itemize{
#'     \item \code{distance} the value of the distance
#'     \item \code{integral} the value of the integral, which is
#'      the \eqn{p}-th power of the distance.
#'     \item \code{error} the estimated relative error of the integral
#'     \item \code{returnCode} the integer return code of the C routine
#'      called by \code{cubature::\link[cubature]{hcubature}()}.
#'      This should be 0 if there is no error.
#'   }
#' If \code{p = Inf}, it returns a list of two items
#'   \itemize{
#'     \item \code{distance} the maximum difference between the two copulas
#'     (respectively, between the two copula densities).
#'     \item \code{u_max} the point at which this difference is attained.
#'  }
#'
#' @examples
#' # Distance between the densities of a Gaussian copula with correlation 0.5
#' # and a Gaussian copula with correlation 0.2
#' BiCopParamDistLp(family = 1, par = 0.5, par_p = 0.2, p = 2, type = "cdf", maxEval = 10)
#' BiCopParamDistLp(family = 1, par = 0.5, par_p = 0.2, p = Inf, type = "cdf")
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
                             p, type, maxEval = 0, truncVal = 0)
{
  if (p < Inf) {
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
                                 lower = c(truncVal, truncVal),
                                 upper = c(1 - truncVal, 1 - truncVal),
                                 vectorInterface = TRUE, maxEval = maxEval)
    result[["distance"]] = result$integral^(1/p)

    return (result[c("distance", "integral", "error", "returnCode")])

  } else if (p == Inf){
    # if (maxEval == 0){maxEval <- 10000}
    # nGrid = floor(sqrt(maxEval))
    # univGrid = seq(0 , nGrid/(nGrid+1), length = nGrid) + 1/(2 * nGrid)
    # u = expand.grid(univGrid , univGrid)
    # switch(
    #   type,
    #   "cdf" = {
    #     diff_ =
    #       VineCopula::BiCopCDF(u1 = u[,1], u2 = u[,2],
    #                              family = family, par = par, par2 = par2) -
    #           VineCopula::BiCopCDF(u1 = u[,1], u2 = u[,2],
    #                                family = family_p, par = par_p, par2 = par2_p)
    #   },
    #
    #   "pdf" = {
    #     diff_ =
    #       VineCopula::BiCopPDF(u1 = u[,1], u2 = u[,2],
    #                              family = family, par = par, par2 = par2) -
    #           VineCopula::BiCopPDF(u1 = u[,1], u2 = u[,2],
    #                                family = family_p, par = par_p, par2 = par2_p)
    #   })
    #
    # result = max(abs(diff_))
    # which_ = which.max(abs(diff_))
    # u_max = as.numeric(u[which_, ])


    switch(
      type,
      "cdf" = {
        fn <- function(u) {
          return ( - (VineCopula::BiCopCDF(u1 = u[1], u2 = u[2],
                                           family = family, par = par, par2 = par2) -
                        VineCopula::BiCopCDF(u1 = u[1], u2 = u[2],
                                             family = family_p, par = par_p, par2 = par2_p)
          )^2 )
        }
        gr <- function(u) {
          return (
            - 2 *
              (unlist(VineCopula::BiCopHfunc(u1 = u[1], u2 = u[2],
                                             family = family, par = par, par2 = par2)) -
                 unlist(VineCopula::BiCopHfunc(u1 = u[1], u2 = u[2],
                                               family = family_p, par = par_p, par2 = par2_p))
              ) * (VineCopula::BiCopCDF(u1 = u[1], u2 = u[2],
                                        family = family, par = par, par2 = par2) -
                     VineCopula::BiCopCDF(u1 = u[1], u2 = u[2],
                                          family = family_p, par = par_p, par2 = par2_p)
              ))
        }
      },

      "pdf" = {
        fn <- function(u) {
          return ( - (
            VineCopula::BiCopPDF(u1 = u[1], u2 = u[2],
                                 family = family, par = par, par2 = par2) -
              VineCopula::BiCopPDF(u1 = u[1], u2 = u[2],
                                   family = family_p, par = par_p, par2 = par2_p) )^2)
        }

        gr <- function(u) {
          return (
            - 2 *
              c(VineCopula::BiCopDeriv(u1 = u[1], u2 = u[2],
                                       family = family, par = par, par2 = par2,
                                       deriv = "u1") -
                  VineCopula::BiCopDeriv(u1 = u[1], u2 = u[2],
                                         family = family_p, par = par_p, par2 = par2_p,
                                         deriv = "u1")
                ,
                VineCopula::BiCopDeriv(u1 = u[1], u2 = u[2],
                                       family = family, par = par, par2 = par2,
                                       deriv = "u2") -
                  VineCopula::BiCopDeriv(u1 = u[1], u2 = u[2],
                                         family = family_p, par = par_p, par2 = par2_p,
                                         deriv = "u2")
              ) * (VineCopula::BiCopPDF(u1 = u[1], u2 = u[2],
                                        family = family, par = par, par2 = par2) -
                     VineCopula::BiCopPDF(u1 = u[1], u2 = u[2],
                                          family = family_p, par = par_p, par2 = par2_p)
              ))
        }
      })

    result = optim(par = c(0.2, 0.2), fn = fn, gr = gr,
                   lower = truncVal, upper = 1 - truncVal, method = "L-BFGS-B")

    return (list(distance = -result$value, u_max = result$par ))
  }
}


