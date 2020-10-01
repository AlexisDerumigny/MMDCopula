
#' Computation of the gradient of the MMD criterion
#' for parametric bivariate copulas models
#'
#' This function computes a stochastic estimate of the gradient of the MMD criterion
#' for parametric estimation of bidimensional copula family.
#' The main arguments are the two vectors of observations, and the copula family.
#' The family is parametrized as in \code{VineCopula::\link[VineCopula]{BiCop}()},
#' using the Kendall's tau instead of the first parameter.
#' This function is used by \code{\link{BiCopEstMMD}()} to perform parameter estimation
#' via MMD minimization.
#'
#'
#' @param u1 vector of observations of the first coordinate, in \eqn{[0,1]}.
#' @param u2 vector of observations of the second coordinate, in \eqn{[0,1]}.
#'
#' @param tau the copula family can be parametrized by the parameter \code{par}
#'   or by Kendall's tau. This function assumes a Kendall tau parametrization.
#'   Thus, the user can choose the value of Kendall tau at
#'   which the stochastic gradient should be computed.
#'
#' @param par if different from \code{NULL},
#'   the user must instead of \code{tau} specify the corresponding parameter \code{par}.
#'   The value of \code{tau} is then ignored.
#'
#' @param par2 value for the second parameter, if any. (Works only for Student copula).
#'
#' @param kernel the kernel used in the MMD distance:
#'   it can be a function taking in parameter \code{(u1, u2, v1, v2, gamma, alpha)}
#'   or a name giving the kernel to use in the list:
#'   \itemize{
#'     \item \code{gaussian}: Gaussian kernel \eqn{k(x,y) = \exp(-\|\frac{x-y}{\gamma}\|_2^2)
#'     }{k(x,y) = exp( - || (x-y) / gamma ||_2^2)}
#'     \item \code{exp.l2}: \eqn{k(x,y) = \exp(-\|\frac{x-y}{\gamma}\|_2)
#'     }{k(x,y) = exp( - || (x-y) / gamma ||_2)}
#'     \item \code{exp.l1}: \eqn{k(x,y) = \exp(-\|\frac{x-y}{\gamma}\|_1)
#'     }{k(x,y) = exp( - || (x-y) / gamma ||_1)}
#'     \item \code{inv.l2}: \eqn{k(x,y) = 1/(1+\|\frac{x-y}{\gamma}\|_2)^\alpha
#'     }{k(x,y) = 1 / ( 1 + || (x-y) / gamma ||_2 )^\alpha}
#'     \item \code{inv.l1}: \eqn{k(x,y) = 1/(1+\|\frac{x-y}{\gamma}\|_1)^\alpha
#'     }{k(x,y) = 1 / ( 1 + || (x-y) / gamma ||_1 )^\alpha}
#'   }
#'  Each of these names can receive the suffix ".KG", such as "gaussian.KG"
#'  to indicates that the kernel \eqn{k(x,y)} is replaced by
#'  \eqn{k(\Phi^{-1}(x) , \Phi^{-1}(y))} where \eqn{\Phi^{-1}} denotes the quantile
#'  function of the standard Normal distribution.
#'
#' @param gamma parameter \eqn{\gamma} to be used in the kernel.
#'
#' @param alpha parameter \eqn{\alpha} to be used in the kernel, if any.
#'
#' @param epsilon the differential of \code{VineCopula::\link[VineCopula]{BiCopTau2Par}()}
#'   is computed thanks to a finite difference with increment \code{epsilon}.
#'
#' @param ndrawings number of replicas of the stochastic estimate of the gradient drawn
#' at each step. The gradient is computed using the average of these replicas.
#'
#' @param family the chosen family of copulas
#'   (see the documentation of the class \code{VineCopula::\link[VineCopula]{BiCop}()}
#'   for the available families).
#'
#' @param method the method of computing the stochastic gradient:
#'   \itemize{
#'     \item \code{MC}: classical Monte-Carlo with \code{ndrawings} replications.
#'     \item \code{QMCV}: usual Monte-Carlo on U with \code{ndrawings} replications,
#'       quasi Monte-Carlo on V.
#'   }
#'
#' @param quasiRNG a function giving the quasi-random points in \eqn{[0,1]^2} or a name giving
#'   the method to use in the list: \itemize{
#'    \item \code{sobol}: use of the Sobol sequence
#'      implemented in \code{randtoolbox::\link[randtoolbox:quasiRNG]{sobol}}
#'    \item \code{halton}: use of the Halton sequence
#'      implemented in \code{randtoolbox::\link[randtoolbox:quasiRNG]{halton}}
#'    \item \code{torus}: use of the Torus sequence
#'      implemented in \code{randtoolbox::\link[randtoolbox:quasiRNG]{torus}}
#'   }
#'
#'
#' @return the value of the gradient.
#'
#' @seealso \code{\link{BiCopEstMMD}()} for the estimation of parametric bivariate copulas by
#' stochastic gradient descent on the MMD criteria.
#'
#' @examples
#' # Simulation from a bivariate Gaussian copula with correlation 0.5.
#' dataSampled = VineCopula::BiCopSim(N = 500, family = 1, par = 0.5)
#'
#' # computation of the gradient of the MMD criteria at different points
#' # Gradient is small at the true parameter
#' BiCopGradMMD(dataSampled[,1], dataSampled[,2], family = 1, par = 0.5)
#' # Gradient is negative when below the parameter
#' BiCopGradMMD(dataSampled[,1], dataSampled[,2], family = 1, par = 0.1)
#' # and positive when above
#' BiCopGradMMD(dataSampled[,1], dataSampled[,2], family = 1, par = 0.8)
#'
#' @export
#'
BiCopGradMMD <- function(
  u1, u2, family, tau, par = NULL, par2=0,
  kernel = "gaussian.KG", gamma=0.95, alpha=1, epsilon=0.0001,
  method = "QMCV", quasiRNG = "sobol", ndrawings=10)
{
  if (!is.null(par))
  {
    tau = VineCopula::BiCopPar2Tau(family = family, par = par, par2 = par2)
  }

  # Choice of the kernel
  if (is.character(kernel))
  {
    kernelFun <- findKernelFunction(kernel)
  } else {
    kernelFun <- kernel
  }

  switch (
    method,
    "MC" = {
      if (par2 == 0){
        return (BiCopGradMMD.MC.1par(
          u1 = u1, u2 = u2, family = family, tau = tau,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          ndrawings = ndrawings)
        )
      } else {
        return(BiCopGradMMD.MC.2par(
          u1 = u1, u2 = u2, family = family, tau = tau, par2 = par2,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          ndrawings = ndrawings)
        )
      }
    },
    "QMCV" = {
      if (is.character(quasiRNG)) {
        switch (
          quasiRNG,

          "sobol" = {quasiRNGFun <- function (n) {
            return (randtoolbox::sobol(n = n, dim = 2))}},

          "halton" = {quasiRNGFun <- function (n) {
            return (randtoolbox::sobol(n = n, dim = 2))}},

          "torus" = {quasiRNGFun <- function (n) {
            return (randtoolbox::torus(n = n, dim = 2))}}
        )
      } else {
        quasiRNGFun <- quasiRNG
      }

      if (par2 == 0){
        return(BiCopGradMMD.QMCV.1par(
          u1 = u1, u2 = u2, family = family, tau = tau,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          quasiRNG = quasiRNGFun, ndrawings = ndrawings)
        )
      } else {
        return(BiCopGradMMD.QMCV.2par.tau(
          u1 = u1, u2 = u2, family = family, tau = tau, par2 = par2,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          quasiRNG = quasiRNGFun, ndrawings = ndrawings)
        )
      }
    }
  )
}


BiCopGradMMD.MC.1par <- function(
  u1, u2, family, tau,
  kernelFun, gamma=0.3, alpha=1, epsilon=0.0001,
  ndrawings=10)
{
  n = length(u1)

  Gfin = 0
  for (j in 1:ndrawings) {
    par = VineCopula::BiCopTau2Par(family, tau)

    U = VineCopula::BiCopSim(N=n, family, par)
    V = VineCopula::BiCopSim(N=n, family, par)

    G = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, deriv="par", log=TRUE)
    G = G*((VineCopula::BiCopTau2Par(family, tau+epsilon) -
              VineCopula::BiCopTau2Par(family, tau)) / epsilon)
    KerCross =
      kernelFun(U[,1], U[,2], V[,1], V[,2], gamma, alpha) -
      kernelFun(u1   , u2   , U[,1], U[,2], gamma, alpha)

    G = mean(G * KerCross)
    Gfin = G/j + Gfin*(j-1)/j
  }

  return(2 * Gfin)
}

BiCopGradMMD.MC.2par <- function(
  u1, u2, family, tau, par2,
  kernelFun, gamma=0.3, alpha=1, epsilon=0.0001,
  ndrawings=10)
{
  n = length(u1)

  Gfin = 0
  for (j in 1:ndrawings) {
    par = VineCopula::BiCopTau2Par(family = family, tau = tau)

    U = VineCopula::BiCopSim(N=n, family, par, par2)
    V = VineCopula::BiCopSim(N=n, family, par, par2)

    G1 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par",  log=TRUE)
    G1 = G1*( (VineCopula::BiCopTau2Par(family, tau+epsilon) -
                 VineCopula::BiCopTau2Par(family, tau)) / epsilon)

    G2 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par2", log=TRUE)
    KerCross =
      kernelFun(U[,1], U[,2], V[,1], V[,2], gamma, alpha) -
      kernelFun(u1   , u2   , U[,1], U[,2], gamma, alpha)

    G1 = mean(G1 * KerCross)
    G2 = mean(G2 * KerCross)
    G = c(G1, G2)
    Gfin = G/j + Gfin*(j-1)/j
  }

  return(2 * Gfin)
}

BiCopGradMMD.QMCV.1par <- function(
  u1, u2, family, tau,
  kernelFun, gamma=0.3, alpha=1, epsilon=0.0001,
  quasiRNG, ndrawings=10)
{
  n = length(u1)
  bivariateGrid = quasiRNG(n = n)

  par = VineCopula::BiCopTau2Par(family,tau)

  V2 = VineCopula::BiCopHinv1(u1 = bivariateGrid[,1], u2 = bivariateGrid[,2],
                              family, par)
  V = cbind(bivariateGrid[,1], V2)

  # We do the Monte Carlo on U
  Gfin = 0
  for (j in 1:ndrawings) {
    par = VineCopula::BiCopTau2Par(family,tau)
    U = VineCopula::BiCopSim(N=n,family,par)

    G = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, deriv="par", log=TRUE)
    G = G*((VineCopula::BiCopTau2Par(family, tau+epsilon) -
              VineCopula::BiCopTau2Par(family, tau)) / epsilon)
    KerCross =
      kernelFun(U[,1], U[,2], V[,1], V[,2], gamma, alpha) -
      kernelFun(u1   , u2   , U[,1], U[,2], gamma, alpha)

    G = mean(G * KerCross)
    Gfin = G/j + Gfin*(j-1)/j
  }

  return(2 * Gfin)
}


BiCopGradMMD.QMCV.2par.tau <- function(
  u1, u2, family, tau, par2,
  kernelFun, gamma=0.3, alpha=1, epsilon=0.0001,
  quasiRNG, ndrawings=10)
{
  n = length(u1)
  bivariateGrid = quasiRNG(n = n)
  par = VineCopula::BiCopTau2Par(family,tau)
  V2 = VineCopula::BiCopHinv1(u1 = bivariateGrid[,1], u2 = bivariateGrid[,2],
                              family, par, par2)
  V = cbind(bivariateGrid[,1], V2)

  # We do the Monte Carlo on U
  Gfin = 0
  for (j in 1:ndrawings) {
    U = VineCopula::BiCopSim(N=n, family, par, par2)

    G1 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par" , log=TRUE)
    # G1 = G1*((VineCopula::BiCopTau2Par(family, tau+epsilon) -
    #             VineCopula::BiCopTau2Par(family, tau)) / epsilon)

    G2 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par2", log=TRUE)
    KerCross =
      kernelFun(U[,1], U[,2], V[,1], V[,2], gamma, alpha) -
      kernelFun(u1   , u2   , U[,1], U[,2], gamma, alpha)

    G1 = mean(G1 * KerCross)
    G2 = mean(G2 * KerCross)
    G = c(G1, G2)

    Gfin = G/j + Gfin*(j-1)/j
  }

  return(2 * Gfin)
}

BiCopGradMMD.QMCV.2par.par <- function(
  u1, u2, family, par, par2,
  kernelFun, gamma=0.3, alpha=1, epsilon=0.0001,
  quasiRNG, ndrawings=10)
{
  n = length(u1)
  bivariateGrid = quasiRNG(n = n)
  V2 = VineCopula::BiCopHinv1(u1 = bivariateGrid[,1], u2 = bivariateGrid[,2],
                              family, par, par2)
  V = cbind(bivariateGrid[,1], V2)

  # We do the Monte Carlo on U
  Gfin = 0
  for (j in 1:ndrawings) {
    U = VineCopula::BiCopSim(N=n, family, par, par2)

    G1 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par" , log=TRUE)
    # G1 = G1*((VineCopula::BiCopTau2Par(family, tau+epsilon) -
    #             VineCopula::BiCopTau2Par(family, tau)) / epsilon)

    G2 = VineCopula::BiCopDeriv(U[,1], U[,2], family, par, par2, deriv="par2", log=TRUE)
    KerCross =
      kernelFun(U[,1], U[,2], V[,1], V[,2], gamma, alpha) -
      kernelFun(u1   , u2   , U[,1], U[,2], gamma, alpha)

    G1 = mean(G1 * KerCross)
    G2 = mean(G2 * KerCross)
    G = c(G1, G2)

    Gfin = G/j + Gfin*(j-1)/j
  }

  return(2 * Gfin)
}

