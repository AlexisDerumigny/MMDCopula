

#' Simulation of Marshall-Olkin copula
#'
#' This functions simulates independent realizations
#' from the Marshall-Olkin copula.
#'
#' @param n number of samples
#' @param alpha parameter of the Marshall-Olkin copula
#'
#' @return an \eqn{n \times 2} matrix containing the samples
#'
#' @seealso \code{\link{BiCopEst.MO}} for the estimation of
#' Marshall-Olkin copulas.
#'
#' @examples
#' # Simulation from a Marshall-Olkin copula with parameter alpha = 0.5
#' BiCopSim.MO(n = 100, alpha = 0.5)
#'
#'
#' @export
#'
BiCopSim.MO <- function(n, alpha) {

  ST = matrix(data = NA, nrow = n, ncol = 2)
  ST[,1] = stats::rbeta(n = n, shape1 = 2-alpha, shape2 = 1)
  ST[,2] = stats::runif(n = n, 0, ST[,1])
  P = alpha/(2-alpha)

  doEqual = stats::runif(n,0,1) < P
  ST[,2] = ifelse(doEqual, ST[,1], ST[,2])

  doExchange = (stats::runif(n,0,1) < 0.5) & !doEqual
  AUX = ST[,2]
  ST[,2] = ifelse(doExchange, ST[,1], ST[,2])
  ST[,1] = ifelse(doExchange, AUX, ST[,1])

  return(ST)
}


#' Estimation of Marshall-Olkin copulas
#'
#' @param u1 vector of observations of the first coordinate, in \eqn{[0,1]}.
#' @param u2 vector of observations of the second coordinate, in \eqn{[0,1]}.
#'
#' @param method a character giving the name of the estimation method, among:
#'   \itemize{
#'     \item \code{curve}: \eqn{\alpha} is estimated by inversion of
#'       the probability measure of the diagonal
#'       \eqn{\{(u,v): u = v\}}{ {(u,v): u = v} }
#'     \item \code{itau}: \eqn{\alpha} is estimated by inversion of Kendall's tau
#'     \item \code{MMD}: \eqn{\alpha} is estimated by MMD optimization
#'   }
#'
#' @param par.start starting parameter of the gradient descent.
#' (only used for \code{method = "MMD"})
#'
#' @param kernel the kernel used in the MMD distance
#' (only used for \code{method = "MMD"}) :
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
#'  Each of these names can receive the suffix ".Phi", such as "gaussian.Phi"
#'  to indicates that the kernel \eqn{k(x,y)} is replaced by
#'  \eqn{k(\Phi^{-1}(x) , \Phi^{-1}(y))} where \eqn{\Phi^{-1}} denotes the quantile
#'  function of the standard Normal distribution.
#'
#' @param gamma parameter \eqn{\gamma} to be used in the kernel.
#' (only used for \code{method = "MMD"})
#'
#' @param alpha parameter \eqn{\alpha} to be used in the kernel, if any.
#' (only used for \code{method = "MMD"})
#'
#' @param ndrawings number of replicas of the stochastic estimate of the gradient
#' drawn at each step. The gradient is computed using the average of these replicas.
#' (only used for \code{method = "MMD"})
#'
#' @param niter number of iterations of the stochastic gradient algorithm.
#' (only used for \code{method = "MMD"})
#'
#' @param naveraging number of full run of the stochastic gradient algorithm
#' that are averaged at the end to give the final estimated parameter.
#' (only used for \code{method = "MMD"})
#'
#' @return the estimated parameter (\code{alpha}) of the Marshall-Olkin copula.
#'
#' @seealso \code{\link{BiCopSim.MO}} for the estimation of
#' Marshall-Olkin copulas.
#' \code{\link{BiCopEstMMD}} for the estimation of other parametric copula families by MMD.
#'
#'
#' @references Alquier, P., Chérief-Abdellatif, B.-E., Derumigny, A., and Fermanian, J.D. (2020).
#' Estimation of copulas via Maximum Mean Discrepancy.
#' ArXiv preprint \href{https://arxiv.org/abs/2010.00408}{arxiv:2010.00408}
#'
#' @examples
#' U <- BiCopSim.MO(n = 1000, alpha = 0.2)
#' estimatedPar <- BiCopEst.MO(u1 = U[,1], u2 = U[,2], method = "MMD", niter = 1, ndrawings = 1)
#' \donttest{
#' estimatedPar <- BiCopEst.MO(u1 = U[,1], u2 = U[,2], method = "MMD")
#' }
#'
#' @export
#'
BiCopEst.MO <- function(
  u1, u2, method,
  par.start = 0.5, kernel = "gaussian.Phi",
  gamma=0.95, alpha=1,
  niter=100, ndrawings=10, naveraging = 1,
  methodMC = "MC")
{
  verifData(u1, u2)

  switch (
    method,

    "curve" = {
      estimator = BiCopEst.MO.curve(u1, u2)
    },

    "itau" = {
      result = BiCopEst.MO.itau(u1, u2)
      return(result)
    },


    "MMD" = {

      # Choice of the kernel
      if (is.character(kernel))
      {
        kernelFun <- findKernelFunction(kernel)
      } else {
        kernelFun <- kernel
      }

      if (methodMC == "MC"){
        estimator = BiCopEst.MO.MMD.MC(
          u1 = u1, u2 = u2, par.start = par.start,
          kernelFun = kernelFun,
          gamma = gamma, alpha = alpha,
          niter = niter, ndrawings = ndrawings, naveraging = naveraging)
      }
      # else if (methodMC == "QMCV"){
      #
      #   if (is.character(quasiRNG)) {
      #     switch (
      #       quasiRNG,
      #
      #       "sobol" = {quasiRNGFun <- function (n) {
      #         return (randtoolbox::sobol(n = n, dim = 2))}},
      #
      #       "halton" = {quasiRNGFun <- function (n) {
      #         return (randtoolbox::halton(n = n, dim = 2))}},
      #
      #       "torus" = {quasiRNGFun <- function (n) {
      #         return (randtoolbox::torus(n = n, dim = 2))}}
      #     )
      #   } else {
      #     quasiRNGFun <- quasiRNG
      #   }
      #
      #   estimator = BiCopEst.MO.MMD.QMCV(
      #     u1 = u1, u2 = u2, par.start = par.start,
      #     kernelFun = kernelFun,
      #     gamma = gamma, alpha = alpha,
      #     niter = niter, ndrawings = ndrawings, naveraging = naveraging,
      #     quasiRNGFun = quasiRNGFun)
      # }

    }

  )

  tau = estimator / (2 - estimator)
  return (list(tau = tau, par = estimator))

}

# Estimation of Marshall Olkin copulas using the method of moments
# based on the probability measure of the diagonal
BiCopEst.MO.curve = function(u1,u2)
{
  proba = mean(u1==u2)
  return(2*proba/(1+proba))
}


# Estimation of Marshall Olkin copulas using the method of moments
# based on the inversion of Kendall's tau
BiCopEst.MO.itau = function(u1,u2)
{
  tau = pcaPP::cor.fk(u1, u2)
  return(list(tau = tau, par = 2*tau/(1+tau)))
}


# Estimation of Marshall-Olkin copulas by MMD
BiCopEst.MO.MMD.MC = function(
  u1, u2, par.start = 0.5, kernelFun,
  gamma = 0.3, alpha = 1,
  niter = 100, ndrawings = 10, naveraging = 1)
{

  n = length(u1)
  estimatorsA = rep(NA, naveraging)

  for (i in 1:naveraging){
    aIter = par.start
    for (i_iter in 1:niter){
      Grad = 0
      for (j in 1:ndrawings){
        U = BiCopSim.MO(n,aIter)
        V = BiCopSim.MO(n,aIter)
        usup = (U[,1]>U[,2])
        uinf = (U[,1]<U[,2])
        ueq = (U[,1]==U[,2])
        dlog = (log(U[,1])-aIter/(1-aIter))*usup +
          (log(U[,2])-aIter/(1-aIter))*uinf +
          (1/aIter-log(U[,1])) * ueq

        grad = mean(2 * ( kernelFun(U[,1], U[,2], V[,1], V[,2],gamma,alpha)
                          - kernelFun( u1,    u2, U[,1], U[,2],gamma,alpha) ) * dlog)
        Grad = grad/j + Grad*(j-1)/j
      }
      aIter = aIter - Grad / sqrt(i_iter)
    }
    for (i_iter in (niter+1):(2*niter)){
      Grad = 0
      for (j in 1:ndrawings){
        U = BiCopSim.MO(n,aIter)
        V = BiCopSim.MO(n,aIter)
        usup = (U[,1]>U[,2])
        uinf = (U[,1]<U[,2])
        ueq = (U[,1]==U[,2])
        dlog = (log(U[,1])-aIter/(1-aIter))*usup +
          (log(U[,2])-aIter/(1-aIter))*uinf +
          (1/aIter-log(U[,1])) * ueq

        grad = mean(2 * ( kernelFun(U[,1], U[,2], V[,1], V[,2],gamma,alpha)
                          - kernelFun( u1,    u2, U[,1], U[,2],gamma,alpha) ) * dlog)
        Grad = grad/j + Grad*(j-1)/j
      }
      aIter = aIter - Grad / (sqrt(i_iter) * i_iter)
    }
    estimatorsA[i] = aIter
  }

  estim = mean( estimatorsA )

  return(estim)
}

#
# # Estimation of Marshall-Olkin copulas by MMD
# BiCopEst.MO.MMD.QMCV = function(
#   u1, u2, par.start=0.5, kernelFun,
#   gamma=0.3, alpha=1,
#   quasiRNGFun, niter=100, ndrawings=10, naveraging = 1)
# {
#
#   n = length(u1)
#   estimatorsA = rep(NA, naveraging)
#
#   for (i in 1:naveraging){
#     aIter = par.start
#     for (i_iter in 1:niter){
#       Grad = 0
#       for (j in 1:ndrawings){
#         U = BiCopSimMO(n,aIter)
#         V = BiCopSimMO(n,aIter)
#         usup = (U[,1]>U[,2])
#         uinf = (U[,1]<U[,2])
#         ueq = (U[,1]==U[,2])
#         dlog = (log(U[,1])-aIter/(1-aIter))*usup +
#           (log(U[,2])-aIter/(1-aIter))*uinf +
#           (1/aIter-log(U[,1])) * ueq
#
#         grad = mean(2 * ( kernelFun(U[,1], U[,2], V[,1], V[,2], kernel,gamma,alpha)
#                           - kernelFun( u1,    u2, U[,1], U[,2], kernel,gamma,alpha) ) * dlog)
#         Grad = grad/j + Grad*(j-1)/j
#       }
#       aIter = aIter - Grad / sqrt(i_iter)
#     }
#     for (i_iter in (niter+1):(2*niter)){
#       Grad = 0
#       for (j in 1:ndrawings){
#         U = BiCopSimMO(n,aIter)
#         V = BiCopSimMO(n,aIter)
#         usup = (U[,1]>U[,2])
#         uinf = (U[,1]<U[,2])
#         ueq = (U[,1]==U[,2])
#         dlog = (log(U[,1])-aIter/(1-aIter))*usup +
#           (log(U[,2])-aIter/(1-aIter))*uinf +
#           (1/aIter-log(U[,1])) * ueq
#
#         grad = mean(2 * ( kernelFun(U[,1], U[,2], V[,1], V[,2], kernel,gamma,alpha)
#                           - kernelFun( u1,    u2, U[,1], U[,2], kernel,gamma,alpha) ) * dlog)
#         Grad = grad/j + Grad*(j-1)/j
#       }
#       aIter = aIter - Grad / (sqrt(i_iter) * i_iter)
#     }
#     estimatorsA[i] = aIter
#   }
#
#   estim = mean( estimatorsA )
#
#   return(estim)
# }


