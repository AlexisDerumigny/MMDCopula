

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
#'
#' @export
#'
BiCopSim.MO <- function(n, alpha) {

  ST = matrix(data = NA, nrow = n, ncol = 2)
  ST[,1] = rbeta(n = n, shape1 = 2-alpha, shape2 = 1)
  ST[,2] = runif(n = n, 0, ST[,1])
  P = alpha/(2-alpha)

  for (i in 1:n)
  {
    if (runif(1,0,1) < P) {
      ST[i,2] = ST[i,1]
    }
    else if (runif(1,0,1) < 0.5) {
      AUX = ST[i,2]; ST[i,2] = ST[i,1]; ST[i,1] = AUX
    }

  }
  return(ST)
}


#' Estimation of Marshall-Olkin copulas
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
#'
#' @export
#'
BiCopEst.MO <- function(
  u1, u2, method,
  par=0.5, kernel = "exp-l2",
  gamma=0.3, alpha=1,
  niter=100, ndrawings=10, naveraging = 1)
{
  verifData(u1, u2)

  switch (
    method,

    "curve" = {
      estimator = BiCopEst.MO.curve(u1, u2)
    },

    "itau" = {
      estimator = BiCopEst.MO.itau(u1, u2)
    },


    "MMD" = {

      # Choice of the kernel
      if (is.character(kernel))
      {
        kernelFun <- findKernelFunction(kernel)
      } else {
        kernelFun <- kernel
      }

      estimator = BiCopEst.MO.MMD.MC(
        u1 = u1, u2 = u2, par = par, kernelFun = kernelFun,
        gamma = gamma, alpha = alpha,
        niter = niter, ndrawings = ndrawings, naveraging = naveraging)
    }

  )

  return (estimator)

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
  return(2*tau/(1+tau))
}


# Estimation of Marshall-Olkin copulas by MMD
BiCopEst.MO.MMD.MC = function(
  u1, u2, par=0.5, kernelFun,
  gamma=0.3, alpha=1,
  niter=100, ndrawings=10, naveraging = 1)
{

  n = length(u1)
  estimatorsA = rep(NA, naveraging)

  for (i in 1:naveraging){
    aIter = par
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

  # # The averaging is done on the Kendall's tau and not on the parameter
  estim = mean( estimatorsA[i] )

  return(estim)
}

#
# # Estimation of Marshall-Olkin copulas by MMD
# BiCopEst.MO.MMD.QMCV = function(
#   u1, u2, par=0.5, kernelFun,
#   gamma=0.3, alpha=1,
#   quasiRNGFun, niter=100, ndrawings=10, naveraging = 1)
# {
#
#   n = length(u1)
#   estimatorsA = rep(NA, naveraging)
#
#   for (i in 1:naveraging){
#     aIter = par
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
#   # # The averaging is done on the Kendall's tau and not on the parameter
#   estim = mean( estimatorsA[i] )
#
#   return(estim)
# }


