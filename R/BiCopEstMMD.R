
#' Estimation of parametric bivariate copulas using
#' stochastic gradient descent on the MMD criteria
#'
#' This function uses computes the MMD-estimator of a bivariate copula family.
#' This computation is done through a stochastic gradient algorithm,
#' that is itself computed by the function \code{\link{BiCopGradMMD}()}.
#' The main arguments are the two vectors of observations, and the copula family.
#' The bidimensional copula families are indexed in the same way as
#' in \code{VineCopula::\link[VineCopula]{BiCop}()} (which computes the MLE estimator).
#'
#' @param tau the copula family can be parametrized by the parameter \code{par}
#'   or by Kendall's tau.
#'   Here, the user can choose the initial value of tau for the stochastic gradient algorithm.
#'
#' @param par if different from \code{NULL}, the parameter \code{tau} is ignored,
#'   and the initial parameter must be given here.
#'   The initial Kendall's tau is then computed thanks to
#'   \code{VineCopula::\link[VineCopula]{BiCopPar2Tau}()}.
#'
#' @param par2 initial value for the second parameter, if any. (Works only for Student copula).
#'
#' @param niter number of iterations of the stochastic gradient algorithm.
#'
#'
#' @inheritParams BiCopGradMMD
#'
#' @return an object of class \code{VineCopula::\link[VineCopula]{BiCop}()}
#' containing the estimated copula.
#'
#'
#' @seealso \code{VineCopula::\link[VineCopula]{BiCopEst}()} for other methods of estimation
#'   such as Maximum Likelihood Estimation or Inversion of Kendall's tau.
#' \code{\link{BiCopGradMMD}()} for the computation of the stochastic gradient.
#' \code{\link{BiCopEst.MO}} for the estimation of Marshall-Olkin copulas by MMD.
#'
#' @references Alquier, P., Chérief-Abdellatif, B.-E., Derumigny, A., and Fermanian, J.D. (2020).
#' Estimation of copulas via Maximum Mean Discrepancy.
#' ArXiv preprint \href{https://arxiv.org/abs/2010.00408}{arxiv:2010.00408}
#'
#'
#' @examples
#' # Estimation of a bivariate Gaussian copula with correlation 0.5.
#' dataSampled = VineCopula::BiCopSim(N = 500, family = 1, par = 0.5)
#' estimator = BiCopEstMMD(u1 = dataSampled[,1], u2 = dataSampled[,2], family = 1, niter=10)
#' estimator$par
#'
#' \donttest{
#' # Estimation of a bivariate Student copula with correlation 0.5 and 5 degrees of freedom
#' dataSampled = VineCopula::BiCopSim(N = 1000, family = 2, par = 0.5, par2 = 5)
#' estimator = BiCopEstMMD(u1 = dataSampled[,1], u2 = dataSampled[,2], family = 2)
#' estimator$par
#' estimator$par2
#'
#'
#' # Comparison with maximum likelihood estimation with and without outliers
#' dataSampled = VineCopula::BiCopSim(N = 500, family = 1, par = 0.5)
#' estimatorMMD = BiCopEstMMD(u1 = dataSampled[,1], u2 = dataSampled[,2], family = 1)
#' estimatorMMD$par
#' estimatorMLE = VineCopula::BiCopEst(u1 = dataSampled[,1], u2 = dataSampled[,2],
#'   family = 1, method = "mle")
#' estimatorMLE$par
#'
#' dataSampled[1:10,1] = 0.999
#' dataSampled[1:10,2] = 0.001
#' estimatorMMD = BiCopEstMMD(u1 = dataSampled[,1], u2 = dataSampled[,2], family = 1)
#' estimatorMMD$par
#' estimatorMLE = VineCopula::BiCopEst(u1 = dataSampled[,1], u2 = dataSampled[,2],
#'   family = 1, method = "mle")
#' estimatorMLE$par
#'
#'
#' # Estimation of a bivariate Gaussian copula with real data
#' data("daxreturns", package = "VineCopula")
#' BiCopEstMMD(u1 = daxreturns[,1], u2 = daxreturns[,2], family = 1)
#' estimator$par
#' }
#'
#'
#' @export
#'
#'
BiCopEstMMD <- function(
  u1, u2,
  family, tau = 0.1, par = NULL, par2 = NULL,
  kernel = "gaussian", gamma=0.23, alpha=1, niter=100, epsilon=0.0001,
  method = "QMCV", quasiRNG = "sobol", ndrawings=10)
{
  verifData(u1, u2)

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
      if (family %in% c(1, 3,13,23,33, 4,14,24,34, 5, 6,16,26,36)){
        if (!is.null(par))
        {
          tau = VineCopula::BiCopPar2Tau(family = family, par = par, par2 = par2)
        }

        tauIter = tau
        for (i_iter in 1:niter){
          tauIter = tauIter -
            BiCopGradMMD.MC.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(i_iter)
        }
        for (i_iter in (niter+1):(2*niter)){
          tauIter = tauIter -
            BiCopGradMMD.MC.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / (sqrt(i_iter) * i_iter)
        }
        estim = VineCopula::BiCop(family, tau = tauIter)

      } else {

        paramIter = c(tau, par2)
        for (i_iter in 1:niter){
          paramIter = paramIter -
            BiCopGradMMD.MC.2par(
              u1 = u1, u2 = u2, family = family, tau = paramIter[1], par2 = paramIter[2],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(i_iter)
        }
        for (i_iter in (niter+1):(2*niter)){
          paramIter = paramIter -
            BiCopGradMMD.MC.2par(
              u1 = u1, u2 = u2, family = family, tau = paramIter[1], par2 = paramIter[2],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / (sqrt(i_iter) * i_iter)
        }
        estim = VineCopula::BiCop(family, tau = paramIter[1], par2 = paramIter[2])
      }

    },

    "QMCV" = {
      if (is.character(quasiRNG)) {
        switch (
          quasiRNG,

          "sobol" = {quasiRNGFun <- function (n) {
            return (randtoolbox::sobol(n = n, dim = 2))}},

          "halton" = {quasiRNGFun <- function (n) {
            return (randtoolbox::halton(n = n, dim = 2))}},

          "torus" = {quasiRNGFun <- function (n) {
            return (randtoolbox::torus(n = n, dim = 2))}}
        )
      } else {
        quasiRNGFun <- quasiRNG
      }

      if (family %in% c(1, 3,13,23,33, 4,14,24,34, 5, 6,16,26,36)){
        if (!is.null(par))
        {
          tau = VineCopula::BiCopPar2Tau(family = family, par = par, par2 = par2)
        }

        tauIter = tau
        for (i_iter in 1:niter){
          tauIter = tauIter -
            BiCopGradMMD.QMCV.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              quasiRNG = quasiRNGFun, ndrawings = ndrawings) / sqrt(i_iter)
          if (family %in% c(3,23,4,24,6,26) & tauIter < 0.001){ tauIter = 0.001}
        }
        for (i_iter in (niter+1):(2*niter)){
          tauIter = tauIter -
            BiCopGradMMD.QMCV.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              quasiRNG = quasiRNGFun, ndrawings = ndrawings) / (sqrt(i_iter) * i_iter)
          if (family %in% c(3,23,4,24,6,26) & tauIter < 0.001){ tauIter = 0.001}
        }
        estim = VineCopula::BiCop(family, tau = tauIter)

      } else if (family == 2) {
        warning("The estimation can be unstable for the Student family.")
        estim = BiCopEstMMD.QMCV.student(
          u1 = u1, u2 = u2,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          quasiRNGFun = quasiRNGFun, ndrawings = ndrawings, niter = niter)

      } else {
        stop("family not implemented yet in BiCopEstMMD.")
      }
    }
  )


  estim$call <- match.call()
  return(estim)
}


BiCopEstMMD.QMCV.student <- function(
  u1, u2, par2 = 8,
  kernelFun, gamma, alpha, epsilon,
  quasiRNGFun, ndrawings, niter)
{
  par = VineCopula::BiCopTau2Par(family = 2, tau = pcaPP::cor.fk(u1, u2))

  paramIter = c(par, par2)

  for (i_iter in 1:niter){
    if        (paramIter[2] < 5) {alpha2 =  1000
    } else if (paramIter[2] < 6) {alpha2 =  4000
    } else                       {alpha2 = 15000 }

    paramIter = paramIter - c(5,alpha2) *
      BiCopGradMMD.QMCV.2par.par(
        u1 = u1, u2 = u2, family = 2, par = paramIter[1], par2 = paramIter[2],
        kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
        quasiRNG = quasiRNGFun, ndrawings = ndrawings) / sqrt(i_iter)
    # print(paramIter)

    if (paramIter[2] < 2.0001) {paramIter[2] = 2.001}
    if (paramIter[1] < -0.9999) {paramIter[1] = -0.9999}
    if (paramIter[1] > 0.9999) {paramIter[1] = 0.9999}
  }

  for (i_iter in (niter+1):(2*niter)){
    if (paramIter[2] < 6) {alpha2 = 1000
    } else {alpha2 = 10000}

    paramIter = paramIter - c(5,alpha2) *
      BiCopGradMMD.QMCV.2par.par(
        u1 = u1, u2 = u2, family = 2, par = paramIter[1], par2 = paramIter[2],
        kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
        quasiRNG = quasiRNGFun, ndrawings = ndrawings) / (sqrt(i_iter) * i_iter)

    if (paramIter[2] < 2.0001) {paramIter[2] = 2.001}
    if (paramIter[1] < -0.9999) {paramIter[1] = -0.9999}
    if (paramIter[1] > 0.9999) {paramIter[1] = 0.9999}
  }

  estim = VineCopula::BiCop(family = 2, par = paramIter[1], par2 = paramIter[2])
  return (estim)
}


