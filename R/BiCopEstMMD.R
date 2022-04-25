
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
#'   If \code{NULL}, a random value is chosen instead.
#'
#' @param par if different from \code{NULL}, the parameter \code{tau} is ignored,
#'   and the initial parameter must be given here.
#'   The initial Kendall's tau is then computed thanks to
#'   \code{VineCopula::\link[VineCopula]{BiCopPar2Tau}()}.
#'
#' @param par2 initial value for the second parameter, if any. (Works only for Student copula).
#'
#' @param gamma parameter \eqn{\gamma} to be used in the kernel.
#' If \code{gamma="default"}, a default value is used.
#'
#' @param niter the stochastic gradient algorithm is composed of two phases:
#' a first "burn-in" phase and a second "averaging" phase.
#' If \code{niter} is of size \code{1}, the same number of iterations is used for
#' both phases of the stochastic gradient algorithm. If \code{niter} is of size \code{2},
#' then \code{niter[1]} iterations are done for the burn-in phase and \code{niter[2]}
#' for the averaging phase.
#'
#' @param C_eta a multiplicative constant controlling for the size of the gradient descent step.
#' The step size is then computed as \code{C_eta / sqrt(i_iter)}
#' where \code{i_iter} is the index of the current iteration of the stochastic gradient algorithm.
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
#' @references Alquier, P., Chérief-Abdellatif, B.-E., Derumigny, A., and Fermanian, J.D. (2022).
#' Estimation of copulas via Maximum Mean Discrepancy.
#' Journal of the American Statistical Association, \doi{10.1080/01621459.2021.2024836}.
#'
#'
#' @examples
#' # Estimation of a bivariate Gaussian copula with correlation 0.5.
#' dataSampled = VineCopula::BiCopSim(N = 500, family = 1, par = 0.5)
#' estimator = BiCopEstMMD(u1 = dataSampled[,1], u2 = dataSampled[,2], family = 1, niter = 10)
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
  family, tau = NULL, par = NULL, par2 = NULL,
  kernel = "gaussian", gamma = "default", alpha = 1,
  niter = 100, C_eta = 1, epsilon = 0.0001,
  method = "QMCV", quasiRNG = "sobol", ndrawings = 10)
{
  # Checking for input validity
  verifData(u1, u2)

  if (!is.finite(alpha)){stop("Finite value for 'alpha' required.")}
  if (!is.finite(epsilon)){stop("Finite value for 'epsilon' required.")}
  if (!is.finite(niter)){stop("Finite value for 'niter' required.")}
  if (!is.finite(C_eta)){stop("Finite value for 'C_eta' required.")}
  if (!is.finite(ndrawings)){stop("Finite value for 'ndrawings' required.")}

  # Choice of the kernel
  if (is.character(kernel))
  {
    kernelFun <- findKernelFunction(kernel)
  } else {
    kernelFun <- kernel
  }

  # Choice of gamma
  if (gamma == "default"){
    if (family == 1 || family == 5){
      if (kernel %in% c("gaussian", "exp-l2", "exp-l1", "inv-l2", "inv-l1")){
        gamma = 0.25
      } else {
        gamma = 0.8
      }
    } else {
      if (kernel %in% c("gaussian", "exp-l2", "exp-l1", "inv-l2", "inv-l1")){
        gamma = 0.1
      } else {
        gamma = 0.4
      }
    }
  } else {
    if (!is.finite(gamma)){stop("Finite value for 'gamma' required.")}
  }

  # If only one number of iterations is given,
  # it is reused for the burn-in phase and the averaging phase
  niter = rep(niter, length.out = 2)
  if (is.null(tau)){
    # We try a random guess
    # tau = wdm::wdm(u1, u2, method = "kendall")
    if (family %in% c(1, 5, 6,16,26,36)){
      tau = stats::runif(n = 1, min = -0.95, max = 0.95)
    } else if (family %in% c(3, 23, 4, 24)){
      tau = stats::runif(n = 1, min = 0.05, max = 0.95)
    } else if (family %in% c(13, 33, 14, 34)){
      tau = stats::runif(n = 1, min = -0.95, max = -0.05)
    }
  }
  switch (
    method,

    "MC" = {
      if (family %in% c(1, 3,13,23,33, 4,14,24,34, 5, 6,16,26,36)){
        if (!is.null(par)) {
          tau = VineCopula::BiCopPar2Tau(family = family, par = par, par2 = par2)
        }

        # 1- Burn-in phase
        tauIter = tau
        for (i_iter in 1:niter[1]){
          tauIter = tauIter -
            C_eta * BiCopGradMMD.MC.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(i_iter)
        }
        # 2- Averaging phase
        tauIter_vec = rep(NA, niter[2]+1)
        tauIter_vec[1] = tauIter
        for (i_iter in 1:niter[2]){
          tauIter_vec[i_iter+1] = tauIter_vec[i_iter] -
            C_eta * BiCopGradMMD.MC.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter_vec[i_iter],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(niter[1] + i_iter)
        }
        estim = VineCopula::BiCop(family, tau = mean(tauIter_vec))

      } else {
        warning("The estimation can be unstable for the 2-parameters families.")
        # 1- Burn-in phase
        paramIter = c(tau, par2)
        for (i_iter in 1:niter[1]){
          paramIter = paramIter -
            C_eta * BiCopGradMMD.MC.2par(
              u1 = u1, u2 = u2, family = family, tau = paramIter[1], par2 = paramIter[2],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(i_iter)
        }
        # 2- Averaging phase
        paramIter_mat = matrix(ncol = 2, nrow = niter[2]+1)
        paramIter_mat[1,] = paramIter
        for (i_iter in 1:niter[2] ){
          paramIter_mat[i_iter+1,] = paramIter_mat[i_iter,] -
            C_eta * BiCopGradMMD.MC.2par(
              u1 = u1, u2 = u2, family = family,
              tau = paramIter_mat[i_iter,1], par2 = paramIter_mat[i_iter,2],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              ndrawings = ndrawings) / sqrt(niter[1] + i_iter)
        }
        estim = VineCopula::BiCop(family, tau = mean(paramIter_mat[,1]),
                                  par2 = mean(paramIter_mat[,2]) )
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
        if (!is.null(par)) {
          tau = VineCopula::BiCopPar2Tau(family = family, par = par, par2 = par2)
        }
        if (family == 1){
          tauMin = -0.999
          tauMax = 0.999
        } else if (family == 5){
          tauMin = -0.96
          tauMax = 0.96
        } else if (family %in% c(3,23,4,24,6,26) ){
          tauMin = 0.001
          tauMax = 0.999
        } else if (family %in% c(13,33,14,34,16,36) ){
          tauMin = -0.999
          tauMax = -0.001
        }

        # 1- Burn-in phase
        tauIter = tau
        for (i_iter in 1:niter[1]){
          tauIter = tauIter -
            C_eta * BiCopGradMMD.QMCV.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter,
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              quasiRNG = quasiRNGFun, ndrawings = ndrawings) / sqrt(i_iter)
          tauIter = max(tauMin, min(tauIter, tauMax)) # Constraining the range of the tau
        }
        # 2- Averaging phase
        tauIter_vec = rep(NA, niter[2]+1)
        tauIter_vec[1] = tauIter
        for (i_iter in 1:niter[2]){
          tauIter_vec[i_iter+1] = tauIter_vec[i_iter] -
            C_eta * BiCopGradMMD.QMCV.1par(
              u1 = u1, u2 = u2, family = family, tau = tauIter_vec[i_iter],
              kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
              quasiRNG = quasiRNGFun, ndrawings = ndrawings) / sqrt(niter[1] + i_iter)
          tauIter_vec[i_iter+1] = max(tauMin, min(tauIter_vec[i_iter+1], tauMax))
        }
        estim = VineCopula::BiCop(family, tau = mean(tauIter_vec))

      } else if (family == 2) {
        warning("The estimation can be unstable for the Student family.")
        estim = BiCopEstMMD.QMCV.student(
          u1 = u1, u2 = u2,
          kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
          quasiRNGFun = quasiRNGFun, ndrawings = ndrawings, niter = niter, C_eta = C_eta)

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
  quasiRNGFun, ndrawings, niter, C_eta)
{
  par = VineCopula::BiCopTau2Par(family = 2, tau = wdm::wdm(u1, u2, method = "kendall"))

  # 1- Burn-in phase
  paramIter = c(par, par2)

  for (i_iter in 1:niter[1]){
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

  # 2- Averaging phase
  paramIter_mat = matrix(ncol = 2, nrow = niter[2]+1)
  paramIter_mat[1,] = paramIter
  for (i_iter in 1:niter[2]){
    if (paramIter_mat[i_iter, 2] < 6) {alpha2 = 1000
    } else {alpha2 = 10000}

    paramIter_mat[i_iter+1,] = paramIter_mat[i_iter,] -
      C_eta * c(5,alpha2) *
      BiCopGradMMD.QMCV.2par.par(
        u1 = u1, u2 = u2, family = 2, par = paramIter[1], par2 = paramIter[2],
        kernelFun = kernelFun, gamma = gamma, alpha = alpha, epsilon = epsilon,
        quasiRNG = quasiRNGFun, ndrawings = ndrawings) / sqrt(niter[1] + i_iter)

    if (paramIter_mat[i_iter+1, 2] < 2.0001) {paramIter_mat[i_iter+1, 2] = 2.001}
    if (paramIter_mat[i_iter+1, 1] < -0.9999) {paramIter_mat[i_iter+1, 1] = -0.9999}
    if (paramIter_mat[i_iter+1, 1] > 0.9999) {paramIter_mat[i_iter+1, 1] = 0.9999}
  }

  estim = VineCopula::BiCop(family = 2,
                            par = mean(paramIter_mat[,1]), par2 = mean(paramIter_mat[,2]))
  return (estim)
}


