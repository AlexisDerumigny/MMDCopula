

#' Confidence intervals for the estimated parameter
#' of a bivariate parametric copula using MMD estimation
#'
#' @param x1 vector of observations of the first coordinate.
#' @param x2 vector of observations of the second coordinate.
#'
#' @param family parametric family of copulas.
#' Supported families are: \itemize{
#'   \item \code{1}: Gaussian copulas
#'   \item \code{3}: Clayton copulas
#'   \item \code{4}: Gumbel copulas
#'   \item \code{5}: Frank copulas
#'   \item \code{MO}: Marshall-Olkin copulas
#'   }
#'
#' @param nResampling number of resampling times.
#'
#' @param subsamplingSize size of the subsample.
#' By default it is \code{length(u1)},
#' i.e. this corresponds to the nonparametric boostrap.
#'
#' @param corrSubSampling this parameter is only used for subsampling-based confidence intervals.
#' If \code{TRUE}, the confidence interval uses the corrected subsample empirical process.
#'
#' @param level the nominal confidence level.
#'
#' @param ... other parameters to be given to \code{\link{BiCopEstMMD}}
#' or \code{\link{BiCopEst.MO}}.
#'
#' @return a list with the confidence intervals CI.Tau for Kendall's tau
#' and CI.Par for the corresponding parameter.
#'
#' @references Alquier, P., Chérief-Abdellatif, B.-E., Derumigny, A., and Fermanian, J.D. (2020).
#' Estimation of copulas via Maximum Mean Discrepancy.
#' ArXiv preprint \href{https://arxiv.org/abs/2010.00408}{arxiv:2010.00408}.
#'
#' Kojadinovic I., and Stemikovskaya, K. (2019)
#' Subsampling (weighted smooth) empirical copula processes.
#' Journal of Multivariate Analysis, 173, 704-723,
#' \doi{10.1016/j.jmva.2019.05.007}.
#'
#' @examples
#' data = VineCopula::BiCopSim(N = 50, family = 1, par = 0.3)
#' result = BiCopConfIntMMD(x1 = data[,1], x2 = data[,2], family = 1,
#'   nResampling = 2, subsamplingSize = 10, niter = 10)
#' \donttest{
#' data_ = VineCopula::BiCopSim(N = 1000, family = 1, par = 0.3)
#' result_ = BiCopConfIntMMD(x1 = data_[,1], x2 = data_[,2], family = 1)
#' result_$CI.Tau
#' result_$CI.Par
#' }
#'
#' @export
#'
BiCopConfIntMMD <- function(
  x1, x2, family,
  nResampling = 100, subsamplingSize = length(x1),
  corrSubSampling = TRUE , level = 0.95, ...)
{

  # Preparation of the arguments
  arguments <- list(...)
  if (is.numeric(family)){
    estFUN = "BiCopEstMMD"
    arguments = c(list("family" = as.integer(family)) , arguments)
  } else if (family == "MO")
  {
    if ("method" %in% names(arguments)){
      method_ == arguments[["method"]]
      arguments[["method"]] <- NULL
    } else {
      # By default, we choose the MMD estimation method
      method_ = "MMD"
    }

    switch (method_,
            "curve" = {
              estFUN = "BiCopEst.MO.curve"
            },
            "itau" = {
              estFUN = "BiCopEst.MO.itau"
            },
            "MMD" = {
              estFUN = "BiCopEst.MO.MMD.MC"
              # Getting the kernel
              if ("kernel" %in% names(arguments)){
                kernel_ = arguments[["kernel"]]
                arguments[["kernel"]] <- NULL
              } else {
                kernel_ = "gaussian.Phi"
              }
              # Converting the name to a function
              if (is.character(kernel_))
              {
                kernelFun <- findKernelFunction(kernel_)
              } else {
                kernelFun <- kernel_
              }
              # Adding it back to the arguments list
              arguments = c(arguments, list(kernelFun = kernelFun))
            }
    )

  } else {
    stop("Unknown family ", family, " in BiCopConfIntMMD.")
  }

  n = length(x1)
  vecPar = rep(NA, nResampling)
  vecTau = rep(NA, nResampling)
  u1 = VineCopula::pobs(x1)
  u2 = VineCopula::pobs(x2)

  # Estimation using the whole sample
  estResult = do.call(estFUN, c(list(u1 = u1, u2 = u2), arguments))
  estPar = estResult$par
  estTau = estResult$tau

  pb = pbapply::startpb(min = 0, max = nResampling)
  if (subsamplingSize == n){

    # NP bootstrap
    for (iResampling in 1:nResampling){
      which_selected = sample.int(n = n, size = n, replace = TRUE)
      u1_st = VineCopula::pobs(x1[which_selected])
      u2_st = VineCopula::pobs(x2[which_selected])

      result = do.call(estFUN, c(list(u1 = u1_st, u2 = u2_st), arguments))
      vecTau[iResampling] = result$tau
      vecPar[iResampling] = result$par

      pbapply::setpb(pb, iResampling)
    }
    pbapply::closepb(pb)

    qLowPar = estPar - stats::quantile(
      vecPar - estPar,
      probs = 1 - (1-level)/2 )

    qHighPar = estPar - stats::quantile(
      vecPar - estPar,
      probs = (1-level)/2 )

    qLowTau = estTau - stats::quantile(
      vecTau - estTau,
      probs = 1 - (1-level)/2 )

    qHighTau = estTau - stats::quantile(
      vecTau - estTau,
      probs = (1-level)/2 )

  } else {

    # Subsampling
    for (iResampling in 1:nResampling){
      which_selected = sample.int(n = n, size = subsamplingSize, replace = TRUE)
      u1_st = VineCopula::pobs(x1[which_selected])
      u2_st = VineCopula::pobs(x2[which_selected])

      result = do.call(estFUN, c(list(u1 = u1_st, u2 = u2_st), arguments))
      vecTau[iResampling] = result$tau
      vecPar[iResampling] = result$par

      pbapply::setpb(pb, iResampling)
    }
    pbapply::closepb(pb)

    if (corrSubSampling){
      qLowPar = estPar - (1 - subsamplingSize / n)^{-1/2} *
        stats::quantile(vecPar - estPar, probs = 1 - (1-level)/2 )

      qHighPar = estPar - (1 - subsamplingSize / n)^{-1/2} *
        stats::quantile(vecPar - estPar, probs = (1-level)/2 )

      qLowTau = estTau - (1 - subsamplingSize / n)^{-1/2} *
        stats::quantile(vecTau - estTau, probs = 1 - (1-level)/2 )

      qHighTau = estTau - (1 - subsamplingSize / n)^{-1/2} *
        stats::quantile(vecTau - estTau, probs = (1-level)/2 )
    } else {
      qLowPar = estPar -
        stats::quantile(vecPar - estPar, probs = 1 - (1-level)/2 )

      qHighPar = estPar -
        stats::quantile(vecPar - estPar, probs = (1-level)/2 )

      qLowTau = estTau -
        stats::quantile(vecTau - estTau, probs = 1 - (1-level)/2 )

      qHighTau = estTau -
        stats::quantile(vecTau - estTau, probs = (1-level)/2 )
    }

  }
  CI.Par = c(qLowPar, qHighPar)
  CI.Par = c(qLowPar, qHighPar)
  names(CI.Par) <- paste(prettyNum(100 * c((1-level)/2 , (1 - (1-level)/2)  ), "%" ) )

  CI.Tau = c(qLowTau, qHighTau)
  CI.Tau = c(qLowTau, qHighTau)
  names(CI.Tau) <- paste(prettyNum(100 * c((1-level)/2 , (1 - (1-level)/2)  ), "%" ) )

  return(list(CI.Par = CI.Par, CI.Tau = CI.Tau
              # , vecPar = vecPar, vecTau = vecTau
              ))
}


