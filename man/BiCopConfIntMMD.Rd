% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiCopConfIntMMD.R
\name{BiCopConfIntMMD}
\alias{BiCopConfIntMMD}
\title{Confidence intervals for the estimated parameter
of a bivariate parametric copula using MMD estimation}
\usage{
BiCopConfIntMMD(
  x1,
  x2,
  family,
  nResampling = 100,
  subsamplingSize = length(x1),
  corrSubSampling = TRUE,
  level = 0.95,
  ...
)
}
\arguments{
\item{x1}{vector of observations of the first coordinate.}

\item{x2}{vector of observations of the second coordinate.}

\item{family}{parametric family of copulas.
Supported families are: \itemize{
  \item \code{1}: Gaussian copulas
  \item \code{3}: Clayton copulas
  \item \code{4}: Gumbel copulas
  \item \code{5}: Frank copulas
  \item \code{MO}: Marshall-Olkin copulas
  }}

\item{nResampling}{number of resampling times.}

\item{subsamplingSize}{size of the subsample.
By default it is \code{length(u1)},
i.e. this corresponds to the nonparametric boostrap.}

\item{corrSubSampling}{this parameter is only used for subsampling-based confidence intervals.
If \code{TRUE}, the confidence interval uses the corrected subsample empirical process.}

\item{level}{the nominal confidence level.}

\item{...}{other parameters to be given to \code{\link{BiCopEstMMD}}
or \code{\link{BiCopEst.MO}}.}
}
\value{
a list with the confidence intervals CI.Tau for Kendall's tau
and CI.Par for the corresponding parameter.
}
\description{
Confidence intervals for the estimated parameter
of a bivariate parametric copula using MMD estimation
}
\examples{
data = VineCopula::BiCopSim(N = 50, family = 1, par = 0.3)
result = BiCopConfIntMMD(x1 = data[,1], x2 = data[,2], family = 1,
  nResampling = 2, subsamplingSize = 10, niter = 10)
\donttest{
data_ = VineCopula::BiCopSim(N = 1000, family = 1, par = 0.3)
result_ = BiCopConfIntMMD(x1 = data_[,1], x2 = data_[,2], family = 1)
result_$CI.Tau
result_$CI.Par
}

}
\references{
Alquier, P., Chérief-Abdellatif, B.-E., Derumigny, A., and Fermanian, J.D. (2022).
Estimation of copulas via Maximum Mean Discrepancy.
Journal of the American Statistical Association, \doi{10.1080/01621459.2021.2024836}.

Kojadinovic I., and Stemikovskaya, K. (2019)
Subsampling (weighted smooth) empirical copula processes.
Journal of Multivariate Analysis, 173, 704-723,
\doi{10.1016/j.jmva.2019.05.007}.
}
