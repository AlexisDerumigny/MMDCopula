Package MMDCopula
===================


**Main functions**

* `BiCopEstMMD`: estimate the parameter of a parametric bivariate copula by MMD minimization.

* `BiCopConfIntMMD`: compute a bootstrap-based or subsampling-based confidence interval for the parameter of a parametric bivariate copula.

* `BiCopGradMMD`: compute the gradient of the MMD criteria. Used in `BiCopEstMMD`.

**Functions for simulation and inference for the Marshall-Olkin copula**

* `BiCopSim.MO`: simulation of observations following a Marshall-Olkin copula.

* `BiCopEst.MO`: estimation of the parameter of a Marshall-Olkin copula.

* `BiCopPar2Tau.MO` and `BiCopTau2Par.MO`: convert between the parameter and the Kendall's tau of a Marshall-Olkin copula.

**Other functions**

* `BiCopParamDistLp`: compute the $L^p$ distance between two parametric copula models.

