
# MMDCopula 0.2.0

* Change email address of maintainer Alexis Derumigny corresponding to new affiliation.
* Vectorized version of `BiCopSim.MO`.
* `BiCopSim.MO` now returns a list with the estimated parameter and Kendall's tau for consistency with the other `BiCopEstMMD` function.
* New function `BiCopMMDConfInt` for bootstrapped- and subsampling-based confidence intervals.
* Renaming of kernels such as `gaussian.KG` to `gaussian.Phi` for coherence.
* Improvement of the stochastic gradient interface.
* Better handling of limiting cases (when Kendall's tau is close to the boundary).

# MMDCopula 0.1.0

* Initial release

