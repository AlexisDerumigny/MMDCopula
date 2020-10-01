## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MMDCopula)
library(VineCopula)
set.seed(1)

## -----------------------------------------------------------------------------

my_data = BiCopSim(N = 500, family = 1, par = 0.5)

estimator_MLE = BiCopEst(u1 = my_data[,1], u2 = my_data[,2], 
                         family = 1, method = "mle")

estimator_MMD = BiCopEstMMD(u1 = my_data[,1], u2 = my_data[,2], family = 1)

print(estimator_MLE)
print(estimator_MMD)


## -----------------------------------------------------------------------------

my_data_contam = my_data

number_outliers = 20
q = 0.001
my_data_contam[1:number_outliers, 1] = runif(n = number_outliers, min = 0, max = q)
my_data_contam[1:number_outliers, 2] = runif(n = number_outliers, min = 1-q, max = 1)

estimator_MLE = BiCopEst(u1 = my_data_contam[,1], u2 = my_data_contam[,2], 
                         family = 1, method = "mle")

estimator_MMD = BiCopEstMMD(u1 = my_data_contam[,1], u2 = my_data_contam[,2], family = 1)
print(estimator_MLE)
print(estimator_MMD)


