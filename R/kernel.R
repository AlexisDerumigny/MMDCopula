

# Takes in parameter the name of the kernel
# and return the appropriate kernel function.
findKernelFunction = function(kernel) {

  switch ( kernel,

    "gaussian" = { kernelFun <- KerMMD.gaussian },

    "gaussian.Phi" = { kernelFun <- KerMMD.gaussian.Phi },

    "exp-l2" = { kernelFun <- KerMMD.exp.l2 },

    "exp-l2.Phi" = { kernelFun <- KerMMD.exp.l2.Phi },

    "exp-l1" = { kernelFun <- KerMMD.exp.l1 },

    "exp-l1.Phi" = { kernelFun <- KerMMD.exp.l1.Phi },

    "inv-l2" = { kernelFun <- KerMMD.inv.l2 },

    "inv-l2.Phi" = { kernelFun <- KerMMD.inv.l2.Phi },

    "inv-l1" = { kernelFun <- KerMMD.inv.l1 },

    "inv-l1.Phi" = { kernelFun <- KerMMD.inv.l1.Phi },

    {
      stop("Unknown kernel: ", kernel, ". ",
           "Possible values are: 'gaussian', 'gaussian.Phi', ",
           "'exp-l2', 'exp-l2.Phi', 'exp-l1', 'exp-l1.Phi', ",
           "'inv-l2', 'inv-l2.Phi', 'inv-l1' and 'inv-l1.Phi'.")
    }
  )

  return (kernelFun)
}



# Gaussian ---------------------------------------------------

KerMMD.gaussian.Phi = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (stats::qnorm(u1) - stats::qnorm(v1))/gamma
  D2 = (stats::qnorm(u2) - stats::qnorm(v2))/gamma

  Ker <- exp(-D1^2-D2^2)
  return(Ker)
}

KerMMD.gaussian = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (u1-v1)/gamma
  D2 = (u2-v2)/gamma

  Ker <- exp(-D1^2-D2^2)
  return(Ker)
}


# Exp.l2 ------------------------------------------------------

KerMMD.exp.l2.Phi = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (stats::qnorm(u1)-stats::qnorm(v1))/gamma
  D2 = (stats::qnorm(u2)-stats::qnorm(v2))/gamma

  Ker <- exp(-sqrt(D1^2+D2^2))
  return(Ker)
}

KerMMD.exp.l2 = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (u1-v1)/gamma
  D2 = (u2-v2)/gamma

  Ker <- exp(-sqrt(D1^2+D2^2))
  return(Ker)
}


# Exp.l1 ------------------------------------------------------

KerMMD.exp.l1.Phi = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = abs(stats::qnorm(u1)-stats::qnorm(v1))/gamma
  D2 = abs(stats::qnorm(u2)-stats::qnorm(v2))/gamma

  Ker <- exp(-D1-D2)
  return(Ker)
}

KerMMD.exp.l1 = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = abs(u1-v1)/gamma
  D2 = abs(u2-v2)/gamma

  Ker <- exp(-D1-D2)
  return(Ker)
}


# Inv.l2 ------------------------------------------------------

KerMMD.inv.l2.Phi = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (stats::qnorm(u1)-stats::qnorm(v1))/gamma
  D2 = (stats::qnorm(u2)-stats::qnorm(v2))/gamma

  Ker <- 1/(1+sqrt(D1^2+D2^2))^alpha
  return(Ker)
}

KerMMD.inv.l2 = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = (u1-v1)/gamma
  D2 = (u2-v2)/gamma

  Ker <- 1/(1+sqrt(D1^2+D2^2))^alpha
  return(Ker)
}


# Inv.l1 ------------------------------------------------------

KerMMD.inv.l1.Phi = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = abs(stats::qnorm(u1)-stats::qnorm(v1))/gamma
  D2 = abs(stats::qnorm(u2)-stats::qnorm(v2))/gamma

  Ker <- 1/(1+D1+D2)^alpha
  return(Ker)
}

KerMMD.inv.l1 = function(u1, u2, v1, v2, gamma=0.3, alpha=1) {
  D1 = abs(u1-v1)/gamma
  D2 = abs(u2-v2)/gamma

  Ker <- 1/(1+D1+D2)^alpha
  return(Ker)
}

