# ============================================================
# 01.CC_functions.R
# ------------------------------------------------------------
# Essential functions for Simulation and Data Application under the
# Continuous Outcome + Continuous Mediator setting.
#
# This script defines:
#   1. Data generation procedure (Surv.sim.cc)
#   2. True effect calculation (CMA.cc)
#   3. Bias-correction estimation (BCA.cc)
#
# Usage:
# - Used in both Simulation and Data Application analyses 
#   when both the outcome and mediator are continuous.
#
# Model setup:
#   X  ~ N(0, 1)
#   U  ~ U(0, 10)
#   ε_M ~ N(0, 2)   # error term for mediator model
#   ε_Y ~ N(0, 2)   # error term for outcome model
#
# Dependencies:
# - Requires 00.basic_functions.R
# ============================================================

source('./functions/00.basic_functions.R')

# ------------------------------------------------------------
# 1. Data Generation Function
# ------------------------------------------------------------
# Generate data for the continuous outcome & mediator setting.
# Inputs:
#   n      : sample size
#   beta   : parameter vector for outcome model
#   theta  : parameter vector for mediator model
# Output: 
#   list with simulated design matrix x = [X, M, U] and outcome y
# Notes:
#   - beta3: interaction effect between exposure X and mediator M
#   - beta4, th2: coefficients for baseline covariates U
# ------------------------------------------------------------
Surv.sim.cc <- function(n, beta, theta) {
  p <- length(beta); q <- length(theta)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- theta[1]; th1 <- theta[2]; th2 <- theta[3:q]
  
  beta4 <- matrix(beta4, nrow = 1) 
  th2 <- matrix(th2, nrow = 1)
  
  # Exposure variable X
  X <- rnorm(n, 0, 1)
  
  # Baseline covariates U
  if (p == 5) {
    U <- runif(n, 0, 1) * 10
  } else if (p > 5) {
    for (i in 1:(p - 4)) {
      assign(paste0("U", i), runif(n, 0, 1) * 10, envir = .GlobalEnv)
    }
    U_names <- paste0("U", 1:(p - 4))
    U <- do.call(cbind, lapply(U_names, function(name) get(name, envir = .GlobalEnv)))
    colnames(U) <- U_names
  }
  
  # Mediator model (continuous mediator)
  M <- rnorm(n, th0 + th1 * X + th2 %*% t(U), sqrt(2))
  
  # Outcome model (continuous outcome)
  Y <- rnorm(n, beta0 + beta1 * X + beta2 * M + beta3 * X * M + beta4 %*% t(U), sqrt(2))
  
  x <- cbind(X, M, U)
  list(x = x, y = Y)
}


# ------------------------------------------------------------
# 2. True Value Calculation Function
# ------------------------------------------------------------
# Compute true values of NDE, NIE, and MP based on model parameters.
# Inputs:
#   beta, th : true parameter vectors
#   u        : covariate values
#   x1, x2   : exposure levels
# Output:
#   list(NIE, NDE, MP)
# ------------------------------------------------------------
CMA.cc <- function(beta, th, u, x1, x2) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2 <- matrix(th2, nrow = 1)
  
  # Compute natural effects
  NDE <- (beta1 + beta3 * (th0 + th1 * x2 + th2 %*% t(u))) * (x1 - x2)
  NIE <- (beta2 * th1 + beta3 * th1 * x1) * (x1 - x2)
  MP <- NIE / (NDE + NIE)
  
  list(NIE = NIE, NDE = NDE, MP = MP)
}


# ------------------------------------------------------------
# 3. Bias-Correction Function for Continuous Outcome & Mediator
# ------------------------------------------------------------
# Apply ratio-form bias-correction (BC1 and BC2) to mediation proportion.
# Inputs:
#   beta, th : parameter vectors
#   u        : covariate data
#   x1, x2   : exposure levels
#   COV      : covariance matrix of estimators
# Output:
#   list of NDE, NIE, MP and their BC variants
# ------------------------------------------------------------
BCA.cc <- function(beta, th, u, x1, x2, COV) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2 <- matrix(th2, nrow = 1)
  u <- unlist(u)  # Convert to numeric vector
  
  # --------------------------------------------------------
  # Matrix construction
  # --------------------------------------------------------
  len_c <- length(beta[5:p])
  
  # Non-symmetric matrices for quadratic forms
  E_D_cc <- t(matrix(c(
    0, 0, 1, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    0, 0, 0, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    0, 0, 0, 0, 0, rep(0, len_c), 1, x2, u,
    rep(c(0, 0, 0, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c)), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_I_cc <- t(matrix(c(
    0, 0, 0, 0, 0, rep(0, len_c), 0, 1, rep(0, len_c),
    0, 0, 0, 0, 0, rep(0, len_c), 0, x1, rep(0, len_c),
    rep(c(0, 0, 0, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c)), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  # --------------------------------------------------------
  # Quadratic form setup
  # --------------------------------------------------------
  eta <- c(beta, th)
  eta_1 <- c(1, eta)
  
  # Derivative matrix d(eta1)/d(eta)
  temp_1 <- c(0, 0, 0, 0, diag(len_c)[1, ], 0, 0, rep(0, len_c))
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_1 <- c(0, 0, 0, 0, diag(len_c)[i, ], 0, 0, rep(0, len_c))
      temp_1 <- c(temp_1, temp2_1)
    }
  }
  
  temp_2 <- c(0, 0, 0, 0, rep(0, len_c), 0, 0, diag(len_c)[1, ])
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_2 <- c(0, 0, 0, 0, rep(0, len_c), 0, 0, diag(len_c)[i, ])
      temp_2 <- c(temp_2, temp2_2)
    }
  }
  
  d_eta1_eta <- t(matrix(c(
    0, 0, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    1, 0, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    0, 1, 0, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    0, 0, 1, 0, rep(0, len_c), 0, 0, rep(0, len_c),
    0, 0, 0, 1, rep(0, len_c), 0, 0, rep(0, len_c),
    temp_1,
    0, 0, 0, 0, rep(0, len_c), 1, 0, rep(0, len_c),
    0, 0, 0, 0, rep(0, len_c), 0, 1, rep(0, len_c),
    temp_2
  ), ncol = 1 + 6 + len_c * 2))
  
  # --------------------------------------------------------
  # NDE, NIE, MP calculation
  # --------------------------------------------------------
  NDE <- (t(eta_1) %*% E_D_cc %*% eta_1) * (x1 - x2)
  NIE <- (t(eta_1) %*% E_I_cc %*% eta_1) * (x1 - x2)
  h1_cc <- t(eta_1) %*% E_I_cc %*% eta_1
  h2_cc <- t(eta_1) %*% (E_I_cc + E_D_cc) %*% eta_1
  MP <- h1_cc / h2_cc
  
  # --------------------------------------------------------
  # Bias-correction derivation
  # --------------------------------------------------------
  h1_cc_df1 <- t(eta_1) %*% (E_I_cc + t(E_I_cc)) %*% d_eta1_eta
  h1_cc_df2 <- t(d_eta1_eta) %*% (E_I_cc + t(E_I_cc)) %*% d_eta1_eta
  h2_cc_df1 <- t(eta_1) %*% (E_I_cc + E_D_cc + t(E_I_cc + E_D_cc)) %*% d_eta1_eta
  h2_cc_df2 <- t(d_eta1_eta) %*% (E_I_cc + E_D_cc + t(E_I_cc + E_D_cc)) %*% d_eta1_eta
  
  MP_bc1 <- ratio_BC1_fnt(COV, h1_cc, h2_cc,
                          h1_cc_df1, h1_cc_df2,
                          h2_cc_df1, h2_cc_df2)
  
  MP_bc2 <- ratio_BC2_fnt(COV, h1_cc, h2_cc,
                          h1_cc_df1, h1_cc_df2,
                          h2_cc_df1, h2_cc_df2)
  
  list(
    NIE = NIE,
    NDE = NDE,
    MP = MP,
    MP_BC1 = MP_bc1,
    MP_BC2 = MP_bc2
  )
}
