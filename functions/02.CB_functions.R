# ============================================================
# 02.CB_functions.R
# ------------------------------------------------------------
# Essential functions for Simulation and Data Application under the
# Continuous Outcome + Binary Mediator setting.
#
# This script defines:
#   1. Data generation procedure (Surv.sim.cb)
#   2. True effect calculation (CMA.cb)
#   3. Bias-correction estimation (BCA.cb)
#
# Usage:
# - Used in both Simulation and Data Application analyses
#   when the outcome is continuous and the mediator is binary.
#
# Model setup:
#   X  ~ N(0, 3)
#   U  ~ U(0, 10)
#   ε_Y ~ N(0, 2)   # error term for the outcome model
#   (Mediator M is binary via a logistic model.)
#
# Dependencies:
# - Requires 00.basic_functions.R 
# ============================================================

source('./functions/00.basic_functions.R')

# ------------------------------------------------------------
# 1. Data Generation Function
# ------------------------------------------------------------
# Generate data for the Continuous Outcome + Binary Mediator setting.
# Inputs:
#   n      : sample size
#   beta   : parameter vector for the outcome model
#   theta  : parameter vector for the mediator (logistic) model
# Output:
#   list with simulated design matrix x = [X, M, U] and outcome y
# Notes:
#   - beta3: interaction effect between exposure X and mediator M
#   - beta4, th2: coefficients for baseline covariates U
# ------------------------------------------------------------
Surv.sim.cb <- function(n, beta, theta) {
  p <- length(beta); q <- length(theta)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- theta[1]; th1 <- theta[2]; th2 <- theta[3:q]
  
  beta4 <- matrix(beta4, nrow = 1) 
  th2 <- matrix(th2, nrow = 1) 
  
  # Exposure variable X
  X <- rnorm(n, 0, sqrt(3))
  
  # Baseline covariates U
  if (p == 5) {
    U <- runif(n, 0, 1.5) * 10
  } else if (p > 5) {
    for (i in 1:(p - 4)) {
      assign(paste0("U", i), runif(n, 0, 1) * 10, envir = .GlobalEnv)
    }
    U_names <- paste0("U", 1:(p - 4))
    U <- do.call(cbind, lapply(U_names, function(name) get(name, envir = .GlobalEnv))) 
    colnames(U) <- U_names
  }
  
  ## Mediation function (binary mediator via logistic link)
  M <- rbinom(n, size = 1, prob = 1 / (1 + exp(-th0 - th1 * X - th2 %*% t(U))))
  
  ## Outcome function (continuous outcome)
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
#   sig2     : (unused) kept for interface compatibility
#   u        : covariate values
#   x1, x2   : exposure levels
# Output:
#   list(NIE, NDE, MP)
# Notes:
#   - Uses the logistic mean for M|X,U in the decomposition.
#   - 'sig2' is not used in computations here (kept for consistency).
# ------------------------------------------------------------
CMA.cb <- function(beta, th, sig2, u, x1, x2) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1) 
  th2 <- matrix(th2, nrow = 1) 
  
  # Compute natural effects
  NDE <- (beta1 + beta3 * (1 / (1 + exp(-th0 - th1 * x2 - th2 %*% t(u))))) * (x1 - x2)
  NIE <- (beta2 + beta3 * x1) * ((1 / (1 + exp(-th0 - th1 * x1 - th2 %*% t(u)))) - (1 / (1 + exp(-th0 - th1 * x2 - th2 %*% t(u)))))
  MP <- NIE / (NDE + NIE)
  
  list(NIE = NIE, NDE = NDE, MP = MP)
}


# ------------------------------------------------------------
# 3. Bias-Correction Function (Continuous Outcome + Binary Mediator)
# ------------------------------------------------------------
# Apply exponential-form bias-correction (BC1 and BC2) to the natural direct effect (NDE) and natural indirect effect (NIE).
# Apply ratio-form bias-correction (BC1 and BC2) to the mediation proportion (MP).
# Inputs:
#   beta, th : parameter vectors
#   u        : covariate data (coerced to numeric vector)
#   x1, x2   : exposure levels
#   COV      : covariance matrix of estimators
# Output:
#   list of NDE, NIE, MP and their BC variants
# ------------------------------------------------------------
BCA.cb <- function(beta, th, u, x1, x2, COV) {
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
  E_D1_cb <- t(matrix(c(
    0, 0, 1, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_D2_cb <- t(matrix(c(
    0, 0, 0, 0, 1, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_D3_cb <- t(matrix(c(
    0, 0, 0, 0, 0, c(rep(0, len_c)), 1, x2, u,
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_I1_cb <- t(matrix(c(
    0, 0, 0, 1, x1, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_I2_cb <- t(matrix(c(
    0, 0, 0, 0, 0, c(rep(0, len_c)), 1, x1, u,
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_I3_cb <- t(matrix(c(
    0, 0, 0, 0, 0, c(rep(0, len_c)), 1, x2, u,
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  # --------------------------------------------------------
  # Quadratic form setup
  # --------------------------------------------------------
  eta <- c(beta, th)
  eta_1 <- c(1, eta)
  
  # Derivative matrix d(eta1)/d(eta)
  temp_1 <- c(0, 0, 0, 0, diag(len_c)[1, ], 0, 0, c(rep(0, len_c)))
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_1 <- c(0, 0, 0, 0, diag(len_c)[i, ], 0, 0, c(rep(0, len_c)))
      temp_1 <- c(temp_1, temp2_1)
    }
  }
  
  temp_2 <- c(0, 0, 0, 0, c(rep(0, len_c)), 0, 0, diag(len_c)[1, ])
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_2 <- c(0, 0, 0, 0, c(rep(0, len_c)), 0, 0, diag(len_c)[i, ])
      temp_2 <- c(temp_2, temp2_2)
    }
  }
  
  d_eta1_eta <- t(matrix(c(
    0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    1, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 1, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 1, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 1, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    temp_1,
    0, 0, 0, 0, c(rep(0, len_c)), 1, 0, c(rep(0, len_c)),
    0, 0, 0, 0, c(rep(0, len_c)), 0, 1, c(rep(0, len_c)),
    temp_2
  ), ncol = 1 + 6 + len_c * 2))
  
  # --------------------------------------------------------
  # NDE, NIE calculation
  # --------------------------------------------------------
  
  # NDE, NIE quadratic form values
  h_cb_D1 <- t(eta_1) %*% E_D1_cb %*% eta_1
  h_cb_D2 <- t(eta_1) %*% E_D2_cb %*% eta_1
  h_cb_D3 <- t(eta_1) %*% E_D3_cb %*% eta_1
  h_cb_I1 <- t(eta_1) %*% E_I1_cb %*% eta_1
  h_cb_I2 <- t(eta_1) %*% E_I2_cb %*% eta_1
  h_cb_I3 <- t(eta_1) %*% E_I3_cb %*% eta_1
  
  # Derivatives of components (for BC)
  h_cb_D1_df1 <- t(eta_1) %*% {E_D1_cb + t(E_D1_cb)} %*% d_eta1_eta
  h_cb_D1_df2 <- t(d_eta1_eta) %*% {E_D1_cb + t(E_D1_cb)} %*% d_eta1_eta
  h_cb_D2_df1 <- t(eta_1) %*% {E_D2_cb + t(E_D2_cb)} %*% d_eta1_eta
  h_cb_D2_df2 <- t(d_eta1_eta) %*% {E_D2_cb + t(E_D2_cb)} %*% d_eta1_eta
  h_cb_D3_df1 <- t(eta_1) %*% {E_D3_cb + t(E_D3_cb)} %*% d_eta1_eta
  h_cb_D3_df2 <- t(d_eta1_eta) %*% {E_D3_cb + t(E_D3_cb)} %*% d_eta1_eta
  
  h_cb_I1_df1 <- t(eta_1) %*% {E_I1_cb + t(E_I1_cb)} %*% d_eta1_eta
  h_cb_I1_df2 <- t(d_eta1_eta) %*% {E_I1_cb + t(E_I1_cb)} %*% d_eta1_eta
  h_cb_I2_df1 <- t(eta_1) %*% {E_I2_cb + t(E_I2_cb)} %*% d_eta1_eta
  h_cb_I2_df2 <- t(d_eta1_eta) %*% {E_I2_cb + t(E_I2_cb)} %*% d_eta1_eta
  h_cb_I3_df1 <- t(eta_1) %*% {E_I3_cb + t(E_I3_cb)} %*% d_eta1_eta
  h_cb_I3_df2 <- t(d_eta1_eta) %*% {E_I3_cb + t(E_I3_cb)} %*% d_eta1_eta
  
  #---------------------
  # (1) NDE
  #---------------------
  
  # Ordinary NDE via logistic mean (exp-ratio)
  exp_ratio_cbD3 <- exponential_ratio_fnt(exp(h_cb_D3))
  NDE <- { h_cb_D1 + h_cb_D2 * exp_ratio_cbD3 } * (x1 - x2); NDE
  
  # Bias correction for exp(h_cb_D3)
  exp_cbD3_bc1 <- exponential_BC1_fnt(COV, h_cb_D3, h_cb_D3_df1, h_cb_D3_df2)
  exp_cbD3_bc2 <- exponential_BC2_fnt(COV, h_cb_D3, h_cb_D3_df1, h_cb_D3_df2)
  exp_ratio_cbD3_bc1 <- as.numeric(exp_cbD3_bc1 / (1 + exp_cbD3_bc1))
  exp_ratio_cbD3_bc2 <- as.numeric(exp_cbD3_bc2 / (1 + exp_cbD3_bc2))
  
  # Bias-corrected NDE
  NDE_bc1 <- { h_cb_D1 + h_cb_D2 * exp_ratio_cbD3_bc1 } * (x1 - x2)
  NDE_bc2 <- { h_cb_D1 + h_cb_D2 * exp_ratio_cbD3_bc2 } * (x1 - x2)
  
  #---------------------
  # (2) NIE
  #---------------------
  
  # Ordinary NIE as a difference of logistic means
  exp_ratio_cbI2 <- exponential_ratio_fnt(exp(h_cb_I2))
  exp_ratio_cbI3 <- exponential_ratio_fnt(exp(h_cb_I3))
  NIE <- h_cb_I1 * { exp_ratio_cbI2 - exp_ratio_cbI3 }; NIE
  
  # Bias correction for exp(h_cb_I2) and exp(h_cb_I3)
  exp_cbI2_bc1 <- exponential_BC1_fnt(COV, h_cb_I2, h_cb_I2_df1, h_cb_I2_df2)
  exp_cbI2_bc2 <- exponential_BC2_fnt(COV, h_cb_I2, h_cb_I2_df1, h_cb_I2_df2)
  exp_cbI3_bc1 <- exponential_BC1_fnt(COV, h_cb_I3, h_cb_I3_df1, h_cb_I3_df2)
  exp_cbI3_bc2 <- exponential_BC2_fnt(COV, h_cb_I3, h_cb_I3_df1, h_cb_I3_df2)
  exp_ratio_cbI2_bc1 <- as.numeric(exp_cbI2_bc1 / (1 + exp_cbI2_bc1))
  exp_ratio_cbI2_bc2 <- as.numeric(exp_cbI2_bc2 / (1 + exp_cbI2_bc2))
  exp_ratio_cbI3_bc1 <- as.numeric(exp_cbI3_bc1 / (1 + exp_cbI3_bc1))
  exp_ratio_cbI3_bc2 <- as.numeric(exp_cbI3_bc2 / (1 + exp_cbI3_bc2))  
  
  # Bias-corrected NIE
  NIE_bc1 <- h_cb_I1 * { exp_ratio_cbI2_bc1 - exp_ratio_cbI3_bc1 }
  NIE_bc2 <- h_cb_I1 * { exp_ratio_cbI2_bc2 - exp_ratio_cbI3_bc2 }
  
  #---------------------
  # (3) MP = NIE / (NDE + NIE)
  #---------------------
  
  # Helper derivatives for ratio-form BC
  NDE_cb_df1 <- function(exp_ratio_cbD3) {
    h_cb_D1_df1 + 
      h_cb_D2_df1 * exp_ratio_cbD3 +
      h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3) -
      h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3^2)
  } 
  
  NDE_cb_df2 <- function(exp_ratio_cbD3) {
    h_cb_D1_df2 + 
      #
      h_cb_D2_df2 * exp_ratio_cbD3 + 
      t(h_cb_D2_df1) %*% h_cb_D3_df1 * exp_ratio_cbD3 - 
      t(h_cb_D2_df1) %*% h_cb_D3_df1 * exp_ratio_cbD3^2 +
      #
      h_cb_D3_df2 * as.numeric(h_cb_D2 * exp_ratio_cbD3) + 
      t(h_cb_D3_df1) %*% h_cb_D2_df1 * exp_ratio_cbD3 + 
      t(h_cb_D3_df1) %*% h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3) - 
      t(h_cb_D3_df1) %*% h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3^2) -
      #
      h_cb_D3_df2 * as.numeric(h_cb_D2 * exp_ratio_cbD3^2) -
      t(h_cb_D3_df1) %*% h_cb_D2_df1 * exp_ratio_cbD3^2 - 
      2 * t(h_cb_D3_df1) %*% h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3^2) +
      2 * t(h_cb_D3_df1) %*% h_cb_D3_df1 * as.numeric(h_cb_D2 * exp_ratio_cbD3^3)
  }
  
  NIE_cb_df1 <- function(exp_ratio_cbI2, exp_ratio_cbI3) {
    h_cb_I1_df1 * exp_ratio_cbI2 - h_cb_I1_df1 * exp_ratio_cbI3 +
      #  
      h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2) -
      h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2^2) -
      #  
      h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3) +
      h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3^2)
  }
  
  NIE_cb_df2 <- function(exp_ratio_cbI2, exp_ratio_cbI3) {
    (h_cb_I1_df2 * exp_ratio_cbI2 + t(h_cb_I1_df1) %*% h_cb_I2_df1 * exp_ratio_cbI2 - t(h_cb_I1_df1) %*% h_cb_I2_df1 * exp_ratio_cbI2^2) -
      (h_cb_I1_df2 * exp_ratio_cbI3 + t(h_cb_I1_df1) %*% h_cb_I3_df1 * exp_ratio_cbI3 - t(h_cb_I1_df1) %*% h_cb_I3_df1 * exp_ratio_cbI3^2) +
      #
      (h_cb_I2_df2 * as.numeric(h_cb_I1 * exp_ratio_cbI2) + 
         t(h_cb_I2_df1) %*% h_cb_I1_df1 * exp_ratio_cbI2 + 
         t(h_cb_I2_df1) %*% h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2) - 
         t(h_cb_I2_df1) %*% h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2^2) -
         t(h_cb_I2_df2) * as.numeric(h_cb_I1 * exp_ratio_cbI2^2) -
         t(h_cb_I2_df1) %*% h_cb_I1_df1 * exp_ratio_cbI2^2 -
         2 * t(h_cb_I2_df1) %*% h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2^2) +
         2 * t(h_cb_I2_df1) %*% h_cb_I2_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI2^3)) -
      #
      (h_cb_I3_df2 * as.numeric(h_cb_I1 * exp_ratio_cbI3) + 
         t(h_cb_I3_df1) %*% h_cb_I1_df1 * exp_ratio_cbI3 +
         t(h_cb_I3_df1) %*% h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3) -
         t(h_cb_I3_df1) %*% h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3^2) -
         t(h_cb_I3_df2) * as.numeric(h_cb_I1 * exp_ratio_cbI3^2) -
         t(h_cb_I3_df1) %*% h_cb_I1_df1 * exp_ratio_cbI3^2 -
         2 * t(h_cb_I3_df1) %*% h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3^2) +
         2 * t(h_cb_I3_df1) %*% h_cb_I3_df1 * as.numeric(h_cb_I1 * exp_ratio_cbI3^3))
  }
  
  # MP and bias-corrected MP inputs
  h1_cb <- NIE
  h2_cb <- NDE + NIE
  MP <- h1_cb / h2_cb
  
  h1_cb_bc1 <- NIE_bc1
  h1_cb_bc2 <- NIE_bc2
  h2_cb_bc1 <- NDE_bc1 + NIE_bc1
  h2_cb_bc2 <- NDE_bc2 + NIE_bc2
  
  # Derivatives evaluated at BC1/BC2-adjusted exp-ratios
  h1_cb_df1_bc1 <- NIE_cb_df1(exp_ratio_cbI2_bc1, exp_ratio_cbI3_bc1)
  h1_cb_df2_bc1 <- NIE_cb_df2(exp_ratio_cbI2_bc1, exp_ratio_cbI3_bc1)
  h2_cb_df1_bc1 <- NDE_cb_df1(exp_ratio_cbD3_bc1) + NIE_cb_df1(exp_ratio_cbI2_bc1, exp_ratio_cbI3_bc1)
  h2_cb_df2_bc1 <- NDE_cb_df2(exp_ratio_cbD3_bc1) + NIE_cb_df2(exp_ratio_cbI2_bc1, exp_ratio_cbI3_bc1)
  
  h1_cb_df1_bc2 <- NIE_cb_df1(exp_ratio_cbI2_bc2, exp_ratio_cbI3_bc2)
  h1_cb_df2_bc2 <- NIE_cb_df2(exp_ratio_cbI2_bc2, exp_ratio_cbI3_bc2)
  h2_cb_df1_bc2 <- NDE_cb_df1(exp_ratio_cbD3_bc2) + NIE_cb_df1(exp_ratio_cbI2_bc2, exp_ratio_cbI3_bc2)
  h2_cb_df2_bc2 <- NDE_cb_df2(exp_ratio_cbD3_bc2) + NIE_cb_df2(exp_ratio_cbI2_bc2, exp_ratio_cbI3_bc2)
  
  # Ratio-form BC results for MP
  MP_bc1_insert_bc1 <- ratio_BC1_fnt(
    COV,
    h1_cb_bc1, h2_cb_bc1, 
    h1_cb_df1_bc1, h1_cb_df2_bc1, 
    h2_cb_df1_bc1, h2_cb_df2_bc1
  )
  MP_bc1_insert_bc2 <- ratio_BC1_fnt(
    COV,
    h1_cb_bc2, h2_cb_bc2, 
    h1_cb_df1_bc2, h1_cb_df2_bc2, 
    h2_cb_df1_bc2, h2_cb_df2_bc2
  )
  
  MP_bc2_insert_bc1 <- ratio_BC2_fnt(
    COV,
    h1_cb_bc1, h2_cb_bc1, 
    h1_cb_df1_bc1, h1_cb_df2_bc1, 
    h2_cb_df1_bc1, h2_cb_df2_bc1
  )
  MP_bc2_insert_bc2 <- ratio_BC2_fnt(
    COV,
    h1_cb_bc2, h2_cb_bc2, 
    h1_cb_df1_bc2, h1_cb_df2_bc2, 
    h2_cb_df1_bc2, h2_cb_df2_bc2
  )
  
  list(
    NIE = NIE, NDE = NDE,
    MP = MP,
    NIE_BC1 = NIE_bc1,
    NIE_BC2 = NIE_bc2,
    NDE_BC1 = NDE_bc1,
    NDE_BC2 = NDE_bc2,
    MP_BC1_insert_bc1 = MP_bc1_insert_bc1,
    MP_BC1_insert_bc2 = MP_bc1_insert_bc2,
    MP_BC2_insert_bc1 = MP_bc2_insert_bc1,
    MP_BC2_insert_bc2 = MP_bc2_insert_bc2
  )
}
