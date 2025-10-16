# ============================================================
# 03.BC_functions.R
# ------------------------------------------------------------
# Essential functions for Simulation and Data Application under the
# Binary Outcome + Continuous Mediator setting.
#
# This script defines:
#   1. Data generation procedure (Surv.sim.bc)
#   2. True effect calculation (CMA.bc)
#   3. Bias-correction estimation (BCA.bc)
#
# Usage:
# - Used in both Simulation and Data Application analyses
#   when the outcome is binary and the mediator is continuous.
#
# Model setup:
#   X  ~ N(0, 0.5)
#   U  ~ U(0, 10)
#   ε_M ~ N(0, 1)   # error term for the mediator model
#   (Outcome Y is binary via a logistic model.)
#
# Dependencies:
# - Requires 00.basic_functions.R 
# ============================================================

source('./functions/00.basic_functions.R')

# ------------------------------------------------------------
# 1. Data Generation Function
# ------------------------------------------------------------
# Generate data for the Binary Outcome + Continuous Mediator setting.
# Inputs:
#   n      : sample size
#   beta   : parameter vector for the outcome (logistic) model
#   theta  : parameter vector for the mediator model
# Output:
#   list with simulated design matrix x = [X, M, U] and outcome y
# Notes:
#   - beta3: interaction effect between exposure X and mediator M
#   - beta4, th2: coefficients for baseline covariates U
# ------------------------------------------------------------
Surv.sim.bc <- function(n, beta, theta) {
  p <- length(beta); q <- length(theta)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- theta[1]; th1 <- theta[2]; th2 <- theta[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2   <- matrix(th2, nrow = 1)
  
  # Exposure variable X
  X <- rnorm(n, 0, sqrt(1 / 2))
  
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
  
  ## Mediator model (continuous mediator)
  M <- rnorm(n, th0 + th1 * X + th2 %*% t(U), sqrt(1))
  
  ## Outcome model (binary outcome via logistic link)
  Y <- rbinom(n, size = 1, prob = 1 / (1 + exp(-beta0 - beta1 * X - beta2 * M - beta3 * X * M - beta4 %*% t(U))))
  
  x <- cbind(X, M, U)
  list(x = x, y = Y)
}


# ------------------------------------------------------------
# 2. True Value Calculation Function
# ------------------------------------------------------------
# Compute true values of NDE, NIE, and MP based on model parameters.
# Inputs:
#   beta, th : true parameter vectors
#   sig2     : (variance of mediator error; used in NDE expression)
#   u        : covariate values
#   x1, x2   : exposure levels
# Output:
#   list(NIE, NDE, MP)
# Notes:
#   - Outcome is logistic; closed forms emerge in exponential (log-link) form.
# ------------------------------------------------------------
CMA.bc <- function(beta, th, sig2, u, x1, x2) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2   <- matrix(th2, nrow = 1)
  
  # Compute natural effects (exponential-form)
  NDE <- exp((beta1 + beta3 * (th0 + th1 * x2 + th2 %*% t(u) + beta2 * sig2)) * (x1 - x2) +
               (beta3^2 * sig2 * (x1^2 - x2^2)) / 2)
  NIE <- exp(th1 * (beta2 + beta3 * x1) * (x1 - x2))
  MP  <- log(NIE) / (log(NDE) + log(NIE))
  
  list(NIE = NIE, NDE = NDE, MP = MP)
}


# ------------------------------------------------------------
# 3. Bias-Correction Function (Binary Outcome + Continuous Mediator)
# ------------------------------------------------------------
# Apply exponential-form bias-correction (BC1 and BC2) to NDE and NIE.
# Apply ratio-form bias-correction (BC1 and BC2) to the mediation proportion (MP).
# Inputs:
#   beta, th : parameter vectors
#   u        : covariate data (coerced to numeric vector)
#   x1, x2   : exposure levels
#   sig2     : variance of mediator error
#   COV      : covariance matrix of estimators
# Output:
#   list of NDE, NIE, MP and their BC variants
# ------------------------------------------------------------
BCA.bc <- function(beta, th, u, x1, x2, sig2, COV) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2   <- matrix(th2, nrow = 1)
  u     <- unlist(u)  # Convert to numeric vector
  
  # --------------------------------------------------------
  # Matrix construction
  # --------------------------------------------------------
  len_c <- length(beta[5:p])
  
  # Non-symmetric matrices for quadratic forms
  E_D_bc <- t(matrix(c(
    0, 0, 1, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, sig2, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, sig2 * (x1 + x2) / 2, c(rep(0, len_c)), 1, x2, u,
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  E_I_bc <- t(matrix(c(
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 1, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, x1, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c)),
    rep(c(0, 0, 0, 0, 0, c(rep(0, len_c)), 0, 0, c(rep(0, len_c))), len_c)
  ), ncol = 1 + 6 + len_c * 2))
  
  # --------------------------------------------------------
  # Quadratic form setup
  # --------------------------------------------------------
  eta   <- c(beta, th)
  eta_1 <- c(1, eta)
  
  # Derivative matrix d(eta1)/d(eta)
  temp_1 <- c(0, 0, 0, 0, diag(len_c)[1, ], 0, 0, c(rep(0, len_c)))
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_1 <- c(0, 0, 0, 0, diag(len_c)[i, ], 0, 0, c(rep(0, len_c)))
      temp_1  <- c(temp_1, temp2_1)
    }
  }
  temp_2 <- c(0, 0, 0, 0, c(rep(0, len_c)), 0, 0, diag(len_c)[1, ])
  if (len_c > 1) {
    for (i in 2:len_c) {
      temp2_2 <- c(0, 0, 0, 0, c(rep(0, len_c)), 0, 0, diag(len_c)[i, ])
      temp_2  <- c(temp_2, temp2_2)
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
  h_bc_D <- t(eta_1) %*% E_D_bc %*% eta_1
  h_bc_I <- t(eta_1) %*% E_I_bc %*% eta_1
  
  # Derivatives of components (for BC)
  h_bc_D_df1 <- t(eta_1) %*% {E_D_bc + t(E_D_bc)} %*% d_eta1_eta
  h_bc_D_df2 <- t(d_eta1_eta) %*% {E_D_bc + t(E_D_bc)} %*% d_eta1_eta
  h_bc_I_df1 <- t(eta_1) %*% {E_I_bc + t(E_I_bc)} %*% d_eta1_eta
  h_bc_I_df2 <- t(d_eta1_eta) %*% {E_I_bc + t(E_I_bc)} %*% d_eta1_eta
  
  #---------------------
  # (1) NDE
  #---------------------
  # Ordinary NDE in exponential form
  NDE <- exp(h_bc_D * (x1 - x2))
  
  # Bias-corrected NDE (exponential-form BC1/BC2)
  NDE_bc1 <- exp(h_bc_D * (x1 - x2)) * {
    1 - 0.5 * (h_bc_D_df1 %*% COV %*% t(h_bc_D_df1)) - 0.5 * tr(h_bc_D_df2 %*% COV)
  }^(x1 - x2)
  NDE_bc2 <- exp(h_bc_D * (x1 - x2)) * {
    1 + 0.5 * (h_bc_D_df1 %*% COV %*% t(h_bc_D_df1)) + 0.5 * tr(h_bc_D_df2 %*% COV)
  }^(-(x1 - x2))
  
  #---------------------
  # (2) NIE
  #---------------------
  # Ordinary NIE in exponential form
  NIE <- exp(h_bc_I * (x1 - x2))
  
  # Bias-corrected NIE (exponential-form BC1/BC2)
  NIE_bc1 <- exp(h_bc_I * (x1 - x2)) * {
    1 - 0.5 * (h_bc_I_df1 %*% COV %*% t(h_bc_I_df1)) - 0.5 * tr(h_bc_I_df2 %*% COV)
  }^(x1 - x2)
  NIE_bc2 <- exp(h_bc_I * (x1 - x2)) * {
    1 + 0.5 * (h_bc_I_df1 %*% COV %*% t(h_bc_I_df1)) + 0.5 * tr(h_bc_I_df2 %*% COV)
  }^(-(x1 - x2))
  
  #---------------------
  # (3) MP = log(NIE) / (log(NIE) + log(NDE))
  #---------------------
  # MP inputs in log-scale (no exponentials inside the ratio)
  h1_bc <- t(eta_1) %*% E_I_bc %*% eta_1                # log(NIE)
  h2_bc <- t(eta_1) %*% {E_I_bc + E_D_bc} %*% eta_1     # log(NIE) + log(NDE)
  MP <- h1_bc / h2_bc
  
  # Derivatives for ratio-form BC
  h1_bc_df1 <- t(eta_1) %*% {E_I_bc + t(E_I_bc)} %*% d_eta1_eta
  h1_bc_df2 <- t(d_eta1_eta) %*% {E_I_bc + t(E_I_bc)} %*% d_eta1_eta
  h2_bc_df1 <- t(eta_1) %*% {E_I_bc + E_D_bc + t(E_I_bc + E_D_bc)} %*% d_eta1_eta
  h2_bc_df2 <- t(d_eta1_eta) %*% {E_I_bc + E_D_bc + t(E_I_bc + E_D_bc)} %*% d_eta1_eta
  
  # Ratio-form BC results for MP
  MP_bc1 <- ratio_BC1_fnt(
    COV,
    h1_bc, h2_bc,
    h1_bc_df1, h1_bc_df2,
    h2_bc_df1, h2_bc_df2
  )
  MP_bc2 <- ratio_BC2_fnt(
    COV,
    h1_bc, h2_bc,
    h1_bc_df1, h1_bc_df2,
    h2_bc_df1, h2_bc_df2
  )
  
  list(
    NIE = NIE, NDE = NDE,
    MP = MP,
    NIE_BC1 = NIE_bc1,
    NIE_BC2 = NIE_bc2,
    NDE_BC1 = NDE_bc1,
    NDE_BC2 = NDE_bc2,
    MP_BC1 = MP_bc1,
    MP_BC2 = MP_bc2
  )
}
