# ============================================================
# 04.BB_functions.R
# ------------------------------------------------------------
# Essential functions for Simulation and Data Application under the
# Binary Outcome + Binary Mediator setting.
#
# This script defines:
#   1. Data generation procedure (Surv.sim.bb)
#   2. True effect calculation (CMA.bb)
#   3. Bias-correction estimation (BCA.bb)
#
# Usage:
# - Used in both Simulation and Data Application analyses
#   when both the outcome and the mediator are binary.
#
# Model setup:
#   X  ~ N(0, 1)
#   U  ~ U(0, 10)
#   (Both M and Y are binary via logistic models.)
#
# Dependencies:
# - Requires 00.basic_functions.R 
# ============================================================

source('./functions/00.basic_functions.R')

# ------------------------------------------------------------
# 1. Data Generation Function
# ------------------------------------------------------------
# Generate data for the Binary Outcome + Binary Mediator setting.
# Inputs:
#   n      : sample size
#   beta   : parameter vector for the outcome (logistic) model
#   theta  : parameter vector for the mediator (logistic) model
# Output:
#   list with simulated design matrix x = [X, M, U] and outcome y
# Notes:
#   - beta3: interaction effect between exposure X and mediator M
#   - beta4, th2: coefficients for baseline covariates U
# ------------------------------------------------------------
Surv.sim.bb <- function(n, beta, theta) {
  p <- length(beta); q <- length(theta)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- theta[1]; th1 <- theta[2]; th2 <- theta[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2   <- matrix(th2, nrow = 1)
  
  # Exposure variable X
  X <- rnorm(n, 0, sqrt(1))
  
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
  
  ## Mediator model (binary mediator via logistic link)
  M <- rbinom(n, size = 1, prob = 1 / (1 + exp(-th0 - th1 * X - th2 %*% t(U))))
  
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
#   u        : covariate values
#   x1, x2   : exposure levels
# Output:
#   list(NIE, NDE, MP)
# Notes:
#   - Logistic models for both mediator and outcome yield closed-form expressions below.
# ------------------------------------------------------------
CMA.bb <- function(beta, th, u, x1, x2) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4, nrow = 1)
  th2   <- matrix(th2, nrow = 1)
  
  NDE <- exp(beta1 * (x1 - x2)) * {
    (1 + exp(beta2 + beta3 * x1 + th0 + th1 * x2 + th2 %*% t(u))) /
      (1 + exp(beta2 + beta3 * x2 + th0 + th1 * x2 + th2 %*% t(u)))
  }
  NIE <- {
    (1 + exp(th0 + th1 * x2 + th2 %*% t(u))) /
      (1 + exp(th0 + th1 * x1 + th2 %*% t(u)))
  } * {
    (1 + exp(beta2 + beta3 * x1 + th0 + th1 * x1 + th2 %*% t(u))) /
      (1 + exp(beta2 + beta3 * x1 + th0 + th1 * x2 + th2 %*% t(u)))
  }
  MP <- log(NIE) / (log(NDE) + log(NIE))
  
  list(NIE = NIE, NDE = NDE, MP = MP)
}


# ------------------------------------------------------------
# 3. Bias-Correction Function (Binary Outcome + Binary Mediator)
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
BCA.bb <- function(beta, th, u, x1, x2, COV) {
  p <- length(beta); q <- length(th)
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]; beta4 <- beta[5:p]
  th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
  
  beta4 <- matrix(beta4,nrow=1) 
  th2   <- matrix(th2,nrow=1) 
  u     <- unlist(u)  # Convert to numeric vector
  
  # --------------------------------------------------------
  # Matrix construction
  # --------------------------------------------------------
  len_c <- length(beta[5:p])
  
  # Non-symmetric matrices for quadratic forms
  E_D1_bb<-t(matrix(c(0,0,1,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  E_D2_bb<-t(matrix(c(0,0,0,1,x1,c(rep(0,len_c)),1,x2,u,
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  E_D3_bb<-t(matrix(c(0,0,0,1,x2,c(rep(0,len_c)),1,x2,u,
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  
  E_I1_bb<-t(matrix(c(0,0,0,0,0,c(rep(0,len_c)),1,x2,u,
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  E_I2_bb<-t(matrix(c(0,0,0,0,0,c(rep(0,len_c)),1,x1,u,
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  E_I3_bb<-t(matrix(c(0,0,0,1,x1,c(rep(0,len_c)),1,x1,u,
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                      rep(c(0,0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c))),len_c)),
                    ncol=1+6+len_c*2))
  
  # --------------------------------------------------------
  # Quadratic form setup
  # --------------------------------------------------------
  eta<-c(beta,th)
  eta_1<-c(1,eta)
  
  # Derivative matrix d(eta1)/d(eta)
  temp_1<-c(0,0,0,0,diag(len_c)[1,],0,0,c(rep(0,len_c)))
  if(len_c>1){
    for(i in 2:len_c){
      temp2_1<-c(0,0,0,0,diag(len_c)[i,],0,0,c(rep(0,len_c)))
      temp_1<-c(temp_1,temp2_1)
    }
  }
  temp_2<-c(0,0,0,0,c(rep(0,len_c)),0,0,diag(len_c)[1,])
  if(len_c>1){
    for(i in 2:len_c){
      temp2_2<-c(0,0,0,0,c(rep(0,len_c)),0,0,diag(len_c)[i,])
      temp_2<-c(temp_2,temp2_2)
    }
  }
  d_eta1_eta<-t(matrix(c(0,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                         1,0,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                         0,1,0,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                         0,0,1,0,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                         0,0,0,1,c(rep(0,len_c)),0,0,c(rep(0,len_c)),
                         temp_1,
                         0,0,0,0,c(rep(0,len_c)),1,0,c(rep(0,len_c)),
                         0,0,0,0,c(rep(0,len_c)),0,1,c(rep(0,len_c)),
                         temp_2),
                       ncol=1+6+len_c*2))
  
  # --------------------------------------------------------
  # NDE, NIE calculation
  # --------------------------------------------------------
  h_bb_D1 <- t(eta_1) %*% E_D1_bb %*% eta_1
  h_bb_D2 <- t(eta_1) %*% E_D2_bb %*% eta_1
  h_bb_D3 <- t(eta_1) %*% E_D3_bb %*% eta_1  
  h_bb_I1 <- t(eta_1) %*% E_I1_bb %*% eta_1  
  h_bb_I2 <- t(eta_1) %*% E_I2_bb %*% eta_1  
  h_bb_I3 <- t(eta_1) %*% E_I3_bb %*% eta_1
  
  # Derivatives of components (for BC)
  h_bb_D1_df1 <- t(eta_1) %*% {E_D1_bb+t(E_D1_bb)} %*% d_eta1_eta
  h_bb_D1_df2 <- t(d_eta1_eta) %*% {E_D1_bb+t(E_D1_bb)} %*% d_eta1_eta
  h_bb_D2_df1 <- t(eta_1) %*% {E_D2_bb+t(E_D2_bb)} %*% d_eta1_eta
  h_bb_D2_df2 <- t(d_eta1_eta) %*% {E_D2_bb+t(E_D2_bb)} %*% d_eta1_eta
  h_bb_D3_df1 <- t(eta_1) %*% {E_D3_bb+t(E_D3_bb)} %*% d_eta1_eta
  h_bb_D3_df2 <- t(d_eta1_eta) %*% {E_D3_bb+t(E_D3_bb)} %*% d_eta1_eta
  h_bb_I1_df1 <- t(eta_1) %*% {E_I1_bb+t(E_I1_bb)} %*% d_eta1_eta
  h_bb_I1_df2 <- t(d_eta1_eta) %*% {E_I1_bb+t(E_I1_bb)} %*% d_eta1_eta
  h_bb_I2_df1 <- t(eta_1) %*% {E_I2_bb+t(E_I2_bb)} %*% d_eta1_eta
  h_bb_I2_df2 <- t(d_eta1_eta) %*% {E_I2_bb+t(E_I2_bb)} %*% d_eta1_eta
  h_bb_I3_df1 <- t(eta_1) %*% {E_I3_bb+t(E_I3_bb)} %*% d_eta1_eta
  h_bb_I3_df2 <- t(d_eta1_eta) %*% {E_I3_bb+t(E_I3_bb)} %*% d_eta1_eta
  
  #---------------------
  # (1) NDE
  #---------------------
  # Ordinary NDE (product of an exponential term and a logistic ratio)
  exp_ratio_bbD2D3<-exponential_ratio_fnt_2(exp(h_bb_D2),exp(h_bb_D3))
  NDE <- exp(h_bb_D1*(x1-x2))*exp_ratio_bbD2D3;NDE
  
  # Bias correction for exp terms in NDE (bb)
  exp_bbD1_bc1_x1x2<-exp(h_bb_D1*(x1-x2))*{ 1 - 
      0.5*(h_bb_D1_df1 %*% COV %*% t(h_bb_D1_df1)) - 
      0.5*tr(h_bb_D1_df2 %*% COV)}^{x1-x2}
  exp_bbD1_bc2_x1x2<-exp(h_bb_D1*(x1-x2))*{ 1 + 
      0.5*(h_bb_D1_df1 %*% COV %*% t(h_bb_D1_df1)) + 
      0.5*tr(h_bb_D1_df2 %*% COV)}^{-(x1-x2)}
  exp_bbD2_bc1<-exponential_BC1_fnt(COV,
                                    h_bb_D2,h_bb_D2_df1,h_bb_D2_df2)
  exp_bbD2_bc2<-exponential_BC2_fnt(COV,
                                    h_bb_D2,h_bb_D2_df1,h_bb_D2_df2)
  exp_bbD3_bc1<-exponential_BC1_fnt(COV,
                                    h_bb_D3,h_bb_D3_df1,h_bb_D3_df2)
  exp_bbD3_bc2<-exponential_BC2_fnt(COV,
                                    h_bb_D3,h_bb_D3_df1,h_bb_D3_df2)
  
  # Bias-corrected NDE
  NDE_bc1 <- exp_bbD1_bc1_x1x2 * {(1+exp_bbD2_bc1)/(1+exp_bbD3_bc1)}
  NDE_bc2 <- exp_bbD1_bc2_x1x2 * {(1+exp_bbD2_bc2)/(1+exp_bbD3_bc2)}
  
  #---------------------
  # (2) NIE
  #---------------------
  # Ordinary NIE (product of two logistic ratios)
  exp_ratio_bbI1I2<-exponential_ratio_fnt_2(exp(h_bb_I1),exp(h_bb_I2))
  exp_ratio_bbI3D2<-exponential_ratio_fnt_2(exp(h_bb_I3),exp(h_bb_D2))
  NIE <- exp_ratio_bbI1I2*exp_ratio_bbI3D2;NIE
  
  # Bias correction for NIE components
  exp_bbI1_bc1<-exponential_BC1_fnt(COV,
                                    h_bb_I1,h_bb_I1_df1,h_bb_I1_df2)
  exp_bbI1_bc2<-exponential_BC2_fnt(COV,
                                    h_bb_I1,h_bb_I1_df1,h_bb_I1_df2)
  exp_bbI2_bc1<-exponential_BC1_fnt(COV,
                                    h_bb_I2,h_bb_I2_df1,h_bb_I2_df2)
  exp_bbI2_bc2<-exponential_BC2_fnt(COV,
                                    h_bb_I2,h_bb_I2_df1,h_bb_I2_df2)
  exp_bbI3_bc1<-exponential_BC1_fnt(COV,
                                    h_bb_I3,h_bb_I3_df1,h_bb_I3_df2)
  exp_bbI3_bc2<-exponential_BC2_fnt(COV,
                                    h_bb_I3,h_bb_I3_df1,h_bb_I3_df2)
  
  # Bias-corrected NIE
  NIE_bc1 <- {(1+exp_bbI1_bc1)/(1+exp_bbI2_bc1)} * {(1+exp_bbI3_bc1)/(1+exp_bbD2_bc1)}
  NIE_bc2 <- {(1+exp_bbI1_bc2)/(1+exp_bbI2_bc2)} * {(1+exp_bbI3_bc2)/(1+exp_bbD2_bc2)}
  
  #---------------------
  # (3) MP = NIE / (NDE + NIE)
  #---------------------
  # Helper ratios and derivatives for log(NDE), log(NIE)
  
  ## log(NDE) derivation
  exp_ratio_bbD1<-exponential_ratio_fnt(exp(h_bb_D1))
  exp_ratio_bbD2<-exponential_ratio_fnt(exp(h_bb_D2))
  exp_ratio_bbD3<-exponential_ratio_fnt(exp(h_bb_D3))
  
  exp_ratio_bbD2_bc1<-as.numeric(exp_bbD2_bc1/(1+exp_bbD2_bc1))
  exp_ratio_bbD2_bc2<-as.numeric(exp_bbD2_bc2/(1+exp_bbD2_bc2))
  exp_ratio_bbD3_bc1<-as.numeric(exp_bbD3_bc1/(1+exp_bbD3_bc1))
  exp_ratio_bbD3_bc2<-as.numeric(exp_bbD3_bc2/(1+exp_bbD3_bc2))
  
  log_NDE_bb_df1<-function(exp_ratio_bbD2,exp_ratio_bbD3){
    h_bb_D1_df1*(x1-x2) + h_bb_D2_df1*exp_ratio_bbD2 - h_bb_D3_df1*exp_ratio_bbD3} 
  log_NDE_bb_df2<-function(exp_ratio_bbD2,exp_ratio_bbD3){
    h_bb_D1_df2*(x1-x2) + h_bb_D2_df2*exp_ratio_bbD2 - h_bb_D3_df2*exp_ratio_bbD3 +
      #
      t(h_bb_D2_df1) %*% h_bb_D2_df1*exp_ratio_bbD2 - 
      t(h_bb_D2_df1) %*% h_bb_D2_df1*(exp_ratio_bbD2)^2 -
      #
      t(h_bb_D3_df1) %*% h_bb_D3_df1*exp_ratio_bbD3 + 
      t(h_bb_D3_df1) %*% h_bb_D3_df1*(exp_ratio_bbD3)^2}
  
  ## log(NIE) derivation
  exp_ratio_bbI1<-exponential_ratio_fnt(exp(h_bb_I1))
  exp_ratio_bbI2<-exponential_ratio_fnt(exp(h_bb_I2))
  exp_ratio_bbI3<-exponential_ratio_fnt(exp(h_bb_I3))
  
  exp_ratio_bbI1_bc1<-as.numeric(exp_bbI1_bc1/(1+exp_bbI1_bc1))
  exp_ratio_bbI1_bc2<-as.numeric(exp_bbI1_bc2/(1+exp_bbI1_bc2))
  exp_ratio_bbI2_bc1<-as.numeric(exp_bbI2_bc1/(1+exp_bbI2_bc1))
  exp_ratio_bbI2_bc2<-as.numeric(exp_bbI2_bc2/(1+exp_bbI2_bc2))
  exp_ratio_bbI3_bc1<-as.numeric(exp_bbI3_bc1/(1+exp_bbI3_bc1))
  exp_ratio_bbI3_bc2<-as.numeric(exp_bbI3_bc2/(1+exp_bbI3_bc2))
  
  log_NIE_bb_df1<-function(exp_ratio_bbI1,exp_ratio_bbI2,exp_ratio_bbI3){
    h_bb_I1_df1*exp_ratio_bbI1 - h_bb_I2_df1*exp_ratio_bbI2 +
      h_bb_I3_df1*exp_ratio_bbI3 - h_bb_D2_df1*exp_ratio_bbD2}
  log_NIE_bb_df2<-function(exp_ratio_bbI1,exp_ratio_bbI2,exp_ratio_bbI3){
    h_bb_I1_df2*exp_ratio_bbI1 - h_bb_I2_df2*exp_ratio_bbI2 + 
      h_bb_I3_df2*exp_ratio_bbI3 - h_bb_D2_df2*exp_ratio_bbD2 +
      #
      t(h_bb_I1_df1) %*% h_bb_I1_df1*exp_ratio_bbI1 - 
      t(h_bb_I1_df1) %*% h_bb_I1_df1*(exp_ratio_bbI1)^2 -
      #
      t(h_bb_I2_df1) %*% h_bb_I2_df1*exp_ratio_bbI2 + 
      t(h_bb_I2_df1) %*% h_bb_I2_df1*(exp_ratio_bbI2)^2 +
      #
      t(h_bb_I3_df1) %*% h_bb_I3_df1*exp_ratio_bbI3 - 
      t(h_bb_I3_df1) %*% h_bb_I3_df1*(exp_ratio_bbI3)^2 -
      #
      t(h_bb_D2_df1) %*% h_bb_D2_df1*exp_ratio_bbD2 + 
      t(h_bb_D2_df1) %*% h_bb_D2_df1*(exp_ratio_bbD2)^2}
  
  # MP and bias-corrected MP inputs
  h1_bb <- log(NIE)
  h2_bb <- log(NDE)+log(NIE)
  MP <- h1_bb/h2_bb ;MP
  
  h1_bb_bc1<-log(NIE_bc1)
  h1_bb_bc2<-log(NIE_bc2)
  h2_bb_bc1<-log(NDE_bc1)+log(NIE_bc1)
  h2_bb_bc2<-log(NDE_bc2)+log(NIE_bc2)
  
  # Derivatives evaluated at BC1/BC2-adjusted ratios
  h1_bb_df1_bc1<-log_NIE_bb_df1(exp_ratio_bbI1_bc1,exp_ratio_bbI2_bc1,exp_ratio_bbI3_bc1)
  h1_bb_df2_bc1<-log_NIE_bb_df2(exp_ratio_bbI1_bc1,exp_ratio_bbI2_bc1,exp_ratio_bbI3_bc1)
  h2_bb_df1_bc1<-log_NDE_bb_df1(exp_ratio_bbD2_bc1,exp_ratio_bbD3_bc1)+
    log_NIE_bb_df1(exp_ratio_bbI1_bc1,exp_ratio_bbI2_bc1,exp_ratio_bbI3_bc1)
  h2_bb_df2_bc1<-log_NDE_bb_df2(exp_ratio_bbD2_bc1,exp_ratio_bbD3_bc1)+
    log_NIE_bb_df2(exp_ratio_bbI1_bc1,exp_ratio_bbI2_bc1,exp_ratio_bbI3_bc1)
  
  h1_bb_df1_bc2<-log_NIE_bb_df1(exp_ratio_bbI1_bc2,exp_ratio_bbI2_bc2,exp_ratio_bbI3_bc2)
  h1_bb_df2_bc2<-log_NIE_bb_df2(exp_ratio_bbI1_bc2,exp_ratio_bbI2_bc2,exp_ratio_bbI3_bc2)
  h2_bb_df1_bc2<-log_NDE_bb_df1(exp_ratio_bbD2_bc2,exp_ratio_bbD3_bc2)+
    log_NIE_bb_df1(exp_ratio_bbI1_bc2,exp_ratio_bbI2_bc2,exp_ratio_bbI3_bc2)
  h2_bb_df2_bc2<-log_NDE_bb_df2(exp_ratio_bbD2_bc2,exp_ratio_bbD3_bc2)+
    log_NIE_bb_df2(exp_ratio_bbI1_bc2,exp_ratio_bbI2_bc2,exp_ratio_bbI3_bc2)
  
  # Ratio-form BC results for MP
  MP_bc1_insert_bc1<-ratio_BC1_fnt(COV,
                                   h1_bb_bc1,h2_bb_bc1, 
                                   h1_bb_df1_bc1,h1_bb_df2_bc1, 
                                   h2_bb_df1_bc1,h2_bb_df2_bc1)
  MP_bc1_insert_bc2<-ratio_BC1_fnt(COV,
                                   h1_bb_bc2,h2_bb_bc2, 
                                   h1_bb_df1_bc2,h1_bb_df2_bc2, 
                                   h2_bb_df1_bc2,h2_bb_df2_bc2)
  
  MP_bc2_insert_bc1<-ratio_BC2_fnt(COV,
                                   h1_bb_bc1,h2_bb_bc1, 
                                   h1_bb_df1_bc1,h1_bb_df2_bc1, 
                                   h2_bb_df1_bc1,h2_bb_df2_bc1)
  MP_bc2_insert_bc2<-ratio_BC2_fnt(COV,
                                   h1_bb_bc2,h2_bb_bc2, 
                                   h1_bb_df1_bc2,h1_bb_df2_bc2, 
                                   h2_bb_df1_bc2,h2_bb_df2_bc2)
  
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
