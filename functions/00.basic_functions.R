# ============================================================
# 00.basic_functions.R
# ------------------------------------------------------------
# This script defines utility functions commonly used in both 
# Simulation and Data Application analyses. 
# It provides bias-correction functions, covariance construction, 
# and data pre-processing utilities for result export.
#
# Note:
# - This file should be sourced at the beginning of all R scripts 
#   within the Simulation and Data Application folders.
# ============================================================

# ============================================================
# Libraries
# ============================================================
# suppressMessages({
#   library(dplyr)
#   library(stringr)
#   library(data.table)
#   library(ggplot2)
#   library(psych)
#   library(medflex)
#   library(mediation)
#   library(boot)
# })

# ------------------------------------------------------------
# Miscellaneous Utilities
# ------------------------------------------------------------

# Custom "not in" operator
`%!in%` <- function(x, y) !(x %in% y)


# ------------------------------------------------------------
# Covariance–Variance Matrix for Linear Regression
# ------------------------------------------------------------
VAR.C <- function(beta, th, u, sig.beta, sig.th) {
  p <- length(beta); q <- length(th); k <- length(u)
  COV <- matrix(0, nrow = p + q, ncol = p + q)
  
  COV[1:p, 1:p] <- sig.beta[c(1:(p - k - 1), p:(p - k)), c(1:(p - k - 1), p:(p - k))]
  COV[(p + 1):(p + q), (p + 1):(p + q)] <- sig.th
  
  COV
}


# ------------------------------------------------------------
# Bias-Correction Functions: Exponential Form
# ------------------------------------------------------------

# Ratio of exponential form
exponential_ratio_fnt <- function(exp_h_fnt) {
  as.numeric(exp_h_fnt / (1 + exp_h_fnt))
}

# Exponential Bias Correction (BC1)
exponential_BC1_fnt <- function(COV, h_fnt, h_fnt_df1, h_fnt_df2) {
  as.numeric(exp(h_fnt) * {
    1 - 0.5 * (h_fnt_df1 %*% COV %*% t(h_fnt_df1)) -
      0.5 * tr(h_fnt_df2 %*% COV)
  })
}

# Exponential Bias Correction (BC2)
exponential_BC2_fnt <- function(COV, h_fnt, h_fnt_df1, h_fnt_df2) {
  as.numeric(exp(h_fnt) * {
    1 + 0.5 * (h_fnt_df1 %*% COV %*% t(h_fnt_df1)) +
      0.5 * tr(h_fnt_df2 %*% COV)
  }^(-1))
}


# ------------------------------------------------------------
# Bias-Correction Functions: Ratio of Exponential Forms
# ------------------------------------------------------------

# Alternative ratio of two exponential terms
exponential_ratio_fnt_2 <- function(exp_h_fnt1, exp_h_fnt2) {
  as.numeric((1 + exp_h_fnt1) / (1 + exp_h_fnt2))
}


# ------------------------------------------------------------
# Bias-Correction Functions: Ratio Form
# ------------------------------------------------------------

# Ratio Bias Correction (BC1)
ratio_BC1_fnt <- function(COV,
                          h1, h2,
                          h1_df1, h1_df2,
                          h2_df1, h2_df2) {
  MP <- h1 / h2
  as.numeric(MP * {
    1 -
      tr(h1_df2 %*% COV) / (2 * h1) +
      (h1_df1 %*% COV %*% t(h2_df1)) / (h1 * h2) +
      tr(h2_df2 %*% COV) / (2 * h2) -
      (h2_df1 %*% COV %*% t(h2_df1)) / (h2^2)
  })
}

# Ratio Bias Correction (BC2)
ratio_BC2_fnt <- function(COV,
                          h1, h2,
                          h1_df1, h1_df2,
                          h2_df1, h2_df2) {
  MP <- h1 / h2
  as.numeric(MP * {
    (h1_df1 %*% COV %*% t(h2_df1)) / (h1 * h2) +
      tr(h2_df2 %*% COV) / (2 * h2) +
      {1 + tr(h1_df2 %*% COV) / (2 * h1) +
          (h2_df1 %*% COV %*% t(h2_df1)) / (h2^2)}^(-1)
  })
}


# ------------------------------------------------------------
# File Preprocessing Utilities
# ------------------------------------------------------------

# Function to reshape and prepare simulation results for CSV export
make_csvfile <- function(data) {
  df <- cbind(
    reshape2::melt(data[-1, c("Case", "N", "Est", "Est1", "Est2")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "Var", "Var1", "Var2")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "Bias", "Bias1", "Bias2")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "RBias", "RBias1", "RBias2")], id.vars = c("Case", "N")) %>%
      arrange(Case, N)
  ) %>% as.data.frame()
  
  df2 <- df[, c(1, 2, 3, 4, 8, 12, 16)]
  colnames(df2) <- c("Case", "N", "Method", "Est", "Var", "Bias", "RBias")
  
  return(df2)
}

# Function to reshape and prepare MP (mediation proportion) results for CSV export
make_csvfile_MP <- function(data) {
  df <- cbind(
    reshape2::melt(data[-1, c("Case", "N", "Est", "Est11", "Est12", "Est21", "Est22")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "Var", "Var11", "Var12", "Var21", "Var22")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "Bias", "Bias11", "Bias12", "Bias21", "Bias22")], id.vars = c("Case", "N")) %>%
      arrange(Case, N),
    reshape2::melt(data[-1, c("Case", "N", "RBias", "RBias11", "RBias12", "RBias21", "RBias22")], id.vars = c("Case", "N")) %>%
      arrange(Case, N)
  ) %>% as.data.frame()
  
  df2 <- df[, c(1, 2, 3, 4, 8, 12, 16)]
  colnames(df2) <- c("Case", "N", "Method", "Est", "Var", "Bias", "RBias")
  
  return(df2)
}

