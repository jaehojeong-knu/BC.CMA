# ============================================================
# 01.CC_Simulation.R
# ------------------------------------------------------------
# Simulation script for Continuous Outcome and Continuous Mediator case
# corresponding to Section 4.1 in the manuscript.
#
# This script performs:
#   1) Parameter setting and true-value calculation
#   2) Monte Carlo simulation for multiple sample sizes
#   3) Estimation of NDE, NIE, MP and their bias-corrected versions
#   4) Summary and visualization-ready data generation
#
# Dependencies:
#   - Requires 00.basic_functions.R and 01.CC_functions.R
# ============================================================

source('./functions/00.basic_functions.R')
source('./functions/021.CC_functions.R')

# ------------------------------------------------------------
# 1) Parameter and True-Value Setting
# ------------------------------------------------------------
# Set exposure levels (x1, x2), covariate value u, and parameter sets.
u  <- 4
x1 <- 1
x2 <- 0

# Three parameter sets roughly targeting MP ≈ 50%, 30%, and 10%
# (Mediator model parameters: th*, Outcome model parameters: beta*)
th1   <- c(-3.00, 0.20, 0.50)
beta1 <- c(0.50, 2.00, 3.50, 1.00, 1.00)

th2   <- c(-3.00, 0.10, 1.00)
beta2 <- c(0.50, 0.50, 3.00, 0.50, 1.00)

th3   <- c(-3.00, 0.10, 1.15)
beta3 <- c(0.50, 0.20, 0.70, 0.50, 1.00)

theta_true <- cbind(th1, th2, th3)
beta_true  <- cbind(beta1, beta2, beta3)

# Compute "true" effects using CMA.cc (analytical target)
TrueNIE <- rep(0, 3)
TrueNDE <- rep(0, 3)
TrueMP  <- rep(0, 3)

for (m in 1:3) {
  i      <- m
  th.T   <- theta_true[, i]
  beta.T <- beta_true[, i]
  True   <- CMA.cc(beta.T, th.T, u = u, x1 = x1, x2 = x2)
  TrueNIE[i] <- True$NIE
  TrueNDE[i] <- True$NDE
  TrueMP[i]  <- True$MP
}
true_value_CC <- cbind(TrueNDE, TrueNIE, TrueMP, t(beta_true), t(theta_true)); true_value_CC



# ------------------------------------------------------------
# 2) Simulation Data Containers
# ------------------------------------------------------------
# Pre-allocate arrays for Monte Carlo results.
CASE.NIE    <- array()
CASE.NDE    <- array()
CASE.MP     <- array()
CASE.MP_BC1 <- array(); CASE.MP_BC2 <- array()


# ------------------------------------------------------------
# 2-1) Simulation Loop (Data Generation → Estimation → Store)
# ------------------------------------------------------------
# For each scenario i and sample size N, run M Monte Carlo iterations.
set.seed(2025)

for (i in c(1, 2, 3)) {
  for (N in c(100, 200, 300, 400, 600)) {
    theta_sim <- theta_true[, i]
    beta_sim  <- beta_true[, i]
    M <- 5000
    
    ############################# Simulation ##############################
    for (m in 1:M) {
      
      # --------------------------------------------------------
      # Data Generation
      # --------------------------------------------------------
      Dat <- Surv.sim.cc(n = N, beta = beta_sim, theta = theta_sim)
      Xc <- Dat$x[, 1]; Mc <- Dat$x[, 2]; Uc <- Dat$x[, 3:ncol(Dat$x)]
      Yc <- Dat$y
      data <- data.frame(Yc, Xc, Mc, Uc)
      
      # Build model formulas depending on the dimensionality of U
      if (is.null(colnames(Uc)) == TRUE) {
        med_formula <- as.formula("Mc ~ Xc + Uc")
        out_formula <- as.formula("Yc ~ Xc + Mc + Xc:Mc + Uc")
      } else if (is.null(colnames(Uc)) == FALSE) {
        med_formula <- as.formula(paste0("Mc ~ Xc + ", paste0(colnames(Uc), collapse = "+")))
        out_formula <- as.formula(paste0("Yc ~ Xc + Mc + Xc:Mc + ", paste0(colnames(Uc), collapse = "+")))
      }
      
      # --------------------------------------------------------
      # Mediation Model (Linear for Continuous Mediator)
      # --------------------------------------------------------
      Med <- lm(med_formula, data = data)
      # Skip this iteration if coefficients contain NA
      if (any(is.na(coef(Med)))) {
        cat(paste0("NA detected in mediation model at iteration ", m, ". Skipping...\n"))
        next
      }
      th        <- summary(Med)$coefficients[, 1]  # theta estimates
      sig_theta <- vcov(Med)                       # Var-Cov of theta
      sig2_med  <- summary(Med)$sigma^2            # Residual variance (mediator)
      
      # --------------------------------------------------------
      # Outcome Model (Linear for Continuous Outcome)
      # --------------------------------------------------------
      Outcome <- lm(out_formula, data = data)
      # Skip this iteration if coefficients contain NA
      if (any(is.na(coef(Outcome)))) {  
        cat(paste0("NA detected in outcome model at iteration ", m, ". Skipping...\n"))
        next
      }
      beta <- summary(Outcome)$coefficients[, 1]
      # Reorder to (Intercept, X, M, X:M, U's) as used in downstream functions
      beta <- c(beta[1], beta[2], beta[3], beta[length(beta)], beta[4:(length(beta) - 1)])
      sig_beta     <- vcov(Outcome)            # Var-Cov of beta
      sig2_outcome <- summary(Outcome)$sigma^2 # Residual variance
      
      # --------------------------------------------------------
      # Joint Covariance–Variance Matrix of (beta, theta)
      # --------------------------------------------------------
      COV <- VAR.C(beta = beta, th = th, u = u, sig.beta = sig_beta, sig.th = sig_theta)
      
      # --------------------------------------------------------
      # Estimation (Ordinary + Bias-Corrected)
      # --------------------------------------------------------
      CASE.CMA       <- BCA.cc(beta = beta, th = th, u = u, x1 = x1, x2 = x2, COV = COV)
      CASE.NIE[m]    <- CASE.CMA$NIE
      CASE.NDE[m]    <- CASE.CMA$NDE
      CASE.MP[m]     <- CASE.CMA$MP
      CASE.MP_BC1[m] <- CASE.CMA$MP_BC1
      CASE.MP_BC2[m] <- CASE.CMA$MP_BC2
      
      print(paste0("End ", i, " / ", N, " / ", m, " / ", M))
    }
    
    # --------------------------------------------------------
    # Save raw simulation outputs for this (i, N)
    # --------------------------------------------------------
    assign(
      paste0("LIST_CC_i", i, "_N", N),
      data.frame(
        NIE    = CASE.NIE,
        NDE    = CASE.NDE,
        MP     = CASE.MP,
        MP_BC1 = CASE.MP_BC1,
        MP_BC2 = CASE.MP_BC2
      )
    )
    
    write.csv(
      get(paste0("LIST_CC_i", i, "_N", N)),
      paste0(
        "./Sec4.Simulation/Simulation_data/01.CC_data/",
        paste0("LIST_CC_i", i, "_N", N), ".csv"
      ),
      row.names = FALSE
    )
  }
}


# ------------------------------------------------------------
# 3) Load Simulation Data & Prepare for Summary
# ------------------------------------------------------------
# Read back the saved CSV files for selected N and cases.
folder_path <- "./Sec4.Simulation/Simulation_data/01.CC_data/"

i_values     <- 1:3
N_values     <- c(200, 300, 600)
CC_data_list <- list()

# Load all (i, N) datasets
for (i in i_values) {
  for (N in N_values) {
    var_name  <- paste0("LIST_CC_i", i, "_N", N)
    file_path <- paste0(folder_path, var_name, ".csv")
    CC_data_list[[var_name]] <- read.csv(file_path)
  }
}


# ------------------------------------------------------------
# 4) Summary Calculations (Bias, Relative Bias, Variance)
# ------------------------------------------------------------
# Output tables for MP.
result_MP <- data.frame(
  Case = 0, N = 0, TrueMP = 0,
  Est = 0, Var = 0, Bias = 0, RBias = 0,
  Est1 = 0, Var1 = 0, Bias1 = 0, RBias1 = 0,
  Est2 = 0, Var2 = 0, Bias2 = 0, RBias2 = 0
)

# Wrapper to compute summary metrics for one dataset
# (i, N are parsed from the object/file names)
Simulation_summary_CC <- function(df, i, N) {
  
  #### MP  
  Est_MP        <- mean(df$MP, na.rm = TRUE)
  Var10_MP      <- 10 * var(df$MP, na.rm = TRUE)
  Bias_MP       <- Est_MP - true_value_CC[i, "TrueMP"]
  Rbias100_MP   <- 100 * (Est_MP - true_value_CC[i, "TrueMP"]) / true_value_CC[i, "TrueMP"]
  
  Est_MP_BC1    <- mean(df$MP_BC1, na.rm = TRUE)
  Var10_MP_BC1  <- 10 * var(df$MP_BC1, na.rm = TRUE)
  Bias_MP_BC1   <- Est_MP_BC1 - true_value_CC[i, "TrueMP"]
  Rbias100_MP_BC1 <- 100 * (Est_MP_BC1 - true_value_CC[i, "TrueMP"]) / true_value_CC[i, "TrueMP"]
  
  Est_MP_BC2    <- mean(df$MP_BC2, na.rm = TRUE)
  Var10_MP_BC2  <- 10 * var(df$MP_BC2, na.rm = TRUE)
  Bias_MP_BC2   <- Est_MP_BC2 - true_value_CC[i, "TrueMP"]
  Rbias100_MP_BC2 <- 100 * (Est_MP_BC2 - true_value_CC[i, "TrueMP"]) / true_value_CC[i, "TrueMP"]
  
  # Assemble summary rows
  result_MP <- data.frame(
    Case   = i, N = N,
    TrueMP = t(true_value_CC[i, c(3)]),
    Est    = round(Est_MP, 5),
    Var    = round(Var10_MP, 5),
    Bias   = round(Bias_MP, 5),
    RBias  = round(Rbias100_MP, 5),
    Est1   = round(Est_MP_BC1, 5),
    Var1   = round(Var10_MP_BC1, 5),
    Bias1  = round(Bias_MP_BC1, 5),
    RBias1 = round(Rbias100_MP_BC1, 5),
    Est2   = round(Est_MP_BC2, 5),
    Var2   = round(Var10_MP_BC2, 5),
    Bias2  = round(Bias_MP_BC2, 5),
    RBias2 = round(Rbias100_MP_BC2, 5)
  )
  return(list(result_MP))
}

# Run summaries across all loaded datasets
names(CC_data_list)
final_list <- Simulation_summary_CC(
  CC_data_list[[1]],
  i = as.numeric(str_extract(names(CC_data_list)[1], "(?<=_i)\\d+")),
  N = as.numeric(str_extract(names(CC_data_list)[1], "(?<=_N)\\d+"))
)
for (i in 1:length(CC_data_list)) {
  temp_list <- Simulation_summary_CC(
    CC_data_list[[i]],
    i = as.numeric(str_extract(names(CC_data_list)[i], "(?<=_i)\\d+")),
    N = as.numeric(str_extract(names(CC_data_list)[i], "(?<=_N)\\d+"))
  )
  final_list[[1]] <- rbind(final_list[[1]], temp_list[[1]])
}

# Reformat to publication-ready CSV
result_MP_CC <- make_csvfile(final_list[[1]])


# ------------------------------------------------------------
# 5) Outlier Inspection and Filtering
# ------------------------------------------------------------
# Identify extreme simulation values and remove them for stability. 
# Thresholds are manually set based on exploratory analysis.
# ------------------------------------------------------------
k <- 1
CC_data_list[[k]] <- CC_data_list[[k]] %>% dplyr::filter(abs(MP) < 1)

k <- 4
CC_data_list[[k]] <- CC_data_list[[k]] %>%
  dplyr::filter(MP < 1 & MP > -1.5)


# ------------------------------------------------------------
# 6) Final Summary and Export
# ------------------------------------------------------------
# Recompute summaries after filtering and export to CSV.
final_list <- Simulation_summary_CC(
  CC_data_list[[1]],
  i = as.numeric(str_extract(names(CC_data_list)[1], "(?<=_i)\\d+")),
  N = as.numeric(str_extract(names(CC_data_list)[1], "(?<=_N)\\d+"))
)

for (i in 1:length(CC_data_list)) {
  temp_list <- Simulation_summary_CC(
    CC_data_list[[i]],
    i = as.numeric(str_extract(names(CC_data_list)[i], "(?<=_i)\\d+")),
    N = as.numeric(str_extract(names(CC_data_list)[i], "(?<=_N)\\d+"))
  )
  final_list[[1]] <- rbind(final_list[[1]], temp_list[[1]])
}

result_MP_CC <- make_csvfile(final_list[[1]])

write.csv(
  result_MP_CC,
  "./Sec4_Simulation/Simulation_result/01.CC_result/Simul_result_CC_MP.csv",
  row.names = FALSE
)
