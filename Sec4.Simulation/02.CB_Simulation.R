# ============================================================
# 02.CB_Simulation.R
# ------------------------------------------------------------
# Simulation script for Continuous Outcome and Binary Mediator case
# corresponding to Section 4.2 in the manuscript.
#
# This script performs:
#   1) Parameter setting and true-value calculation
#   2) Monte Carlo simulation for multiple sample sizes
#   3) Estimation of NDE, NIE, MP and their bias-corrected versions
#   4) Summary and visualization-ready data generation
#
# Dependencies:
#   - Requires 00.basic_functions.R and 02.CB_functions.R
# ============================================================

source('./functions/00.basic_functions.R')
source('./functions/02.CB_functions.R')

# ------------------------------------------------------------
# 1) Parameter and True-Value Setting
# ------------------------------------------------------------
# Set exposure levels (x1, x2), covariate value u, and parameter sets.
u  <- 4
x1 <- 1
x2 <- 0

# Three parameter sets roughly targeting MP ≈ 50%, 30%, and 10%
# (Mediator model parameters: th*, Outcome model parameters: beta*)
th1   <- c(-3.50, 2.20, 1.10)
beta1 <- c(0.50, 0.10, 2.00, 1.00, 1.00)

th2   <- c(-3.00, 2.20, 1.10)
beta2 <- c(0.50, 0.20, 1.50, 1.00, 1.00)

th3   <- c(-2.20, 2.20, 1.10)
beta3 <- c(0.50, 0.20, 0.60, 2.20, 1.00)

theta_true <- cbind(th1, th2, th3)
beta_true  <- cbind(beta1, beta2, beta3)

# Compute "true" effects using CMA.cb (analytical target)
TrueNIE <- rep(0, 3)
TrueNDE <- rep(0, 3)
TrueMP  <- rep(0, 3)

for (m in 1:3) {
  i      <- m
  th.T   <- theta_true[, i]
  beta.T <- beta_true[, i]
  True   <- CMA.cb(beta.T, th.T, sig2 = 1, u, x1 = x1, x2 = x2)
  TrueNIE[i] <- True$NIE
  TrueNDE[i] <- True$NDE
  TrueMP[i]  <- True$MP
}
true_value_CB <- cbind(TrueNDE, TrueNIE, TrueMP, t(beta_true), t(theta_true)); true_value_CB


# ------------------------------------------------------------
# 2) Simulation Data Containers
# ------------------------------------------------------------
# Pre-allocate arrays for Monte Carlo results.
CASE.NIE     <- array()
CASE.NIE_BC1 <- array(); CASE.NIE_BC2 <- array()
CASE.NDE     <- array()
CASE.NDE_BC1 <- array(); CASE.NDE_BC2 <- array()
CASE.MP      <- array()
CASE.MP_BC11 <- array(); CASE.MP_BC12 <- array()
CASE.MP_BC21 <- array(); CASE.MP_BC22 <- array()


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
      Dat <- Surv.sim.cb(n = N, beta = beta_sim, theta = theta_sim)
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
      # Mediation Model (Logistic for Binary Mediator)
      # --------------------------------------------------------
      Med <- glm(med_formula, family = binomial, data = data)
      # Skip this iteration if coefficients contain NA
      if (any(is.na(coef(Med)))) {
        cat(paste0("NA detected in mediation model at iteration ", m, ". Skipping...\n"))
        next
      }
      th        <- summary(Med)$coefficients[, 1]  # theta estimates
      sig_theta <- vcov(Med)                       # Var-Cov of theta
      
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
      sig2_outcome <- summary(Outcome)$sigma^2 # Residual variance (not directly used here)
      
      # --------------------------------------------------------
      # Joint Covariance–Variance Matrix of (beta, theta)
      # --------------------------------------------------------
      COV <- VAR.C(beta = beta, th = th, u = u, sig.beta = sig_beta, sig.th = sig_theta)
      
      # --------------------------------------------------------
      # Estimation (Ordinary + Bias-Corrected)
      # --------------------------------------------------------
      CASE.CMA       <- BCA.cb(beta = beta, th = th, u = u, x1 = x1, x2 = x2, COV = COV)
      CASE.NIE[m]     <- CASE.CMA$NIE
      CASE.NIE_BC1[m] <- CASE.CMA$NIE_BC1
      CASE.NIE_BC2[m] <- CASE.CMA$NIE_BC2
      CASE.NDE[m]     <- CASE.CMA$NDE
      CASE.NDE_BC1[m] <- CASE.CMA$NDE_BC1
      CASE.NDE_BC2[m] <- CASE.CMA$NDE_BC2
      CASE.MP[m]      <- CASE.CMA$MP
      CASE.MP_BC11[m] <- CASE.CMA$MP_BC1_insert_bc1
      CASE.MP_BC12[m] <- CASE.CMA$MP_BC1_insert_bc2
      CASE.MP_BC21[m] <- CASE.CMA$MP_BC2_insert_bc1
      CASE.MP_BC22[m] <- CASE.CMA$MP_BC2_insert_bc2
      
      print(paste0("End ", i, " / ", N, " / ", m, " / ", M))
    }
    
    # --------------------------------------------------------
    # Save raw simulation outputs for this (i, N)
    # --------------------------------------------------------
    assign(
      paste0("LIST_CB_i", i, "_N", N),
      data.frame(
        NDE      = CASE.NDE,
        NDE_BC1  = CASE.NDE_BC1,
        NDE_BC2  = CASE.NDE_BC2,
        NIE      = CASE.NIE,
        NIE_BC1  = CASE.NIE_BC1,
        NIE_BC2  = CASE.NIE_BC2,
        MP       = CASE.MP,
        MP_BC11  = CASE.MP_BC11,
        MP_BC12  = CASE.MP_BC12,
        MP_BC21  = CASE.MP_BC21,
        MP_BC22  = CASE.MP_BC22
      )
    )
    write.csv(
      get(paste0("LIST_CB_i", i, "_N", N)),
      paste0(
        "./Sec4.Simulation/Simulation_data/02.CB_data/",
        paste0("LIST_CB_i", i, "_N", N), ".csv"
      ),
      row.names = FALSE
    )
  }
}


# ------------------------------------------------------------
# 3) Load Simulation Data & Prepare for Summary
# ------------------------------------------------------------
# Read back the saved CSV files for selected N and cases.
folder_path <- "./Sec4.Simulation/Simulation_data/02.CB_data/"

i_values     <- 1:3
N_values     <- c(200, 300, 600)
CB_data_list <- list()

# Load all (i, N) datasets
for (i in i_values) {
  for (N in N_values) {
    var_name  <- paste0("LIST_CB_i", i, "_N", N)
    file_path <- paste0(folder_path, var_name, ".csv")
    CB_data_list[[var_name]] <- read.csv(file_path)
  }
}


# ------------------------------------------------------------
# 4) Summary Calculations (Bias, Relative Bias, Variance)
# ------------------------------------------------------------
# Output tables for NDE, NIE, and MP (with BC variants).
result_NDE <- data.frame(
  Case = 0, N = 0, TrueNDE = 0,
  Est = 0, Var = 0, Bias = 0, RBias = 0,
  Est1 = 0, Var1 = 0, Bias1 = 0, RBias1 = 0,
  Est2 = 0, Var2 = 0, Bias2 = 0, RBias2 = 0
)
result_NIE <- data.frame(
  Case = 0, N = 0, TrueNIE = 0,
  Est = 0, Var = 0, Bias = 0, RBias = 0,
  Est1 = 0, Var1 = 0, Bias1 = 0, RBias1 = 0,
  Est2 = 0, Var2 = 0, Bias2 = 0, RBias2 = 0
)
result_MP <- data.frame(
  Case = 0, N = 0, TrueMP = 0,
  Est = 0, Var = 0, Bias = 0, RBias = 0,
  Est11 = 0, Var11 = 0, Bias11 = 0, RBias11 = 0,
  Est12 = 0, Var12 = 0, Bias12 = 0, RBias12 = 0, 
  Est21 = 0, Var21 = 0, Bias21 = 0, RBias21 = 0,
  Est22 = 0, Var22 = 0, Bias22 = 0, RBias22 = 0
)

# Wrapper to compute summary metrics for one dataset
# (i, N are parsed from the object/file names)
Simulation_summary_CB <- function(df, i, N) {
  
  #### NDE
  Est_NDE       <- mean(df$NDE, na.rm = TRUE)
  Var10_NDE     <- 10 * var(df$NDE, na.rm = TRUE)
  Bias_NDE      <- Est_NDE - true_value_CB[i, "TrueNDE"]
  Rbias100_NDE  <- 100 * (Est_NDE - true_value_CB[i, "TrueNDE"]) / true_value_CB[i, "TrueNDE"]
  
  Est_NDE_BC1   <- mean(df$NDE_BC1, na.rm = TRUE)
  Var10_NDE_BC1 <- 10 * var(df$NDE_BC1, na.rm = TRUE)
  Bias_NDE_BC1  <- Est_NDE_BC1 - true_value_CB[i, "TrueNDE"]
  Rbias100_NDE_BC1 <- 100 * (Est_NDE_BC1 - true_value_CB[i, "TrueNDE"]) / true_value_CB[i, "TrueNDE"]
  
  Est_NDE_BC2   <- mean(df$NDE_BC2, na.rm = TRUE)
  Var10_NDE_BC2 <- 10 * var(df$NDE_BC2, na.rm = TRUE)
  Bias_NDE_BC2  <- Est_NDE_BC2 - true_value_CB[i, "TrueNDE"]
  Rbias100_NDE_BC2 <- 100 * (Est_NDE_BC2 - true_value_CB[i, "TrueNDE"]) / true_value_CB[i, "TrueNDE"]
  
  #### NIE
  Est_NIE       <- mean(df$NIE, na.rm = TRUE)
  Var10_NIE     <- 10 * var(df$NIE, na.rm = TRUE)
  Bias_NIE      <- Est_NIE - true_value_CB[i, "TrueNIE"]
  Rbias100_NIE  <- 100 * (Est_NIE - true_value_CB[i, "TrueNIE"]) / true_value_CB[i, "TrueNIE"]
  
  Est_NIE_BC1   <- mean(df$NIE_BC1, na.rm = TRUE)
  Var10_NIE_BC1 <- 10 * var(df$NIE_BC1, na.rm = TRUE)
  Bias_NIE_BC1  <- Est_NIE_BC1 - true_value_CB[i, "TrueNIE"]
  Rbias100_NIE_BC1 <- 100 * (Est_NIE_BC1 - true_value_CB[i, "TrueNIE"]) / true_value_CB[i, "TrueNIE"]
  
  Est_NIE_BC2   <- mean(df$NIE_BC2, na.rm = TRUE)
  Var10_NIE_BC2 <- 10 * var(df$NIE_BC2, na.rm = TRUE)
  Bias_NIE_BC2  <- Est_NIE_BC2 - true_value_CB[i, "TrueNIE"]
  Rbias100_NIE_BC2 <- 100 * (Est_NIE_BC2 - true_value_CB[i, "TrueNIE"]) / true_value_CB[i, "TrueNIE"]
  
  #### MP
  Est_MP       <- mean(df$MP, na.rm = TRUE)
  Var10_MP     <- 10 * var(df$MP, na.rm = TRUE)
  Bias_MP      <- Est_MP - true_value_CB[i, "TrueMP"]
  Rbias100_MP  <- 100 * (Est_MP - true_value_CB[i, "TrueMP"]) / true_value_CB[i, "TrueMP"]
  
  Est_MP_BC11   <- mean(df$MP_BC11, na.rm = TRUE)
  Var10_MP_BC11 <- 10 * var(df$MP_BC11, na.rm = TRUE)
  Bias_MP_BC11  <- Est_MP_BC11 - true_value_CB[i, "TrueMP"]
  Rbias100_MP_BC11 <- 100 * (Est_MP_BC11 - true_value_CB[i, "TrueMP"]) / true_value_CB[i, "TrueMP"]
  
  Est_MP_BC12   <- mean(df$MP_BC12, na.rm = TRUE)
  Var10_MP_BC12 <- 10 * var(df$MP_BC12, na.rm = TRUE)
  Bias_MP_BC12  <- Est_MP_BC12 - true_value_CB[i, "TrueMP"]
  Rbias100_MP_BC12 <- 100 * (Est_MP_BC12 - true_value_CB[i, "TrueMP"]) / true_value_CB[i, "TrueMP"]
  
  Est_MP_BC21   <- mean(df$MP_BC21, na.rm = TRUE)
  Var10_MP_BC21 <- 10 * var(df$MP_BC21, na.rm = TRUE)
  Bias_MP_BC21  <- Est_MP_BC21 - true_value_CB[i, "TrueMP"]
  Rbias100_MP_BC21 <- 100 * (Est_MP_BC21 - true_value_CB[i, "TrueMP"]) / true_value_CB[i, "TrueMP"]
  
  Est_MP_BC22   <- mean(df$MP_BC22, na.rm = TRUE)
  Var10_MP_BC22 <- 10 * var(df$MP_BC22, na.rm = TRUE)
  Bias_MP_BC22  <- Est_MP_BC22 - true_value_CB[i, "TrueMP"]
  Rbias100_MP_BC22 <- 100 * (Est_MP_BC22 - true_value_CB[i, "TrueMP"]) / true_value_CB[i, "TrueMP"]
  
  # Assemble summary rows
  result_NDE <- data.frame(
    Case   = i, N = N,
    TrueNDE = t(true_value_CB[i, c(1)]),
    Est    = round(Est_NDE, 5),
    Var    = round(Var10_NDE, 5),
    Bias   = round(Bias_NDE, 5),
    RBias  = round(Rbias100_NDE, 5),
    Est1   = round(Est_NDE_BC1, 5),
    Var1   = round(Var10_NDE_BC1, 5),
    Bias1  = round(Bias_NDE_BC1, 5),
    RBias1 = round(Rbias100_NDE_BC1, 5),
    Est2   = round(Est_NDE_BC2, 5),
    Var2   = round(Var10_NDE_BC2, 5),
    Bias2  = round(Bias_NDE_BC2, 5),
    RBias2 = round(Rbias100_NDE_BC2, 5)
  )
  
  result_NIE <- data.frame(
    Case   = i, N = N,
    TrueNIE = t(true_value_CB[i, c(2)]),
    Est    = round(Est_NIE, 5),
    Var    = round(Var10_NIE, 5),
    Bias   = round(Bias_NIE, 5),
    RBias  = round(Rbias100_NIE, 5),
    Est1   = round(Est_NIE_BC1, 5),
    Var1   = round(Var10_NIE_BC1, 5),
    Bias1  = round(Bias_NIE_BC1, 5),
    RBias1 = round(Rbias100_NIE_BC1, 5),
    Est2   = round(Est_NIE_BC2, 5),
    Var2   = round(Var10_NIE_BC2, 5),
    Bias2  = round(Bias_NIE_BC2, 5),
    RBias2 = round(Rbias100_NIE_BC2, 5)
  )
  
  result_MP <- data.frame(
    Case   = i, N = N,
    TrueMP = t(true_value_CB[i, c(3)]),
    Est    = round(Est_MP, 5),
    Var    = round(Var10_MP, 5),
    Bias   = round(Bias_MP, 5),
    RBias  = round(Rbias100_MP, 5),
    Est11  = round(Est_MP_BC11, 5),
    Var11  = round(Var10_MP_BC11, 5),
    Bias11 = round(Bias_MP_BC11, 5),
    RBias11 = round(Rbias100_MP_BC11, 5),
    Est12  = round(Est_MP_BC12, 5),
    Var12  = round(Var10_MP_BC12, 5),
    Bias12 = round(Bias_MP_BC12, 5),
    RBias12 = round(Rbias100_MP_BC12, 5),
    Est21  = round(Est_MP_BC21, 5),
    Var21  = round(Var10_MP_BC21, 5),
    Bias21 = round(Bias_MP_BC21, 5),
    RBias21 = round(Rbias100_MP_BC21, 5),
    Est22  = round(Est_MP_BC22, 5),
    Var22  = round(Var10_MP_BC22, 5),
    Bias22 = round(Bias_MP_BC22, 5),
    RBias22 = round(Rbias100_MP_BC22, 5)
  )
  
  return(list(result_NDE, result_NIE, result_MP))
}

# Run summaries across all loaded datasets
names(CB_data_list)
final_list <- Simulation_summary_CB(
  CB_data_list[[1]],
  i = as.numeric(str_extract(names(CB_data_list)[1], "(?<=_i)\\d+")),
  N = as.numeric(str_extract(names(CB_data_list)[1], "(?<=_N)\\d+"))
)
for (i in 1:length(CB_data_list)) {
  temp_list <- Simulation_summary_CB(
    CB_data_list[[i]],
    i = as.numeric(str_extract(names(CB_data_list)[i], "(?<=_i)\\d+")),
    N = as.numeric(str_extract(names(CB_data_list)[i], "(?<=_N)\\d+"))
  )
  final_list[[1]] <- rbind(final_list[[1]], temp_list[[1]])
  final_list[[2]] <- rbind(final_list[[2]], temp_list[[2]])
  final_list[[3]] <- rbind(final_list[[3]], temp_list[[3]])
}

# Reformat to publication-ready CSVs
result_NDE_CB <- make_csvfile(final_list[[1]])
result_NIE_CB <- make_csvfile(final_list[[2]])
result_MP_CB  <- make_csvfile_MP(final_list[[3]]) %>% dplyr::filter(Method %!in% c("Est11", "Est21"))

write.csv(result_NDE_CB, paste0("./Sec4.Simulation/Simulation_result/02.CB_result/Simul_result_CB_NDE.csv"))
write.csv(result_NIE_CB, paste0("./Sec4.Simulation/Simulation_result/02.CB_result/Simul_result_CB_NIE.csv"))
write.csv(result_MP_CB,  paste0("./Sec4.Simulation/Simulation_result/02.CB_result/Simul_result_CB_MP.csv"))


# ------------------------------------------------------------
# 5) Outlier Inspection and Filtering
# ------------------------------------------------------------
# Identify extreme simulation values and remove them for stability. 
# Thresholds are manually set based on exploratory analysis.
# ------------------------------------------------------------
# Not applied here. Add your own thresholds/filters if needed.
# ------------------------------------------------------------
