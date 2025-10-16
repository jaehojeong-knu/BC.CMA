# ============================================================
# 01.BC_Data_application.R
# ------------------------------------------------------------
# Data application for the Binary Outcome + Continuous Mediator (BC) case
# corresponding to Section 5.1 of the manuscript.
#
# Dataset:
#   - We use the UPBdata dataset bundled with the {medflex} package.
#     No external data files are required.
#
# Reference:
#   Steen, J., Loeys, T., Moerkerke, B., & Vansteelandt, S. (2017).
#   "Medflex: an R package for flexible mediation analysis using natural
#   effect models." Journal of Statistical Software, 76, 1–46.
#
# This script performs:
#   1) Data preparation and modeling formulas
#   2) Point estimation of NDE, NIE, MP and bias-corrected variants
#   3) Bootstrap inference (percentile 95% CI) for all quantities
#
# Dependencies:
#   - Requires 00.basic_functions.R and 03.BC_functions.R
# ============================================================

source('./functions/00.basic_functions.R')
source('./functions/03.BC_functions.R')

# ------------------------------------------------------------
# 1) Load & Prepare Data
# ------------------------------------------------------------
# Variables kept:
#   Y (binary):  UPB      - unwanted pursuit behavior (0/1)
#   X (binary):  attbin   - attitude (1 = higher-than-mean, 0 = lower-than-mean)
#   M (cont.) :  negaff   - negative affectivity (standardized)
#   C (conf.) :  age (numeric), educ (1–3), gender (factor: 1/2)
# ------------------------------------------------------------

data("UPBdata", package = "medflex")

modeling_data <- UPBdata %>%
  dplyr::select(UPB, attbin, negaff, age, educ, gender)

# Outcome & exposure to factors (explicit binary)
modeling_data$UPB    <- as.factor(modeling_data$UPB)     # 0/1
modeling_data$attbin <- as.factor(modeling_data$attbin)  # 0/1

# Confounders
# age: numeric (keep as is; visualize / summarize if needed)
# educ: convert to numeric scale 1(L)–3(H)
modeling_data$educ <- as.numeric(modeling_data$educ)

# gender: ensure factor coding (original levels -> numeric -> factor)
# NOTE: level meanings follow the package coding; 1/2 labels are retained.
modeling_data$gender <- as.numeric(modeling_data$gender)
modeling_data$gender <- as.factor(modeling_data$gender)

# ------------------------------------------------------------
# 2) Modeling Setup (BC case: Binary Y, Continuous M)
# ------------------------------------------------------------
# Conditional effects at selected covariate combinations:
#   - gender: {1, 2}
#   - age   : sample quantiles (Q1, median, Q3)
#   - educ  : {1, 2, 3}
# This yields a 2 (gender) × 3 (age) × 3 (educ) grid per "3×3 table by gender".
# ------------------------------------------------------------

med_formula <- as.formula("negaff ~ attbin + gender + age + educ")
out_formula <- as.formula("UPB ~ attbin + negaff + attbin:negaff + gender + age + educ")

age_grid   <- c(summary(modeling_data$age)[2],  # Q1
                summary(modeling_data$age)[3],  # median
                summary(modeling_data$age)[5])  # Q3
educ_grid  <- c(1, 2, 3)
gender_grid <- c(1, 2)

u_combi <- expand.grid(gender = gender_grid, age = age_grid, educ = educ_grid)

# ------------------------------------------------------------
# 3) Point Estimation
# ------------------------------------------------------------
# For each covariate combination u, we compute:
#   - Ordinary: NDE, NIE, MP
#   - Bias-corrected: NDE_BC1/BC2, NIE_BC1/BC2, MP_BC1/BC2
# ------------------------------------------------------------

# Fit mediator model (linear) and outcome model (logistic)
Med <- lm(med_formula, data = modeling_data)
th  <- summary(Med)$coefficients[, 1]      # theta estimates
sig_theta <- vcov(Med)                      # Var-Cov(theta)
sig2_med  <- summary(Med)$sigma^2           # residual variance of mediator model

Outcome <- glm(out_formula, family = binomial, data = modeling_data)
beta <- summary(Outcome)$coefficients[, 1]  # beta estimates
# Reorder to (Intercept, X, M, X:M, U's) as used by BC functions
beta <- c(beta[1], beta[2], beta[3], beta[length(beta)], beta[4:(length(beta) - 1)])
sig_beta <- vcov(Outcome)                    # Var-Cov(beta)

# Joint covariance for (beta, theta); 'u' here provides the confounder design
conf <- modeling_data[, 4:(length(beta) - 1)]
COV  <- VAR.C(beta = beta, th = th, u = conf, sig.beta = sig_beta, sig.th = sig_theta)

# Container for point estimates across all (gender, age, educ)
Result_Ord <- data.frame(
  gender = 0, age = 0, educ = 0,
  NDE = 0, NIE = 0, MP = 0,
  NDE_BC1 = 0, NDE_BC2 = 0,
  NIE_BC1 = 0, NIE_BC2 = 0,
  MP_BC1 = 0, MP_BC2 = 0
)

for (i in 1:nrow(u_combi)) {
  u_unit <- u_combi[i, ]
  
  est <- CMA.bc(beta = beta, th = th, sig2 = sig2_med,
                u = u_unit, x1 = 1, x2 = 0)
  bc  <- BCA.bc(beta = beta, th = th,
                u = u_unit, x1 = 1, x2 = 0,
                sig2 = sig2_med, COV = COV)
  
  Result_Ord[i, ] <- unlist(c(
    u_unit$gender, u_unit$age, u_unit$educ,
    bc["NDE"], bc["NIE"], bc["MP"],
    bc["NDE_BC1"], bc["NDE_BC2"],
    bc["NIE_BC1"], bc["NIE_BC2"],
    bc["MP_BC1"], bc["MP_BC2"]
  ))
}

Result_Ord <- Result_Ord %>% arrange(gender, age, educ)

# ------------------------------------------------------------
# 4) Bootstrap Inference (percentile 95% CI)
# ------------------------------------------------------------
# We bootstrap the full procedure (refit Med/Outcome; recompute COV, CMA, BCA)
# for each covariate combination u. Returns percentile 95% CIs.
# ------------------------------------------------------------

n <- nrow(modeling_data)  # 385
Result_boot <- data.frame()

for (i in 1:nrow(u_combi)) {
  u_unit <- u_combi[i, ]
  
  mediation_stat <- function(data, indices) {
    data_b <- data[indices, ]
    
    # Refit mediator model
    Med_b <- lm(med_formula, data = data_b)
    th_b  <- summary(Med_b)$coefficients[, 1]
    sig_theta_b <- vcov(Med_b)
    sig2_med_b  <- summary(Med_b)$sigma^2
    
    # Refit outcome model
    Outcome_b <- glm(out_formula, family = binomial, data = data_b)
    beta_b <- summary(Outcome_b)$coefficients[, 1]
    beta_b <- c(beta_b[1], beta_b[2], beta_b[3], beta_b[length(beta_b)], beta_b[4:(length(beta_b) - 1)])
    sig_beta_b <- vcov(Outcome_b)
    
    # Joint covariance
    conf_b <- data_b[, 4:(length(beta_b) - 1)]
    COV_b  <- VAR.C(beta = beta_b, th = th_b, u = conf_b,
                    sig.beta = sig_beta_b, sig.th = sig_theta_b)
    
    # Recompute ordinary & bias-corrected effects at u_unit
    est_b <- CMA.bc(beta = beta_b, th = th_b, sig2 = sig2_med_b,
                    u = u_unit, x1 = 1, x2 = 0)
    bc_b  <- BCA.bc(beta = beta_b, th = th_b,
                    u = u_unit, x1 = 1, x2 = 0,
                    sig2 = sig2_med_b, COV = COV_b)
    
    # Return vector: (NDE, NIE, MP, NDE_BC1, NDE_BC2, NIE_BC1, NIE_BC2, MP_BC1, MP_BC2)
    return(unlist(c(est_b["NDE"], est_b["NIE"], est_b["MP"],
                    bc_b["NDE_BC1"], bc_b["NDE_BC2"],
                    bc_b["NIE_BC1"], bc_b["NIE_BC2"],
                    bc_b["MP_BC1"], bc_b["MP_BC2"])))
  }
  
  set.seed(i)
  boot_res <- boot(data = modeling_data, statistic = mediation_stat, R = 1000)
  
  # Percentile bootstrap 95% CIs
  ci_NDE     <- boot.ci(boot_res, index = 1, type = "perc")$perc[4:5]
  ci_NIE     <- boot.ci(boot_res, index = 2, type = "perc")$perc[4:5]
  ci_MP      <- boot.ci(boot_res, index = 3, type = "perc")$perc[4:5]
  ci_NDE_BC1 <- boot.ci(boot_res, index = 4, type = "perc")$perc[4:5]
  ci_NDE_BC2 <- boot.ci(boot_res, index = 5, type = "perc")$perc[4:5]
  ci_NIE_BC1 <- boot.ci(boot_res, index = 6, type = "perc")$perc[4:5]
  ci_NIE_BC2 <- boot.ci(boot_res, index = 7, type = "perc")$perc[4:5]
  ci_MP_BC1  <- boot.ci(boot_res, index = 8, type = "perc")$perc[4:5]
  ci_MP_BC2  <- boot.ci(boot_res, index = 9, type = "perc")$perc[4:5]
  
  # Collect point estimates (t0) and CIs
  temp <- data.frame(
    gender = u_unit$gender,
    age    = u_unit$age,
    educ   = u_unit$educ,
    
    NDE       = boot_res$t0[1],
    NDE_CI_L  = ci_NDE[1],     NDE_CI_U  = ci_NDE[2],
    NIE       = boot_res$t0[2],
    NIE_CI_L  = ci_NIE[1],     NIE_CI_U  = ci_NIE[2],
    MP        = boot_res$t0[3],
    MP_CI_L   = ci_MP[1],      MP_CI_U   = ci_MP[2],
    
    NDE_BC1      = boot_res$t0[4],
    NDE_BC1_CI_L = ci_NDE_BC1[1], NDE_BC1_CI_U = ci_NDE_BC1[2],
    NDE_BC2      = boot_res$t0[5],
    NDE_BC2_CI_L = ci_NDE_BC2[1], NDE_BC2_CI_U = ci_NDE_BC2[2],
    
    NIE_BC1      = boot_res$t0[6],
    NIE_BC1_CI_L = ci_NIE_BC1[1], NIE_BC1_CI_U = ci_NIE_BC1[2],
    NIE_BC2      = boot_res$t0[7],
    NIE_BC2_CI_L = ci_NIE_BC2[1], NIE_BC2_CI_U = ci_NIE_BC2[2],
    
    MP_BC1      = boot_res$t0[8],
    MP_BC1_CI_L = ci_MP_BC1[1], MP_BC1_CI_U = ci_MP_BC1[2],
    MP_BC2      = boot_res$t0[9],
    MP_BC2_CI_L = ci_MP_BC2[1], MP_BC2_CI_U = ci_MP_BC2[2]
  )
  
  Result_boot <- rbind(Result_boot, temp)
  message(sprintf("%d / %d completed", i, nrow(u_combi)))
}

rownames(Result_boot) <- NULL
Result_boot <- Result_boot %>% arrange(gender, age, educ)

# ------------------------------------------------------------
# (Optional) Save outputs
# ------------------------------------------------------------
save(Result_boot,file="./Sec5.Data_application/Data_application_result/01.BC_application.RData")
