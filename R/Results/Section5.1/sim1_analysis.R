###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####            R-code for reproducing Tables 2-4 in Section 5.1.              ####
###################################################################################

# The following R-code contains functions for calculating TPR, FPR, FDR, and MSE.

# Each ".rda" file contains the modified BIC values, optimal estimators('pi', 'Bk', 'sigma'), 
# and the estimators of mixing proportion 'w'.

### Data ###
source("./Data/simulation_data.R")

### Functions ###
source("./Scripts/Functions/simulation_summary_functions.R")

### Model1 ###
## results ##
load("./Results/Section5.1/Model1/sim1_mod1_flexmix.rda")
load("./Results/Section5.1/Model1/sim1_mod1_oracle.rda")
load("./Results/Section5.1/Model1/sim1_mod1_mvFMR.rda")
load("./Results/Section5.1/Model1/sim1_mod1_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model1/sim1_mod1_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model1/sim1_mod1_PGM_mvFMR_MCP.rda")
load("./Results/Section5.1/Model1/sim1_mod1_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model1/sim1_mod1_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model1/sim1_mod1_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.1.1(1, 500)$true
true_1000 <- data_generate_5.1.1(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_indepX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_flexmix_corX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_corX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_corX_500 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_flexmix_indepX_500 <- do.call(rbind, lapply(sim1_mod1_flexmix$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_flexmix_corX_500 <- do.call(rbind, lapply(sim1_mod1_flexmix$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim1_mod1_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim1_mod1_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim1_mod1_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim1_mod1_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim1_mod1_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim1_mod1_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim1_mod1_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim1_mod1_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim1_mod1_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim1_mod1_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod1_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod1_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model1_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model1_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model1_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model1_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model1_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model1_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model1_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model1_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model1_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model1_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model1_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model1_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model1_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model1_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model1_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model1_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model1_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model1_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model1_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model1_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model1_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model1_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model1_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model1_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model1_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model1_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model1_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model1_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model1_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model1_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model1_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model1_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model1_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model1_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod1_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod1_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$TPR, model1_table$indepX$n500$sd$TPR)
mod1_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$FPR, model1_table$indepX$n500$sd$FPR)
mod1_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$FDR, model1_table$indepX$n500$sd$FDR)
mod1_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$MSE_pi, model1_table$indepX$n500$sd$MSE_pi)
mod1_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$MSE_Bk, model1_table$indepX$n500$sd$MSE_Bk)
mod1_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model1_table$indepX$n500$mean$MSE_sigma, model1_table$indepX$n500$sd$MSE_sigma)

mod1_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$TPR, model1_table$indepX$n1000$sd$TPR)
mod1_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$FPR, model1_table$indepX$n1000$sd$FPR)
mod1_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$FDR, model1_table$indepX$n1000$sd$FDR)
mod1_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$MSE_pi, model1_table$indepX$n1000$sd$MSE_pi)
mod1_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$MSE_Bk, model1_table$indepX$n1000$sd$MSE_Bk)
mod1_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model1_table$indepX$n1000$mean$MSE_sigma, model1_table$indepX$n1000$sd$MSE_sigma)

mod1_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$TPR, model1_table$corX$n500$sd$TPR)
mod1_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$FPR, model1_table$corX$n500$sd$FPR)
mod1_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$FDR, model1_table$corX$n500$sd$FDR)
mod1_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$MSE_pi, model1_table$corX$n500$sd$MSE_pi)
mod1_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$MSE_Bk, model1_table$corX$n500$sd$MSE_Bk)
mod1_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model1_table$corX$n500$mean$MSE_sigma, model1_table$corX$n500$sd$MSE_sigma)

mod1_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$TPR, model1_table$corX$n1000$sd$TPR)
mod1_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$FPR, model1_table$corX$n1000$sd$FPR)
mod1_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$FDR, model1_table$corX$n1000$sd$FDR)
mod1_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$MSE_pi, model1_table$corX$n1000$sd$MSE_pi)
mod1_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$MSE_Bk, model1_table$corX$n1000$sd$MSE_Bk)
mod1_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model1_table$corX$n1000$mean$MSE_sigma, model1_table$corX$n1000$sd$MSE_sigma)
# save(mod1_table, file="./Results/Section5.1/Model1/mod1_table.rda")


### Model2 ###
## results ##
load("./Results/Section5.1/Model2/sim1_mod2_flexmix.rda")
load("./Results/Section5.1/Model2/sim1_mod2_oracle.rda")
load("./Results/Section5.1/Model2/sim1_mod2_mvFMR.rda")
load("./Results/Section5.1/Model2/sim1_mod2_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model2/sim1_mod2_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model2/sim1_mod2_PGM_mvFMR_MCP.rda")
load("./Results/Section5.1/Model2/sim1_mod2_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model2/sim1_mod2_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model2/sim1_mod2_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.1.2(1, 500)$true
true_1000 <- data_generate_5.1.2(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_indepX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_flexmix_corX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_corX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_corX_500 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_flexmix_indepX_500 <- do.call(rbind, lapply(sim1_mod2_flexmix$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_flexmix_corX_500 <- do.call(rbind, lapply(sim1_mod2_flexmix$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim1_mod2_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim1_mod2_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim1_mod2_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim1_mod2_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim1_mod2_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim1_mod2_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim1_mod2_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim1_mod2_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim1_mod2_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim1_mod2_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod2_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod2_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model2_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model2_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model2_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model2_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model2_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model2_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model2_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model2_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model2_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model2_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model2_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model2_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model2_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model2_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model2_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model2_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model2_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model2_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model2_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model2_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model2_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model2_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model2_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model2_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model2_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model2_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model2_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model2_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model2_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model2_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model2_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model2_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model2_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model2_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod2_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod2_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$TPR, model2_table$indepX$n500$sd$TPR)
mod2_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$FPR, model2_table$indepX$n500$sd$FPR)
mod2_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$FDR, model2_table$indepX$n500$sd$FDR)
mod2_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$MSE_pi, model2_table$indepX$n500$sd$MSE_pi)
mod2_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$MSE_Bk, model2_table$indepX$n500$sd$MSE_Bk)
mod2_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model2_table$indepX$n500$mean$MSE_sigma, model2_table$indepX$n500$sd$MSE_sigma)

mod2_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$TPR, model2_table$indepX$n1000$sd$TPR)
mod2_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$FPR, model2_table$indepX$n1000$sd$FPR)
mod2_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$FDR, model2_table$indepX$n1000$sd$FDR)
mod2_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$MSE_pi, model2_table$indepX$n1000$sd$MSE_pi)
mod2_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$MSE_Bk, model2_table$indepX$n1000$sd$MSE_Bk)
mod2_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model2_table$indepX$n1000$mean$MSE_sigma, model2_table$indepX$n1000$sd$MSE_sigma)

mod2_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$TPR, model2_table$corX$n500$sd$TPR)
mod2_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$FPR, model2_table$corX$n500$sd$FPR)
mod2_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$FDR, model2_table$corX$n500$sd$FDR)
mod2_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$MSE_pi, model2_table$corX$n500$sd$MSE_pi)
mod2_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$MSE_Bk, model2_table$corX$n500$sd$MSE_Bk)
mod2_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model2_table$corX$n500$mean$MSE_sigma, model2_table$corX$n500$sd$MSE_sigma)

mod2_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$TPR, model2_table$corX$n1000$sd$TPR)
mod2_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$FPR, model2_table$corX$n1000$sd$FPR)
mod2_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$FDR, model2_table$corX$n1000$sd$FDR)
mod2_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$MSE_pi, model2_table$corX$n1000$sd$MSE_pi)
mod2_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$MSE_Bk, model2_table$corX$n1000$sd$MSE_Bk)
mod2_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model2_table$corX$n1000$mean$MSE_sigma, model2_table$corX$n1000$sd$MSE_sigma)
# save(mod2_table, file="./Results/Section5.1/Model2/mod2_table.rda")


### Model3 ###
## results ##
load("./Results/Section5.1/Model3/sim1_mod3_flexmix.rda")
load("./Results/Section5.1/Model3/sim1_mod3_oracle.rda")
load("./Results/Section5.1/Model3/sim1_mod3_mvFMR.rda")
load("./Results/Section5.1/Model3/sim1_mod3_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model3/sim1_mod3_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model3/sim1_mod3_PGM_mvFMR_MCP.rda")
load("./Results/Section5.1/Model3/sim1_mod3_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.1/Model3/sim1_mod3_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.1/Model3/sim1_mod3_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.1.3(1, 500)$true
true_1000 <- data_generate_5.1.3(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_indepX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_indepX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_flexmix_corX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_corX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_corX_500 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_flexmix_indepX_500 <- do.call(rbind, lapply(sim1_mod3_flexmix$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_flexmix_corX_500 <- do.call(rbind, lapply(sim1_mod3_flexmix$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim1_mod3_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim1_mod3_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim1_mod3_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim1_mod3_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim1_mod3_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim1_mod3_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim1_mod3_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim1_mod3_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim1_mod3_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim1_mod3_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod3_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim1_mod3_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model3_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model3_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model3_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model3_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model3_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model3_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model3_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model3_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model3_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model3_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model3_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model3_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model3_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model3_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model3_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model3_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model3_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model3_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model3_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model3_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model3_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model3_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model3_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model3_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model3_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model3_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model3_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model3_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model3_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model3_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model3_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model3_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model3_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model3_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod3_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod3_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$TPR, model3_table$indepX$n500$sd$TPR)
mod3_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$FPR, model3_table$indepX$n500$sd$FPR)
mod3_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$FDR, model3_table$indepX$n500$sd$FDR)
mod3_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$MSE_pi, model3_table$indepX$n500$sd$MSE_pi)
mod3_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$MSE_Bk, model3_table$indepX$n500$sd$MSE_Bk)
mod3_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model3_table$indepX$n500$mean$MSE_sigma, model3_table$indepX$n500$sd$MSE_sigma)

mod3_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$TPR, model3_table$indepX$n1000$sd$TPR)
mod3_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$FPR, model3_table$indepX$n1000$sd$FPR)
mod3_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$FDR, model3_table$indepX$n1000$sd$FDR)
mod3_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$MSE_pi, model3_table$indepX$n1000$sd$MSE_pi)
mod3_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$MSE_Bk, model3_table$indepX$n1000$sd$MSE_Bk)
mod3_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model3_table$indepX$n1000$mean$MSE_sigma, model3_table$indepX$n1000$sd$MSE_sigma)

mod3_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$TPR, model3_table$corX$n500$sd$TPR)
mod3_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$FPR, model3_table$corX$n500$sd$FPR)
mod3_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$FDR, model3_table$corX$n500$sd$FDR)
mod3_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$MSE_pi, model3_table$corX$n500$sd$MSE_pi)
mod3_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$MSE_Bk, model3_table$corX$n500$sd$MSE_Bk)
mod3_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model3_table$corX$n500$mean$MSE_sigma, model3_table$corX$n500$sd$MSE_sigma)

mod3_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$TPR, model3_table$corX$n1000$sd$TPR)
mod3_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$FPR, model3_table$corX$n1000$sd$FPR)
mod3_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$FDR, model3_table$corX$n1000$sd$FDR)
mod3_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$MSE_pi, model3_table$corX$n1000$sd$MSE_pi)
mod3_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$MSE_Bk, model3_table$corX$n1000$sd$MSE_Bk)
mod3_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model3_table$corX$n1000$mean$MSE_sigma, model3_table$corX$n1000$sd$MSE_sigma)
# save(mod3_table, file="./Results/Section5.1/Model3/mod3_table.rda")