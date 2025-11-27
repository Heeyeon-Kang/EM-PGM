###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####            R-code for reproducing Tables 6-8 in Section 5.2.              ####
###################################################################################

# The following R-code contains functions for calculating TPR, FPR, FDR, and MSE.

# Each ".rda" file contains the modified BIC values, optimal estimators('pi', 'Bk', 'sigma'), 
# and the estimators of mixing proportion 'w'.


### Data ###
source("./Data/simulation_data.R")

### Functions ###
source("./Scripts/Functions/simulation_summary_functions.R")

### Model 4 ###
## results ##
load("./Results/Section5.2/Model4/sim2_mod4_flexmix.rda")
load("./Results/Section5.2/Model4/sim2_mod4_oracle.rda")
load("./Results/Section5.2/Model4/sim2_mod4_mvFMR.rda")
load("./Results/Section5.2/Model4/sim2_mod4_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model4/sim2_mod4_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model4/sim2_mod4_PGM_mvFMR_MCP.rda")
load("./Results/Section5.2/Model4/sim2_mod4_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model4/sim2_mod4_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model4/sim2_mod4_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.2.1(1, 500)$true
true_1000 <- data_generate_5.2.1(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_indepX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_indepX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_flexmix_corX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_corX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_corX_500 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_flexmix_indepX_500 <- do.call(rbind, lapply(sim2_mod4_flexmix$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_flexmix_corX_500 <- do.call(rbind, lapply(sim2_mod4_flexmix$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim2_mod4_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim2_mod4_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim2_mod4_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim2_mod4_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim2_mod4_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim2_mod4_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim2_mod4_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim2_mod4_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim2_mod4_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim2_mod4_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod4_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod4_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model4_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model4_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model4_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model4_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model4_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model4_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model4_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model4_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model4_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model4_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model4_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model4_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model4_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model4_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model4_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model4_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model4_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model4_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model4_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model4_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model4_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model4_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model4_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model4_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model4_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model4_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model4_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model4_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model4_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model4_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model4_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model4_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model4_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model4_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod4_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod4_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$TPR, model4_table$indepX$n500$sd$TPR)
mod4_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$FPR, model4_table$indepX$n500$sd$FPR)
mod4_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$FDR, model4_table$indepX$n500$sd$FDR)
mod4_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$MSE_pi, model4_table$indepX$n500$sd$MSE_pi)
mod4_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$MSE_Bk, model4_table$indepX$n500$sd$MSE_Bk)
mod4_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model4_table$indepX$n500$mean$MSE_sigma, model4_table$indepX$n500$sd$MSE_sigma)

mod4_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$TPR, model4_table$indepX$n1000$sd$TPR)
mod4_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$FPR, model4_table$indepX$n1000$sd$FPR)
mod4_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$FDR, model4_table$indepX$n1000$sd$FDR)
mod4_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$MSE_pi, model4_table$indepX$n1000$sd$MSE_pi)
mod4_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$MSE_Bk, model4_table$indepX$n1000$sd$MSE_Bk)
mod4_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model4_table$indepX$n1000$mean$MSE_sigma, model4_table$indepX$n1000$sd$MSE_sigma)

mod4_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$TPR, model4_table$corX$n500$sd$TPR)
mod4_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$FPR, model4_table$corX$n500$sd$FPR)
mod4_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$FDR, model4_table$corX$n500$sd$FDR)
mod4_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$MSE_pi, model4_table$corX$n500$sd$MSE_pi)
mod4_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$MSE_Bk, model4_table$corX$n500$sd$MSE_Bk)
mod4_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model4_table$corX$n500$mean$MSE_sigma, model4_table$corX$n500$sd$MSE_sigma)

mod4_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$TPR, model4_table$corX$n1000$sd$TPR)
mod4_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$FPR, model4_table$corX$n1000$sd$FPR)
mod4_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$FDR, model4_table$corX$n1000$sd$FDR)
mod4_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$MSE_pi, model4_table$corX$n1000$sd$MSE_pi)
mod4_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$MSE_Bk, model4_table$corX$n1000$sd$MSE_Bk)
mod4_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model4_table$corX$n1000$mean$MSE_sigma, model4_table$corX$n1000$sd$MSE_sigma)
# save(mod4_table, file="./Results/Section5.2/Model4/mod4_table.rda")


### Model5 ###
## results ##
load("./Results/Section5.2/Model5/sim2_mod5_flexmix.rda")
load("./Results/Section5.2/Model5/sim2_mod5_oracle.rda")
load("./Results/Section5.2/Model5/sim2_mod5_mvFMR.rda")
load("./Results/Section5.2/Model5/sim2_mod5_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model5/sim2_mod5_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model5/sim2_mod5_PGM_mvFMR_MCP.rda")
load("./Results/Section5.2/Model5/sim2_mod5_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model5/sim2_mod5_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model5/sim2_mod5_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.2.2(1, 500)$true
true_1000 <- data_generate_5.2.2(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_indepX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_indepX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_flexmix_corX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_flexmix_corX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_flexmix_corX_500 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_flexmix_indepX_500 <- do.call(rbind, lapply(sim2_mod5_flexmix$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_flexmix_corX_500 <- do.call(rbind, lapply(sim2_mod5_flexmix$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim2_mod5_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim2_mod5_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim2_mod5_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim2_mod5_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim2_mod5_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim2_mod5_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim2_mod5_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim2_mod5_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim2_mod5_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim2_mod5_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod5_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod5_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model5_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model5_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model5_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model5_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model5_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model5_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model5_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model5_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model5_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model5_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model5_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model5_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model5_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model5_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model5_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model5_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model5_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model5_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model5_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model5_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model5_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model5_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model5_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model5_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model5_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model5_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model5_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model5_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model5_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model5_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model5_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model5_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model5_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model5_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod5_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod5_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$TPR, model5_table$indepX$n500$sd$TPR)
mod5_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$FPR, model5_table$indepX$n500$sd$FPR)
mod5_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$FDR, model5_table$indepX$n500$sd$FDR)
mod5_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$MSE_pi, model5_table$indepX$n500$sd$MSE_pi)
mod5_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$MSE_Bk, model5_table$indepX$n500$sd$MSE_Bk)
mod5_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model5_table$indepX$n500$mean$MSE_sigma, model5_table$indepX$n500$sd$MSE_sigma)

mod5_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$TPR, model5_table$indepX$n1000$sd$TPR)
mod5_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$FPR, model5_table$indepX$n1000$sd$FPR)
mod5_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$FDR, model5_table$indepX$n1000$sd$FDR)
mod5_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$MSE_pi, model5_table$indepX$n1000$sd$MSE_pi)
mod5_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$MSE_Bk, model5_table$indepX$n1000$sd$MSE_Bk)
mod5_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model5_table$indepX$n1000$mean$MSE_sigma, model5_table$indepX$n1000$sd$MSE_sigma)

mod5_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$TPR, model5_table$corX$n500$sd$TPR)
mod5_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$FPR, model5_table$corX$n500$sd$FPR)
mod5_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$FDR, model5_table$corX$n500$sd$FDR)
mod5_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$MSE_pi, model5_table$corX$n500$sd$MSE_pi)
mod5_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$MSE_Bk, model5_table$corX$n500$sd$MSE_Bk)
mod5_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model5_table$corX$n500$mean$MSE_sigma, model5_table$corX$n500$sd$MSE_sigma)

mod5_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$TPR, model5_table$corX$n1000$sd$TPR)
mod5_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$FPR, model5_table$corX$n1000$sd$FPR)
mod5_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$FDR, model5_table$corX$n1000$sd$FDR)
mod5_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$MSE_pi, model5_table$corX$n1000$sd$MSE_pi)
mod5_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$MSE_Bk, model5_table$corX$n1000$sd$MSE_Bk)
mod5_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model5_table$corX$n1000$mean$MSE_sigma, model5_table$corX$n1000$sd$MSE_sigma)
# save(mod5_table, file="./Results/Section5.2/Model5/mod5_table.rda")


### Model6 ###
## results ##
load("./Results/Section5.2/Model6/sim2_mod6_flexmix.rda")
load("./Results/Section5.2/Model6/sim2_mod6_oracle.rda")
load("./Results/Section5.2/Model6/sim2_mod6_mvFMR.rda")
load("./Results/Section5.2/Model6/sim2_mod6_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model6/sim2_mod6_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model6/sim2_mod6_PGM_mvFMR_MCP.rda")
load("./Results/Section5.2/Model6/sim2_mod6_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.2/Model6/sim2_mod6_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.2/Model6/sim2_mod6_ADMM_mvFMR_MCP.rda")
true_500 <- data_generate_5.2.3(1, 500)$true
true_1000 <- data_generate_5.2.3(1, 1000)$true
K <- 2

## flexmix ##
TPR_flexmix_indepX_500 <- rep(0, 500)
FPR_flexmix_indepX_500 <- rep(0, 500)
FDR_flexmix_indepX_500 <- rep(0, 500)
TPR_flexmix_corX_500 <- rep(0, 500)
FPR_flexmix_corX_500 <- rep(0, 500)
FDR_flexmix_corX_500 <- rep(0, 500)
MSE_flexmix_indepX_500 <- data.frame(pi = rep(0, 500), Bk = rep(0, 500), sigma = rep(0, 500))
MSE_flexmix_corX_500 <- data.frame(pi = rep(0, 500), Bk = rep(0, 500), sigma = rep(0, 500))

TPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_indepX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_flexmix_corX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_flexmix_corX_1000 <- unlist(lapply(sim2_mod6_flexmix$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_flexmix_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_flexmix$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_flexmix_corX_1000 <- do.call(rbind, lapply(sim2_mod6_flexmix$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))


## oracle ##
TPR_oracle_indepX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_indepX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_indepX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_oracle_corX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_oracle_corX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_oracle_corX_500 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_oracle_indepX_500 <- do.call(rbind, lapply(sim2_mod6_oracle$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_oracle_corX_500 <- do.call(rbind, lapply(sim2_mod6_oracle$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_indepX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_indepX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_oracle_corX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_oracle_corX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_oracle_corX_1000 <- unlist(lapply(sim2_mod6_oracle$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_oracle_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_oracle$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_oracle_corX_1000 <- do.call(rbind, lapply(sim2_mod6_oracle$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## mvFMR ##
TPR_unpen_indepX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_indepX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_indepX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n500, function(x) FDR_no_penalty(true_500, x)))
TPR_unpen_corX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n500, function(x) TPR_no_penalty(true_500, x)))
FPR_unpen_corX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n500, function(x) FPR_no_penalty(true_500, x)))
FDR_unpen_corX_500 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n500, function(x) FDR_no_penalty(true_500, x)))
MSE_unpen_indepX_500 <- do.call(rbind, lapply(sim2_mod6_mvFMR$optimal$indepX$n500, function(x) MSE_no_penalty(true_500, x)))
MSE_unpen_corX_500 <- do.call(rbind, lapply(sim2_mod6_mvFMR$optimal$corX$n500, function(x) MSE_no_penalty(true_500, x)))

TPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_indepX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_indepX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$indepX$n1000, function(x) FDR_no_penalty(true_1000, x)))
TPR_unpen_corX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n1000, function(x) TPR_no_penalty(true_1000, x)))
FPR_unpen_corX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n1000, function(x) FPR_no_penalty(true_1000, x)))
FDR_unpen_corX_1000 <- unlist(lapply(sim2_mod6_mvFMR$optimal$corX$n1000, function(x) FDR_no_penalty(true_1000, x)))
MSE_unpen_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_mvFMR$optimal$indepX$n1000, function(x) MSE_no_penalty(true_1000, x)))
MSE_unpen_corX_1000 <- do.call(rbind, lapply(sim2_mod6_mvFMR$optimal$corX$n1000, function(x) MSE_no_penalty(true_1000, x)))

## PGM_mvFMR_LASSO ##
TPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_SCAD ##
TPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## PGM_mvFMR_MCP ##
TPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_PGM_MCP_corX_500 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_PGM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_PGM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_PGM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_PGM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_PGM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod6_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_LASSO ##
TPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_LASSO_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_LASSO_indepX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_LASSO_corX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_LASSO_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_LASSO_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_LASSO_corX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_SCAD ##
TPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_SCAD_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_SCAD_indepX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_SCAD_corX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_SCAD_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_SCAD_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_SCAD_corX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## ADMM_mvFMR_MCP ##
TPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_indepX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) FDR(true_500, x)))
TPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) TPR(true_500, x)))
FPR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FPR(true_500, x)))
FDR_ADMM_MCP_corX_500 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) FDR(true_500, x)))
MSE_ADMM_MCP_indepX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) MSE(true_500, x)))
MSE_ADMM_MCP_corX_500 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) MSE(true_500, x)))

TPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_indepX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) FDR(true_1000, x)))
TPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) TPR(true_1000, x)))
FPR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FPR(true_1000, x)))
FDR_ADMM_MCP_corX_1000 <- unlist(lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) FDR(true_1000, x)))
MSE_ADMM_MCP_indepX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) MSE(true_1000, x)))
MSE_ADMM_MCP_corX_1000 <- do.call(rbind, lapply(sim2_mod6_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) MSE(true_1000, x)))

## Table ##
guide <- data.frame(TPR=rep(NA, 9), FPR=rep(NA, 9), FDR=rep(NA, 9), MSE_pi=rep(NA, 9), MSE_Bk=rep(NA, 9), MSE_sigma=rep(NA, 9))
rownames(guide) <- c("Flexmix", "Oracle", "mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model6_table <- list(indepX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                                 n1000=list(mean=data.frame(guide), sd=data.frame(guide))),
                     corX=list(n500=list(mean=data.frame(guide), sd=data.frame(guide)),
                               n1000=list(mean=data.frame(guide), sd=data.frame(guide))))
model6_table$indepX$n500$mean$TPR <- c(mean(TPR_flexmix_indepX_500), mean(TPR_oracle_indepX_500), mean(TPR_unpen_indepX_500),
                                       mean(TPR_PGM_LASSO_indepX_500), mean(TPR_ADMM_LASSO_indepX_500), mean(TPR_PGM_SCAD_indepX_500),
                                       mean(TPR_ADMM_SCAD_indepX_500), mean(TPR_PGM_MCP_indepX_500), mean(TPR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$mean$FPR <- c(mean(FPR_flexmix_indepX_500), mean(FPR_oracle_indepX_500), mean(FPR_unpen_indepX_500),
                                       mean(FPR_PGM_LASSO_indepX_500), mean(FPR_ADMM_LASSO_indepX_500), mean(FPR_PGM_SCAD_indepX_500),
                                       mean(FPR_ADMM_SCAD_indepX_500), mean(FPR_PGM_MCP_indepX_500), mean(FPR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$mean$FDR <- c(mean(FDR_flexmix_indepX_500), mean(FDR_oracle_indepX_500), mean(FDR_unpen_indepX_500),
                                       mean(FDR_PGM_LASSO_indepX_500), mean(FDR_ADMM_LASSO_indepX_500), mean(FDR_PGM_SCAD_indepX_500),
                                       mean(FDR_ADMM_SCAD_indepX_500), mean(FDR_PGM_MCP_indepX_500), mean(FDR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_500[,1]), mean(MSE_oracle_indepX_500[,1]), mean(MSE_unpen_indepX_500[,1]),
                                          mean(MSE_PGM_LASSO_indepX_500[,1]), mean(MSE_ADMM_LASSO_indepX_500[,1]), mean(MSE_PGM_SCAD_indepX_500[,1]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,1]), mean(MSE_PGM_MCP_indepX_500[,1]), mean(MSE_ADMM_MCP_indepX_500[,1]))
model6_table$indepX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_500[,2]), mean(MSE_oracle_indepX_500[,2]), mean(MSE_unpen_indepX_500[,2]),
                                          mean(MSE_PGM_LASSO_indepX_500[,2]), mean(MSE_ADMM_LASSO_indepX_500[,2]), mean(MSE_PGM_SCAD_indepX_500[,2]),
                                          mean(MSE_ADMM_SCAD_indepX_500[,2]), mean(MSE_PGM_MCP_indepX_500[,2]), mean(MSE_ADMM_MCP_indepX_500[,2]))
model6_table$indepX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_500[,3]), mean(MSE_oracle_indepX_500[,3]), mean(MSE_unpen_indepX_500[,3]),
                                             mean(MSE_PGM_LASSO_indepX_500[,3]), mean(MSE_ADMM_LASSO_indepX_500[,3]), mean(MSE_PGM_SCAD_indepX_500[,3]),
                                             mean(MSE_ADMM_SCAD_indepX_500[,3]), mean(MSE_PGM_MCP_indepX_500[,3]), mean(MSE_ADMM_MCP_indepX_500[,3]))

model6_table$indepX$n500$sd$TPR <- c(sd(TPR_flexmix_indepX_500), sd(TPR_oracle_indepX_500), sd(TPR_unpen_indepX_500),
                                     sd(TPR_PGM_LASSO_indepX_500), sd(TPR_ADMM_LASSO_indepX_500), sd(TPR_PGM_SCAD_indepX_500),
                                     sd(TPR_ADMM_SCAD_indepX_500), sd(TPR_PGM_MCP_indepX_500), sd(TPR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$sd$FPR <- c(sd(FPR_flexmix_indepX_500), sd(FPR_oracle_indepX_500), sd(FPR_unpen_indepX_500),
                                     sd(FPR_PGM_LASSO_indepX_500), sd(FPR_ADMM_LASSO_indepX_500), sd(FPR_PGM_SCAD_indepX_500),
                                     sd(FPR_ADMM_SCAD_indepX_500), sd(FPR_PGM_MCP_indepX_500), sd(FPR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$sd$FDR <- c(sd(FDR_flexmix_indepX_500), sd(FDR_oracle_indepX_500), sd(FDR_unpen_indepX_500),
                                     sd(FDR_PGM_LASSO_indepX_500), sd(FDR_ADMM_LASSO_indepX_500), sd(FDR_PGM_SCAD_indepX_500),
                                     sd(FDR_ADMM_SCAD_indepX_500), sd(FDR_PGM_MCP_indepX_500), sd(FDR_ADMM_MCP_indepX_500))
model6_table$indepX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_500[,1]), sd(MSE_oracle_indepX_500[,1]), sd(MSE_unpen_indepX_500[,1]),
                                        sd(MSE_PGM_LASSO_indepX_500[,1]), sd(MSE_ADMM_LASSO_indepX_500[,1]), sd(MSE_PGM_SCAD_indepX_500[,1]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,1]), sd(MSE_PGM_MCP_indepX_500[,1]), sd(MSE_ADMM_MCP_indepX_500[,1]))
model6_table$indepX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_500[,2]), sd(MSE_oracle_indepX_500[,2]), sd(MSE_unpen_indepX_500[,2]),
                                        sd(MSE_PGM_LASSO_indepX_500[,2]), sd(MSE_ADMM_LASSO_indepX_500[,2]), sd(MSE_PGM_SCAD_indepX_500[,2]),
                                        sd(MSE_ADMM_SCAD_indepX_500[,2]), sd(MSE_PGM_MCP_indepX_500[,2]), sd(MSE_ADMM_MCP_indepX_500[,2]))
model6_table$indepX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_500[,3]), sd(MSE_oracle_indepX_500[,3]), sd(MSE_unpen_indepX_500[,3]),
                                           sd(MSE_PGM_LASSO_indepX_500[,3]), sd(MSE_ADMM_LASSO_indepX_500[,3]), sd(MSE_PGM_SCAD_indepX_500[,3]),
                                           sd(MSE_ADMM_SCAD_indepX_500[,3]), sd(MSE_PGM_MCP_indepX_500[,3]), sd(MSE_ADMM_MCP_indepX_500[,3]))

model6_table$indepX$n1000$mean$TPR <- c(mean(TPR_flexmix_indepX_1000), mean(TPR_oracle_indepX_1000), mean(TPR_unpen_indepX_1000),
                                        mean(TPR_PGM_LASSO_indepX_1000), mean(TPR_ADMM_LASSO_indepX_1000), mean(TPR_PGM_SCAD_indepX_1000),
                                        mean(TPR_ADMM_SCAD_indepX_1000), mean(TPR_PGM_MCP_indepX_1000), mean(TPR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$mean$FPR <- c(mean(FPR_flexmix_indepX_1000), mean(FPR_oracle_indepX_1000), mean(FPR_unpen_indepX_1000),
                                        mean(FPR_PGM_LASSO_indepX_1000), mean(FPR_ADMM_LASSO_indepX_1000), mean(FPR_PGM_SCAD_indepX_1000),
                                        mean(FPR_ADMM_SCAD_indepX_1000), mean(FPR_PGM_MCP_indepX_1000), mean(FPR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$mean$FDR <- c(mean(FDR_flexmix_indepX_1000), mean(FDR_oracle_indepX_1000), mean(FDR_unpen_indepX_1000),
                                        mean(FDR_PGM_LASSO_indepX_1000), mean(FDR_ADMM_LASSO_indepX_1000), mean(FDR_PGM_SCAD_indepX_1000),
                                        mean(FDR_ADMM_SCAD_indepX_1000), mean(FDR_PGM_MCP_indepX_1000), mean(FDR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_indepX_1000[,1]), mean(MSE_oracle_indepX_1000[,1]), mean(MSE_unpen_indepX_1000[,1]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,1]), mean(MSE_ADMM_LASSO_indepX_1000[,1]), mean(MSE_PGM_SCAD_indepX_1000[,1]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,1]), mean(MSE_PGM_MCP_indepX_1000[,1]), mean(MSE_ADMM_MCP_indepX_1000[,1]))
model6_table$indepX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_indepX_1000[,2]), mean(MSE_oracle_indepX_1000[,2]), mean(MSE_unpen_indepX_1000[,2]),
                                           mean(MSE_PGM_LASSO_indepX_1000[,2]), mean(MSE_ADMM_LASSO_indepX_1000[,2]), mean(MSE_PGM_SCAD_indepX_1000[,2]),
                                           mean(MSE_ADMM_SCAD_indepX_1000[,2]), mean(MSE_PGM_MCP_indepX_1000[,2]), mean(MSE_ADMM_MCP_indepX_1000[,2]))
model6_table$indepX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_indepX_1000[,3]), mean(MSE_oracle_indepX_1000[,3]), mean(MSE_unpen_indepX_1000[,3]),
                                              mean(MSE_PGM_LASSO_indepX_1000[,3]), mean(MSE_ADMM_LASSO_indepX_1000[,3]), mean(MSE_PGM_SCAD_indepX_1000[,3]),
                                              mean(MSE_ADMM_SCAD_indepX_1000[,3]), mean(MSE_PGM_MCP_indepX_1000[,3]), mean(MSE_ADMM_MCP_indepX_1000[,3]))

model6_table$indepX$n1000$sd$TPR <- c(sd(TPR_flexmix_indepX_1000), sd(TPR_oracle_indepX_1000), sd(TPR_unpen_indepX_1000),
                                      sd(TPR_PGM_LASSO_indepX_1000), sd(TPR_ADMM_LASSO_indepX_1000), sd(TPR_PGM_SCAD_indepX_1000),
                                      sd(TPR_ADMM_SCAD_indepX_1000), sd(TPR_PGM_MCP_indepX_1000), sd(TPR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$sd$FPR <- c(sd(FPR_flexmix_indepX_1000), sd(FPR_oracle_indepX_1000), sd(FPR_unpen_indepX_1000),
                                      sd(FPR_PGM_LASSO_indepX_1000), sd(FPR_ADMM_LASSO_indepX_1000), sd(FPR_PGM_SCAD_indepX_1000),
                                      sd(FPR_ADMM_SCAD_indepX_1000), sd(FPR_PGM_MCP_indepX_1000), sd(FPR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$sd$FDR <- c(sd(FDR_flexmix_indepX_1000), sd(FDR_oracle_indepX_1000), sd(FDR_unpen_indepX_1000),
                                      sd(FDR_PGM_LASSO_indepX_1000), sd(FDR_ADMM_LASSO_indepX_1000), sd(FDR_PGM_SCAD_indepX_1000),
                                      sd(FDR_ADMM_SCAD_indepX_1000), sd(FDR_PGM_MCP_indepX_1000), sd(FDR_ADMM_MCP_indepX_1000))
model6_table$indepX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_indepX_1000[,1]), sd(MSE_oracle_indepX_1000[,1]), sd(MSE_unpen_indepX_1000[,1]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,1]), sd(MSE_ADMM_LASSO_indepX_1000[,1]), sd(MSE_PGM_SCAD_indepX_1000[,1]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,1]), sd(MSE_PGM_MCP_indepX_1000[,1]), sd(MSE_ADMM_MCP_indepX_1000[,1]))
model6_table$indepX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_indepX_1000[,2]), sd(MSE_oracle_indepX_1000[,2]), sd(MSE_unpen_indepX_1000[,2]),
                                         sd(MSE_PGM_LASSO_indepX_1000[,2]), sd(MSE_ADMM_LASSO_indepX_1000[,2]), sd(MSE_PGM_SCAD_indepX_1000[,2]),
                                         sd(MSE_ADMM_SCAD_indepX_1000[,2]), sd(MSE_PGM_MCP_indepX_1000[,2]), sd(MSE_ADMM_MCP_indepX_1000[,2]))
model6_table$indepX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_indepX_1000[,3]), sd(MSE_oracle_indepX_1000[,3]), sd(MSE_unpen_indepX_1000[,3]),
                                            sd(MSE_PGM_LASSO_indepX_1000[,3]), sd(MSE_ADMM_LASSO_indepX_1000[,3]), sd(MSE_PGM_SCAD_indepX_1000[,3]),
                                            sd(MSE_ADMM_SCAD_indepX_1000[,3]), sd(MSE_PGM_MCP_indepX_1000[,3]), sd(MSE_ADMM_MCP_indepX_1000[,3]))

model6_table$corX$n500$mean$TPR <- c(mean(TPR_flexmix_corX_500), mean(TPR_oracle_corX_500), mean(TPR_unpen_corX_500),
                                     mean(TPR_PGM_LASSO_corX_500), mean(TPR_ADMM_LASSO_corX_500), mean(TPR_PGM_SCAD_corX_500),
                                     mean(TPR_ADMM_SCAD_corX_500), mean(TPR_PGM_MCP_corX_500), mean(TPR_ADMM_MCP_corX_500))
model6_table$corX$n500$mean$FPR <- c(mean(FPR_flexmix_corX_500), mean(FPR_oracle_corX_500), mean(FPR_unpen_corX_500),
                                     mean(FPR_PGM_LASSO_corX_500), mean(FPR_ADMM_LASSO_corX_500), mean(FPR_PGM_SCAD_corX_500),
                                     mean(FPR_ADMM_SCAD_corX_500), mean(FPR_PGM_MCP_corX_500), mean(FPR_ADMM_MCP_corX_500))
model6_table$corX$n500$mean$FDR <- c(mean(FDR_flexmix_corX_500), mean(FDR_oracle_corX_500), mean(FDR_unpen_corX_500),
                                     mean(FDR_PGM_LASSO_corX_500), mean(FDR_ADMM_LASSO_corX_500), mean(FDR_PGM_SCAD_corX_500),
                                     mean(FDR_ADMM_SCAD_corX_500), mean(FDR_PGM_MCP_corX_500), mean(FDR_ADMM_MCP_corX_500))
model6_table$corX$n500$mean$MSE_pi <- c(mean(MSE_flexmix_corX_500[,1]), mean(MSE_oracle_corX_500[,1]), mean(MSE_unpen_corX_500[,1]),
                                        mean(MSE_PGM_LASSO_corX_500[,1]), mean(MSE_ADMM_LASSO_corX_500[,1]), mean(MSE_PGM_SCAD_corX_500[,1]),
                                        mean(MSE_ADMM_SCAD_corX_500[,1]), mean(MSE_PGM_MCP_corX_500[,1]), mean(MSE_ADMM_MCP_corX_500[,1]))
model6_table$corX$n500$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_500[,2]), mean(MSE_oracle_corX_500[,2]), mean(MSE_unpen_corX_500[,2]),
                                        mean(MSE_PGM_LASSO_corX_500[,2]), mean(MSE_ADMM_LASSO_corX_500[,2]), mean(MSE_PGM_SCAD_corX_500[,2]),
                                        mean(MSE_ADMM_SCAD_corX_500[,2]), mean(MSE_PGM_MCP_corX_500[,2]), mean(MSE_ADMM_MCP_corX_500[,2]))
model6_table$corX$n500$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_500[,3]), mean(MSE_oracle_corX_500[,3]), mean(MSE_unpen_corX_500[,3]),
                                           mean(MSE_PGM_LASSO_corX_500[,3]), mean(MSE_ADMM_LASSO_corX_500[,3]), mean(MSE_PGM_SCAD_corX_500[,3]),
                                           mean(MSE_ADMM_SCAD_corX_500[,3]), mean(MSE_PGM_MCP_corX_500[,3]), mean(MSE_ADMM_MCP_corX_500[,3]))

model6_table$corX$n500$sd$TPR <- c(sd(TPR_flexmix_corX_500), sd(TPR_oracle_corX_500), sd(TPR_unpen_corX_500),
                                   sd(TPR_PGM_LASSO_corX_500), sd(TPR_ADMM_LASSO_corX_500), sd(TPR_PGM_SCAD_corX_500),
                                   sd(TPR_ADMM_SCAD_corX_500), sd(TPR_PGM_MCP_corX_500), sd(TPR_ADMM_MCP_corX_500))
model6_table$corX$n500$sd$FPR <- c(sd(FPR_flexmix_corX_500), sd(FPR_oracle_corX_500), sd(FPR_unpen_corX_500),
                                   sd(FPR_PGM_LASSO_corX_500), sd(FPR_ADMM_LASSO_corX_500), sd(FPR_PGM_SCAD_corX_500),
                                   sd(FPR_ADMM_SCAD_corX_500), sd(FPR_PGM_MCP_corX_500), sd(FPR_ADMM_MCP_corX_500))
model6_table$corX$n500$sd$FDR <- c(sd(FDR_flexmix_corX_500), sd(FDR_oracle_corX_500), sd(FDR_unpen_corX_500),
                                   sd(FDR_PGM_LASSO_corX_500), sd(FDR_ADMM_LASSO_corX_500), sd(FDR_PGM_SCAD_corX_500),
                                   sd(FDR_ADMM_SCAD_corX_500), sd(FDR_PGM_MCP_corX_500), sd(FDR_ADMM_MCP_corX_500))
model6_table$corX$n500$sd$MSE_pi <- c(sd(MSE_flexmix_corX_500[,1]), sd(MSE_oracle_corX_500[,1]), sd(MSE_unpen_corX_500[,1]),
                                      sd(MSE_PGM_LASSO_corX_500[,1]), sd(MSE_ADMM_LASSO_corX_500[,1]), sd(MSE_PGM_SCAD_corX_500[,1]),
                                      sd(MSE_ADMM_SCAD_corX_500[,1]), sd(MSE_PGM_MCP_corX_500[,1]), sd(MSE_ADMM_MCP_corX_500[,1]))
model6_table$corX$n500$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_500[,2]), sd(MSE_oracle_corX_500[,2]), sd(MSE_unpen_corX_500[,2]),
                                      sd(MSE_PGM_LASSO_corX_500[,2]), sd(MSE_ADMM_LASSO_corX_500[,2]), sd(MSE_PGM_SCAD_corX_500[,2]),
                                      sd(MSE_ADMM_SCAD_corX_500[,2]), sd(MSE_PGM_MCP_corX_500[,2]), sd(MSE_ADMM_MCP_corX_500[,2]))
model6_table$corX$n500$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_500[,3]), sd(MSE_oracle_corX_500[,3]), sd(MSE_unpen_corX_500[,3]),
                                         sd(MSE_PGM_LASSO_corX_500[,3]), sd(MSE_ADMM_LASSO_corX_500[,3]), sd(MSE_PGM_SCAD_corX_500[,3]),
                                         sd(MSE_ADMM_SCAD_corX_500[,3]), sd(MSE_PGM_MCP_corX_500[,3]), sd(MSE_ADMM_MCP_corX_500[,3]))

model6_table$corX$n1000$mean$TPR <- c(mean(TPR_flexmix_corX_1000), mean(TPR_oracle_corX_1000), mean(TPR_unpen_corX_1000),
                                      mean(TPR_PGM_LASSO_corX_1000), mean(TPR_ADMM_LASSO_corX_1000), mean(TPR_PGM_SCAD_corX_1000),
                                      mean(TPR_ADMM_SCAD_corX_1000), mean(TPR_PGM_MCP_corX_1000), mean(TPR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$mean$FPR <- c(mean(FPR_flexmix_corX_1000), mean(FPR_oracle_corX_1000), mean(FPR_unpen_corX_1000),
                                      mean(FPR_PGM_LASSO_corX_1000), mean(FPR_ADMM_LASSO_corX_1000), mean(FPR_PGM_SCAD_corX_1000),
                                      mean(FPR_ADMM_SCAD_corX_1000), mean(FPR_PGM_MCP_corX_1000), mean(FPR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$mean$FDR <- c(mean(FDR_flexmix_corX_1000), mean(FDR_oracle_corX_1000), mean(FDR_unpen_corX_1000),
                                      mean(FDR_PGM_LASSO_corX_1000), mean(FDR_ADMM_LASSO_corX_1000), mean(FDR_PGM_SCAD_corX_1000),
                                      mean(FDR_ADMM_SCAD_corX_1000), mean(FDR_PGM_MCP_corX_1000), mean(FDR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$mean$MSE_pi <- c(mean(MSE_flexmix_corX_1000[,1]), mean(MSE_oracle_corX_1000[,1]), mean(MSE_unpen_corX_1000[,1]),
                                         mean(MSE_PGM_LASSO_corX_1000[,1]), mean(MSE_ADMM_LASSO_corX_1000[,1]), mean(MSE_PGM_SCAD_corX_1000[,1]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,1]), mean(MSE_PGM_MCP_corX_1000[,1]), mean(MSE_ADMM_MCP_corX_1000[,1]))
model6_table$corX$n1000$mean$MSE_Bk <- c(mean(MSE_flexmix_corX_1000[,2]), mean(MSE_oracle_corX_1000[,2]), mean(MSE_unpen_corX_1000[,2]),
                                         mean(MSE_PGM_LASSO_corX_1000[,2]), mean(MSE_ADMM_LASSO_corX_1000[,2]), mean(MSE_PGM_SCAD_corX_1000[,2]),
                                         mean(MSE_ADMM_SCAD_corX_1000[,2]), mean(MSE_PGM_MCP_corX_1000[,2]), mean(MSE_ADMM_MCP_corX_1000[,2]))
model6_table$corX$n1000$mean$MSE_sigma <- c(mean(MSE_flexmix_corX_1000[,3]), mean(MSE_oracle_corX_1000[,3]), mean(MSE_unpen_corX_1000[,3]),
                                            mean(MSE_PGM_LASSO_corX_1000[,3]), mean(MSE_ADMM_LASSO_corX_1000[,3]), mean(MSE_PGM_SCAD_corX_1000[,3]),
                                            mean(MSE_ADMM_SCAD_corX_1000[,3]), mean(MSE_PGM_MCP_corX_1000[,3]), mean(MSE_ADMM_MCP_corX_1000[,3]))

model6_table$corX$n1000$sd$TPR <- c(sd(TPR_flexmix_corX_1000), sd(TPR_oracle_corX_1000), sd(TPR_unpen_corX_1000),
                                    sd(TPR_PGM_LASSO_corX_1000), sd(TPR_ADMM_LASSO_corX_1000), sd(TPR_PGM_SCAD_corX_1000),
                                    sd(TPR_ADMM_SCAD_corX_1000), sd(TPR_PGM_MCP_corX_1000), sd(TPR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$sd$FPR <- c(sd(FPR_flexmix_corX_1000), sd(FPR_oracle_corX_1000), sd(FPR_unpen_corX_1000),
                                    sd(FPR_PGM_LASSO_corX_1000), sd(FPR_ADMM_LASSO_corX_1000), sd(FPR_PGM_SCAD_corX_1000),
                                    sd(FPR_ADMM_SCAD_corX_1000), sd(FPR_PGM_MCP_corX_1000), sd(FPR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$sd$FDR <- c(sd(FDR_flexmix_corX_1000), sd(FDR_oracle_corX_1000), sd(FDR_unpen_corX_1000),
                                    sd(FDR_PGM_LASSO_corX_1000), sd(FDR_ADMM_LASSO_corX_1000), sd(FDR_PGM_SCAD_corX_1000),
                                    sd(FDR_ADMM_SCAD_corX_1000), sd(FDR_PGM_MCP_corX_1000), sd(FDR_ADMM_MCP_corX_1000))
model6_table$corX$n1000$sd$MSE_pi <- c(sd(MSE_flexmix_corX_1000[,1]), sd(MSE_oracle_corX_1000[,1]), sd(MSE_unpen_corX_1000[,1]),
                                       sd(MSE_PGM_LASSO_corX_1000[,1]), sd(MSE_ADMM_LASSO_corX_1000[,1]), sd(MSE_PGM_SCAD_corX_1000[,1]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,1]), sd(MSE_PGM_MCP_corX_1000[,1]), sd(MSE_ADMM_MCP_corX_1000[,1]))
model6_table$corX$n1000$sd$MSE_Bk <- c(sd(MSE_flexmix_corX_1000[,2]), sd(MSE_oracle_corX_1000[,2]), sd(MSE_unpen_corX_1000[,2]),
                                       sd(MSE_PGM_LASSO_corX_1000[,2]), sd(MSE_ADMM_LASSO_corX_1000[,2]), sd(MSE_PGM_SCAD_corX_1000[,2]),
                                       sd(MSE_ADMM_SCAD_corX_1000[,2]), sd(MSE_PGM_MCP_corX_1000[,2]), sd(MSE_ADMM_MCP_corX_1000[,2]))
model6_table$corX$n1000$sd$MSE_sigma <- c(sd(MSE_flexmix_corX_1000[,3]), sd(MSE_oracle_corX_1000[,3]), sd(MSE_unpen_corX_1000[,3]),
                                          sd(MSE_PGM_LASSO_corX_1000[,3]), sd(MSE_ADMM_LASSO_corX_1000[,3]), sd(MSE_PGM_SCAD_corX_1000[,3]),
                                          sd(MSE_ADMM_SCAD_corX_1000[,3]), sd(MSE_PGM_MCP_corX_1000[,3]), sd(MSE_ADMM_MCP_corX_1000[,3]))

mod6_table <- list(indepX=list(n500=data.frame(guide), n1000=data.frame(guide)),
                   corX=list(n500=data.frame(guide), n1000=data.frame(guide)))
mod6_table$indepX$n500$TPR <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$TPR, model6_table$indepX$n500$sd$TPR)
mod6_table$indepX$n500$FPR <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$FPR, model6_table$indepX$n500$sd$FPR)
mod6_table$indepX$n500$FDR <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$FDR, model6_table$indepX$n500$sd$FDR)
mod6_table$indepX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$MSE_pi, model6_table$indepX$n500$sd$MSE_pi)
mod6_table$indepX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$MSE_Bk, model6_table$indepX$n500$sd$MSE_Bk)
mod6_table$indepX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model6_table$indepX$n500$mean$MSE_sigma, model6_table$indepX$n500$sd$MSE_sigma)

mod6_table$indepX$n1000$TPR <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$TPR, model6_table$indepX$n1000$sd$TPR)
mod6_table$indepX$n1000$FPR <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$FPR, model6_table$indepX$n1000$sd$FPR)
mod6_table$indepX$n1000$FDR <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$FDR, model6_table$indepX$n1000$sd$FDR)
mod6_table$indepX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$MSE_pi, model6_table$indepX$n1000$sd$MSE_pi)
mod6_table$indepX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$MSE_Bk, model6_table$indepX$n1000$sd$MSE_Bk)
mod6_table$indepX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model6_table$indepX$n1000$mean$MSE_sigma, model6_table$indepX$n1000$sd$MSE_sigma)

mod6_table$corX$n500$TPR <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$TPR, model6_table$corX$n500$sd$TPR)
mod6_table$corX$n500$FPR <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$FPR, model6_table$corX$n500$sd$FPR)
mod6_table$corX$n500$FDR <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$FDR, model6_table$corX$n500$sd$FDR)
mod6_table$corX$n500$MSE_pi <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$MSE_pi, model6_table$corX$n500$sd$MSE_pi)
mod6_table$corX$n500$MSE_Bk <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$MSE_Bk, model6_table$corX$n500$sd$MSE_Bk)
mod6_table$corX$n500$MSE_sigma <- sprintf("%.4f(%.4f)", model6_table$corX$n500$mean$MSE_sigma, model6_table$corX$n500$sd$MSE_sigma)

mod6_table$corX$n1000$TPR <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$TPR, model6_table$corX$n1000$sd$TPR)
mod6_table$corX$n1000$FPR <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$FPR, model6_table$corX$n1000$sd$FPR)
mod6_table$corX$n1000$FDR <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$FDR, model6_table$corX$n1000$sd$FDR)
mod6_table$corX$n1000$MSE_pi <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$MSE_pi, model6_table$corX$n1000$sd$MSE_pi)
mod6_table$corX$n1000$MSE_Bk <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$MSE_Bk, model6_table$corX$n1000$sd$MSE_Bk)
mod6_table$corX$n1000$MSE_sigma <- sprintf("%.4f(%.4f)", model6_table$corX$n1000$mean$MSE_sigma, model6_table$corX$n1000$sd$MSE_sigma)
# save(mod6_table, file="./Results/Section5.2/Model6/mod6_table.rda")