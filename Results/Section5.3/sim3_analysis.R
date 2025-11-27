###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####           R-code for reproducing Tables 10-11 in Section 5.3.             ####
###################################################################################

# The following R-code contains functions for calculating predictive log-likelihood,
# which is the twice negative log-likelihood.

# Each ".rda" file contains the modified BIC values, optimal estimators('pi', 'Bk', 'sigma'), 
# and the estimators of mixing proportion 'w'.

### Data ###
source("./Data/simulation_data.R")
source("./Data/simulation_seed_number.R")

### Functions ###
source("./Scripts/Functions/functions.R")
source("./Scripts/Functions/simulation_summary_functions.R")

### Model 7 ###
## results ##
load("./Results/Section5.3/Model7/sim3_mod7_mvFMR.rda")
load("./Results/Section5.3/Model7/sim3_mod7_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.3/Model7/sim3_mod7_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.3/Model7/sim3_mod7_PGM_mvFMR_MCP.rda")
load("./Results/Section5.3/Model7/sim3_mod7_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.3/Model7/sim3_mod7_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.3/Model7/sim3_mod7_ADMM_mvFMR_MCP.rda")

## mvFMR ##
dens_unpen_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                sim3_mod7_mvFMR$optimal$indepX$n500[[x]]))
pred_ll_unpen_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_mvFMR$optimal$indepX$n500, 
                                                                    dens_unpen_indepX_500, x, 500))

dens_unpen_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                 sim3_mod7_mvFMR$optimal$indepX$n1000[[x]]))
pred_ll_unpen_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_mvFMR$optimal$indepX$n1000, 
                                                                     dens_unpen_indepX_1000, x, 1000))

dens_unpen_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                              sim3_mod7_mvFMR$optimal$corX$n500[[x]]))
pred_ll_unpen_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_mvFMR$optimal$corX$n500, 
                                                                  dens_unpen_corX_500, x, 500))

dens_unpen_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                               sim3_mod7_mvFMR$optimal$corX$n1000[[x]]))
pred_ll_unpen_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_mvFMR$optimal$corX$n1000, 
                                                                   dens_unpen_corX_1000, x, 1000))

K_unpen_indepX_500 <- sapply(sim3_mod7_mvFMR$optimal$indepX$n500, function(x) length(x$Bk))
K_unpen_indepX_1000 <- sapply(sim3_mod7_mvFMR$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_unpen_corX_500 <- sapply(sim3_mod7_mvFMR$optimal$corX$n500, function(x) length(x$Bk))
K_unpen_corX_1000 <- sapply(sim3_mod7_mvFMR$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_LASSO ##
dens_PGM_mvFMR_LASSO_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                          sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_LASSO_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n500, 
                                                                              dens_PGM_LASSO_indepX_500, x, 500))

dens_PGM_mvFMR_LASSO_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                           sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_LASSO_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n1000,
                                                                               dens_PGM_LASSO_indepX_1000, x, 1000))

dens_PGM_mvFMR_LASSO_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                        sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_LASSO_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n500, 
                                                                            dens_PGM_LASSO_corX_500, x, 500))

dens_PGM_mvFMR_LASSO_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                         sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_LASSO_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n1000, 
                                                                             dens_PGM_LASSO_corX_1000, x, 1000))

K_PGM_mvFMR_LASSO_indepX_500 <- sapply(sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_LASSO_indepX_1000 <- sapply(sim3_mod7_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_LASSO_corX_500 <- sapply(sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_LASSO_corX_1000 <- sapply(sim3_mod7_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_SCAD ##
dens_PGM_mvFMR_SCAD_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                         sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_SCAD_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n500, 
                                                                             dens_PGM_SCAD_indepX_500, x, 500))

dens_PGM_mvFMR_SCAD_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                          sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_SCAD_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n1000,
                                                                              dens_PGM_SCAD_indepX_1000, x, 1000))

dens_PGM_mvFMR_SCAD_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                       sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_SCAD_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n500, 
                                                                           dens_PGM_SCAD_corX_500, x, 500))

dens_PGM_mvFMR_SCAD_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                        sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_SCAD_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n1000, 
                                                                            dens_PGM_SCAD_corX_1000, x, 1000))

K_PGM_mvFMR_SCAD_indepX_500 <- sapply(sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_SCAD_indepX_1000 <- sapply(sim3_mod7_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_SCAD_corX_500 <- sapply(sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_SCAD_corX_1000 <- sapply(sim3_mod7_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_MCP ##
dens_PGM_mvFMR_MCP_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                        sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_MCP_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n500, 
                                                                            dens_PGM_MCP_indepX_500, x, 500))

dens_PGM_mvFMR_MCP_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                         sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_MCP_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n1000,
                                                                             dens_PGM_MCP_indepX_1000, x, 1000))

dens_PGM_mvFMR_MCP_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                      sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_MCP_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n500, 
                                                                          dens_PGM_MCP_corX_500, x, 500))

dens_PGM_mvFMR_MCP_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                       sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_MCP_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n1000, 
                                                                           dens_PGM_MCP_corX_1000, x, 1000))

K_PGM_mvFMR_MCP_indepX_500 <- sapply(sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_MCP_indepX_1000 <- sapply(sim3_mod7_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_MCP_corX_500 <- sapply(sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_MCP_corX_1000 <- sapply(sim3_mod7_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_LASSO ##
dens_ADMM_mvFMR_LASSO_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                           sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_LASSO_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n500, 
                                                                               dens_ADMM_LASSO_indepX_500, x, 500))

dens_ADMM_mvFMR_LASSO_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                            sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_LASSO_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n1000,
                                                                                dens_ADMM_LASSO_indepX_1000, x, 1000))

dens_ADMM_mvFMR_LASSO_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                         sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_LASSO_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n500, 
                                                                             dens_ADMM_LASSO_corX_500, x, 500))

dens_ADMM_mvFMR_LASSO_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                          sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_LASSO_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n1000, 
                                                                              dens_ADMM_LASSO_corX_1000, x, 1000))

K_ADMM_mvFMR_LASSO_indepX_500 <- sapply(sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_LASSO_indepX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_LASSO_corX_500 <- sapply(sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_LASSO_corX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_SCAD ##
dens_ADMM_mvFMR_SCAD_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                          sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_SCAD_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n500, 
                                                                              dens_ADMM_SCAD_indepX_500, x, 500))

dens_ADMM_mvFMR_SCAD_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                           sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_SCAD_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n1000,
                                                                               dens_ADMM_SCAD_indepX_1000, x, 1000))

dens_ADMM_mvFMR_SCAD_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                        sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_SCAD_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n500, 
                                                                            dens_ADMM_SCAD_corX_500, x, 500))

dens_ADMM_mvFMR_SCAD_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                         sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_SCAD_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n1000, 
                                                                             dens_ADMM_SCAD_corX_1000, x, 1000))

K_ADMM_mvFMR_SCAD_indepX_500 <- sapply(sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_SCAD_indepX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_SCAD_corX_500 <- sapply(sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_SCAD_corX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_MCP ##
dens_ADMM_mvFMR_MCP_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, TRUE), 
                                                                         sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_MCP_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n500, 
                                                                             dens_ADMM_MCP_indepX_500, x, 500))

dens_ADMM_mvFMR_MCP_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, TRUE), 
                                                                          sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_MCP_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n1000,
                                                                              dens_ADMM_MCP_indepX_1000, x, 1000))

dens_ADMM_mvFMR_MCP_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 500, FALSE), 
                                                                       sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_MCP_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n500, 
                                                                           dens_ADMM_MCP_corX_500, x, 500))

dens_ADMM_mvFMR_MCP_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.1(seed_number_5.3.1[x], 1000, FALSE), 
                                                                        sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_MCP_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n1000, 
                                                                            dens_ADMM_MCP_corX_1000, x, 1000))

K_ADMM_mvFMR_MCP_indepX_500 <- sapply(sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_MCP_indepX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_MCP_corX_500 <- sapply(sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_MCP_corX_1000 <- sapply(sim3_mod7_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) length(x$Bk))

## Table ##
model7_table_mean <- list(indepX=data.frame(n500=vector(length=7), n1000=vector(length=7)), 
                          corX=data.frame(n500=vector(length=7), n1000=vector(length=7)))

rownames(model7_table_mean$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                        "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model7_table_mean$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                      "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model7_table_mean$indepX$n500 <- c(mean(pred_ll_UNPEN_indepX_500), mean(pred_ll_PGM_LASSO_indepX_500), 
                                   mean(pred_ll_ADMM_LASSO_indepX_500), mean(pred_ll_PGM_SCAD_indepX_500),
                                   mean(pred_ll_ADMM_SCAD_indepX_500), mean(pred_ll_PGM_MCP_indepX_500),
                                   mean(pred_ll_ADMM_MCP_indepX_500))
model7_table_mean$indepX$n1000 <- c(mean(pred_ll_UNPEN_indepX_1000), mean(pred_ll_PGM_LASSO_indepX_1000),
                                    mean(pred_ll_ADMM_LASSO_indepX_1000), mean(pred_ll_PGM_SCAD_indepX_1000),
                                    mean(pred_ll_ADMM_SCAD_indepX_1000), mean(pred_ll_PGM_MCP_indepX_1000), 
                                    mean(pred_ll_ADMM_MCP_indepX_1000))
model7_table_mean$corX$n500 <- c(mean(pred_ll_UNPEN_corX_500), mean(pred_ll_PGM_LASSO_corX_500), 
                                 mean(pred_ll_ADMM_LASSO_corX_500), mean(pred_ll_PGM_SCAD_corX_500),
                                 mean(pred_ll_ADMM_SCAD_corX_500), mean(pred_ll_PGM_MCP_corX_500),
                                 mean(pred_ll_ADMM_MCP_corX_500))
model7_table_mean$corX$n1000 <- c(mean(pred_ll_UNPEN_corX_1000), mean(pred_ll_PGM_LASSO_corX_1000), 
                                  mean(pred_ll_ADMM_LASSO_corX_1000), mean(pred_ll_PGM_SCAD_corX_1000),
                                  mean(pred_ll_ADMM_SCAD_corX_1000), mean(pred_ll_PGM_MCP_corX_1000),
                                  mean(pred_ll_ADMM_MCP_corX_1000))

model7_table_sd <- list(indepX=data.frame(n500=vector(length=7), n1000=vector(length=7)), 
                        corX=data.frame(n500=vector(length=7), n1000=vector(length=7)))
rownames(model7_table_sd$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                      "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model7_table_sd$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                    "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model7_table_sd$indepX$n500 <- c(sd(pred_ll_UNPEN_indepX_500), sd(pred_ll_PGM_LASSO_indepX_500), 
                                 sd(pred_ll_ADMM_LASSO_indepX_500), sd(pred_ll_PGM_SCAD_indepX_500),
                                 sd(pred_ll_ADMM_SCAD_indepX_500), sd(pred_ll_PGM_MCP_indepX_500),
                                 sd(pred_ll_ADMM_MCP_indepX_500))
model7_table_sd$indepX$n1000 <- c(sd(pred_ll_UNPEN_indepX_1000), sd(pred_ll_PGM_LASSO_indepX_1000), 
                                  sd(pred_ll_ADMM_LASSO_indepX_1000), sd(pred_ll_PGM_SCAD_indepX_1000),
                                  sd(pred_ll_ADMM_SCAD_indepX_1000), sd(pred_ll_PGM_MCP_indepX_1000), 
                                  sd(pred_ll_ADMM_MCP_indepX_1000))
model7_table_sd$corX$n500 <- c(sd(pred_ll_UNPEN_corX_500), sd(pred_ll_PGM_LASSO_corX_500), 
                               sd(pred_ll_ADMM_LASSO_corX_500), sd(pred_ll_PGM_SCAD_corX_500),
                               sd(pred_ll_ADMM_SCAD_corX_500), sd(pred_ll_PGM_MCP_corX_500), 
                               sd(pred_ll_ADMM_MCP_corX_500))
model7_table_sd$corX$n1000 <- c(sd(pred_ll_UNPEN_corX_1000), sd(pred_ll_PGM_LASSO_corX_1000), 
                                sd(pred_ll_ADMM_LASSO_corX_1000), sd(pred_ll_PGM_SCAD_corX_1000),
                                sd(pred_ll_ADMM_SCAD_corX_1000), sd(pred_ll_PGM_MCP_corX_1000),
                                sd(pred_ll_ADMM_MCP_corX_1000))

model7_table <- list(indepX=data.frame(n500=rep(NA,7), n1000=rep(NA,7)), 
                     corX=data.frame(n500=rep(NA,7), n1000=rep(NA,7)))
rownames(model7_table$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM",
                                   "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model7_table$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                 "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model7_table$indepX$n500 <- sprintf("%.4f(%.4f)", model7_table_mean$indepX$n500, model7_table_sd$indepX$n500)
model7_table$indepX$n1000 <- sprintf("%.4f(%.4f)", model7_table_mean$indepX$n1000, model7_table_sd$indepX$n1000)
model7_table$corX$n500 <- sprintf("%.4f(%.4f)", model7_table_mean$corX$n500, model7_table_sd$corX$n500)
model7_table$corX$n1000 <- sprintf("%.4f(%.4f)", model7_table_mean$corX$n1000, model7_table_sd$corX$n1000)
# mod7_table <- model7_table
# save(mod7_table, file="./Results/Section5.3/Model7/mod7_table.rda")

mod7_K_table <- list(indepX=list(n500=data.frame(mvFMR=K_unpen_indepX_500, 
                                                 PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_indepX_500,
                                                 ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_indepX_500,
                                                 PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_indepX_500,
                                                 ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_indepX_500,
                                                 PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_indepX_500,
                                                 ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_indepX_500),
                                 n1000=data.frame(mvFMR=K_unpen_indepX_1000, 
                                                  PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_indepX_1000,
                                                  ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_indepX_1000,
                                                  PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_indepX_1000,
                                                  ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_indepX_1000,
                                                  PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_indepX_1000,
                                                  ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_indepX_1000)),
                     corX=list(n500=data.frame(mvFMR=K_unpen_corX_500, 
                                               PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_corX_500,
                                               ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_corX_500,
                                               PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_corX_500,
                                               ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_corX_500,
                                               PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_corX_500,
                                               ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_corX_500),
                               n1000=data.frame(mvFMR=K_unpen_corX_1000, 
                                                PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_corX_1000,
                                                ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_corX_1000,
                                                PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_corX_1000,
                                                ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_corX_1000,
                                                PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_corX_1000,
                                                ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_corX_1000)))
# save(mod7_K_table, file="./Results/Section5.3/Model7/mod7_K_table.rda")


### Model 8 ###
## results ##
load("./Results/Section5.3/Model8/sim3_mod8_mvFMR.rda")
load("./Results/Section5.3/Model8/sim3_mod8_PGM_mvFMR_LASSO.rda")
load("./Results/Section5.3/Model8/sim3_mod8_PGM_mvFMR_SCAD.rda")
load("./Results/Section5.3/Model8/sim3_mod8_PGM_mvFMR_MCP.rda")
load("./Results/Section5.3/Model8/sim3_mod8_ADMM_mvFMR_LASSO.rda")
load("./Results/Section5.3/Model8/sim3_mod8_ADMM_mvFMR_SCAD.rda")
load("./Results/Section5.3/Model8/sim3_mod8_ADMM_mvFMR_MCP.rda")

## mvFMR ##
dens_unpen_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                sim3_mod8_mvFMR$optimal$indepX$n500[[x]]))
pred_ll_unpen_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_mvFMR$optimal$indepX$n500, 
                                                                    dens_unpen_indepX_500, x, 500))

dens_unpen_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                 sim3_mod8_mvFMR$optimal$indepX$n1000[[x]]))
pred_ll_unpen_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_mvFMR$optimal$indepX$n1000, 
                                                                     dens_unpen_indepX_1000, x, 1000))

dens_unpen_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                              sim3_mod8_mvFMR$optimal$corX$n500[[x]]))
pred_ll_unpen_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_mvFMR$optimal$corX$n500, 
                                                                  dens_unpen_corX_500, x, 500))

dens_unpen_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                               sim3_mod8_mvFMR$optimal$corX$n1000[[x]]))
pred_ll_unpen_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_mvFMR$optimal$corX$n1000, 
                                                                   dens_unpen_corX_1000, x, 1000))

K_unpen_indepX_500 <- sapply(sim3_mod8_mvFMR$optimal$indepX$n500, function(x) length(x$Bk))
K_unpen_indepX_1000 <- sapply(sim3_mod8_mvFMR$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_unpen_corX_500 <- sapply(sim3_mod8_mvFMR$optimal$corX$n500, function(x) length(x$Bk))
K_unpen_corX_1000 <- sapply(sim3_mod8_mvFMR$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_LASSO ##
dens_PGM_mvFMR_LASSO_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                          sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_LASSO_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n500, 
                                                                              dens_PGM_LASSO_indepX_500, x, 500))

dens_PGM_mvFMR_LASSO_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                           sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_LASSO_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n1000,
                                                                               dens_PGM_LASSO_indepX_1000, x, 1000))

dens_PGM_mvFMR_LASSO_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                        sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_LASSO_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n500, 
                                                                            dens_PGM_LASSO_corX_500, x, 500))

dens_PGM_mvFMR_LASSO_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                         sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_LASSO_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n1000, 
                                                                             dens_PGM_LASSO_corX_1000, x, 1000))

K_PGM_mvFMR_LASSO_indepX_500 <- sapply(sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_LASSO_indepX_1000 <- sapply(sim3_mod8_PGM_mvFMR_LASSO$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_LASSO_corX_500 <- sapply(sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_LASSO_corX_1000 <- sapply(sim3_mod8_PGM_mvFMR_LASSO$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_SCAD ##
dens_PGM_mvFMR_SCAD_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                         sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_SCAD_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n500, 
                                                                             dens_PGM_SCAD_indepX_500, x, 500))

dens_PGM_mvFMR_SCAD_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                          sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_SCAD_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n1000,
                                                                              dens_PGM_SCAD_indepX_1000, x, 1000))

dens_PGM_mvFMR_SCAD_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                       sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_SCAD_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n500, 
                                                                           dens_PGM_SCAD_corX_500, x, 500))

dens_PGM_mvFMR_SCAD_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                        sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_SCAD_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n1000, 
                                                                            dens_PGM_SCAD_corX_1000, x, 1000))

K_PGM_mvFMR_SCAD_indepX_500 <- sapply(sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_SCAD_indepX_1000 <- sapply(sim3_mod8_PGM_mvFMR_SCAD$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_SCAD_corX_500 <- sapply(sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_SCAD_corX_1000 <- sapply(sim3_mod8_PGM_mvFMR_SCAD$optimal$corX$n1000, function(x) length(x$Bk))

## PGM_mvFMR_MCP ##
dens_PGM_mvFMR_MCP_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                        sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n500[[x]]))
pred_ll_PGM_mvFMR_MCP_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n500, 
                                                                            dens_PGM_MCP_indepX_500, x, 500))

dens_PGM_mvFMR_MCP_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                         sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n1000[[x]]))
pred_ll_PGM_mvFMR_MCP_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n1000,
                                                                             dens_PGM_MCP_indepX_1000, x, 1000))

dens_PGM_mvFMR_MCP_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                      sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n500[[x]]))
pred_ll_PGM_mvFMR_MCP_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n500, 
                                                                          dens_PGM_MCP_corX_500, x, 500))

dens_PGM_mvFMR_MCP_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                       sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n1000[[x]]))
pred_ll_PGM_mvFMR_MCP_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n1000, 
                                                                           dens_PGM_MCP_corX_1000, x, 1000))

K_PGM_mvFMR_MCP_indepX_500 <- sapply(sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_MCP_indepX_1000 <- sapply(sim3_mod8_PGM_mvFMR_MCP$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_PGM_mvFMR_MCP_corX_500 <- sapply(sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n500, function(x) length(x$Bk))
K_PGM_mvFMR_MCP_corX_1000 <- sapply(sim3_mod8_PGM_mvFMR_MCP$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_LASSO ##
dens_ADMM_mvFMR_LASSO_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                           sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_LASSO_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n500, 
                                                                               dens_ADMM_LASSO_indepX_500, x, 500))

dens_ADMM_mvFMR_LASSO_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                            sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_LASSO_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n1000,
                                                                                dens_ADMM_LASSO_indepX_1000, x, 1000))

dens_ADMM_mvFMR_LASSO_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                         sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_LASSO_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n500, 
                                                                             dens_ADMM_LASSO_corX_500, x, 500))

dens_ADMM_mvFMR_LASSO_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                          sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_LASSO_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n1000, 
                                                                              dens_ADMM_LASSO_corX_1000, x, 1000))

K_ADMM_mvFMR_LASSO_indepX_500 <- sapply(sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_LASSO_indepX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_LASSO$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_LASSO_corX_500 <- sapply(sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_LASSO_corX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_LASSO$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_SCAD ##
dens_ADMM_mvFMR_SCAD_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                          sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_SCAD_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n500, 
                                                                              dens_ADMM_SCAD_indepX_500, x, 500))

dens_ADMM_mvFMR_SCAD_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                           sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_SCAD_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n1000,
                                                                               dens_ADMM_SCAD_indepX_1000, x, 1000))

dens_ADMM_mvFMR_SCAD_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                        sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_SCAD_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n500, 
                                                                            dens_ADMM_SCAD_corX_500, x, 500))

dens_ADMM_mvFMR_SCAD_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                         sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_SCAD_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n1000, 
                                                                             dens_ADMM_SCAD_corX_1000, x, 1000))

K_ADMM_mvFMR_SCAD_indepX_500 <- sapply(sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_SCAD_indepX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_SCAD$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_SCAD_corX_500 <- sapply(sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_SCAD_corX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_SCAD$optimal$corX$n1000, function(x) length(x$Bk))

## ADMM_mvFMR_MCP ##
dens_ADMM_mvFMR_MCP_indepX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, TRUE), 
                                                                         sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n500[[x]]))
pred_ll_ADMM_mvFMR_MCP_indepX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n500, 
                                                                             dens_ADMM_MCP_indepX_500, x, 500))

dens_ADMM_mvFMR_MCP_indepX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, TRUE), 
                                                                          sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n1000[[x]]))
pred_ll_ADMM_mvFMR_MCP_indepX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n1000,
                                                                              dens_ADMM_MCP_indepX_1000, x, 1000))

dens_ADMM_mvFMR_MCP_corX_500 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 500, FALSE), 
                                                                       sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n500[[x]]))
pred_ll_ADMM_mvFMR_MCP_corX_500 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n500, 
                                                                           dens_ADMM_MCP_corX_500, x, 500))

dens_ADMM_mvFMR_MCP_corX_1000 <- lapply(1:500, function(x) pred_ll_dens(data_generate_5.3.2(seed_number_5.3.2[x], 1000, FALSE), 
                                                                        sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n1000[[x]]))
pred_ll_ADMM_mvFMR_MCP_corX_1000 <- sapply(1:500, function(x) predictive_ll(sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n1000, 
                                                                            dens_ADMM_MCP_corX_1000, x, 1000))

K_ADMM_mvFMR_MCP_indepX_500 <- sapply(sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_MCP_indepX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_MCP$optimal$indepX$n1000, function(x) length(x$Bk)) 
K_ADMM_mvFMR_MCP_corX_500 <- sapply(sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n500, function(x) length(x$Bk))
K_ADMM_mvFMR_MCP_corX_1000 <- sapply(sim3_mod8_ADMM_mvFMR_MCP$optimal$corX$n1000, function(x) length(x$Bk))

## Table ##
model8_table_mean <- list(indepX=data.frame(n500=vector(length=7), n1000=vector(length=7)), 
                          corX=data.frame(n500=vector(length=7), n1000=vector(length=7)))

rownames(model8_table_mean$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                        "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model8_table_mean$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                      "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model8_table_mean$indepX$n500 <- c(mean(pred_ll_UNPEN_indepX_500), mean(pred_ll_PGM_LASSO_indepX_500), 
                                   mean(pred_ll_ADMM_LASSO_indepX_500), mean(pred_ll_PGM_SCAD_indepX_500),
                                   mean(pred_ll_ADMM_SCAD_indepX_500), mean(pred_ll_PGM_MCP_indepX_500),
                                   mean(pred_ll_ADMM_MCP_indepX_500))
model8_table_mean$indepX$n1000 <- c(mean(pred_ll_UNPEN_indepX_1000), mean(pred_ll_PGM_LASSO_indepX_1000),
                                    mean(pred_ll_ADMM_LASSO_indepX_1000), mean(pred_ll_PGM_SCAD_indepX_1000),
                                    mean(pred_ll_ADMM_SCAD_indepX_1000), mean(pred_ll_PGM_MCP_indepX_1000), 
                                    mean(pred_ll_ADMM_MCP_indepX_1000))
model8_table_mean$corX$n500 <- c(mean(pred_ll_UNPEN_corX_500), mean(pred_ll_PGM_LASSO_corX_500), 
                                 mean(pred_ll_ADMM_LASSO_corX_500), mean(pred_ll_PGM_SCAD_corX_500),
                                 mean(pred_ll_ADMM_SCAD_corX_500), mean(pred_ll_PGM_MCP_corX_500),
                                 mean(pred_ll_ADMM_MCP_corX_500))
model8_table_mean$corX$n1000 <- c(mean(pred_ll_UNPEN_corX_1000), mean(pred_ll_PGM_LASSO_corX_1000), 
                                  mean(pred_ll_ADMM_LASSO_corX_1000), mean(pred_ll_PGM_SCAD_corX_1000),
                                  mean(pred_ll_ADMM_SCAD_corX_1000), mean(pred_ll_PGM_MCP_corX_1000),
                                  mean(pred_ll_ADMM_MCP_corX_1000))

model8_table_sd <- list(indepX=data.frame(n500=vector(length=7), n1000=vector(length=7)), 
                        corX=data.frame(n500=vector(length=7), n1000=vector(length=7)))
rownames(model8_table_sd$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                      "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model8_table_sd$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                    "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model8_table_sd$indepX$n500 <- c(sd(pred_ll_UNPEN_indepX_500), sd(pred_ll_PGM_LASSO_indepX_500), 
                                 sd(pred_ll_ADMM_LASSO_indepX_500), sd(pred_ll_PGM_SCAD_indepX_500),
                                 sd(pred_ll_ADMM_SCAD_indepX_500), sd(pred_ll_PGM_MCP_indepX_500),
                                 sd(pred_ll_ADMM_MCP_indepX_500))
model8_table_sd$indepX$n1000 <- c(sd(pred_ll_UNPEN_indepX_1000), sd(pred_ll_PGM_LASSO_indepX_1000), 
                                  sd(pred_ll_ADMM_LASSO_indepX_1000), sd(pred_ll_PGM_SCAD_indepX_1000),
                                  sd(pred_ll_ADMM_SCAD_indepX_1000), sd(pred_ll_PGM_MCP_indepX_1000), 
                                  sd(pred_ll_ADMM_MCP_indepX_1000))
model8_table_sd$corX$n500 <- c(sd(pred_ll_UNPEN_corX_500), sd(pred_ll_PGM_LASSO_corX_500), 
                               sd(pred_ll_ADMM_LASSO_corX_500), sd(pred_ll_PGM_SCAD_corX_500),
                               sd(pred_ll_ADMM_SCAD_corX_500), sd(pred_ll_PGM_MCP_corX_500), 
                               sd(pred_ll_ADMM_MCP_corX_500))
model8_table_sd$corX$n1000 <- c(sd(pred_ll_UNPEN_corX_1000), sd(pred_ll_PGM_LASSO_corX_1000), 
                                sd(pred_ll_ADMM_LASSO_corX_1000), sd(pred_ll_PGM_SCAD_corX_1000),
                                sd(pred_ll_ADMM_SCAD_corX_1000), sd(pred_ll_PGM_MCP_corX_1000),
                                sd(pred_ll_ADMM_MCP_corX_1000))

model8_table <- list(indepX=data.frame(n500=rep(NA,7), n1000=rep(NA,7)), 
                     corX=data.frame(n500=rep(NA,7), n1000=rep(NA,7)))
rownames(model8_table$indepX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM",
                                   "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
rownames(model8_table$corX) <- c("mvFMR", "mvFMR-L.PGM", "mvFMR-L.ADMM", "mvFMR-S.PGM", 
                                 "mvFMR-S.ADMM", "mvFMR-M.PGM", "mvFMR-M.ADMM")
model8_table$indepX$n500 <- sprintf("%.4f(%.4f)", model8_table_mean$indepX$n500, model8_table_sd$indepX$n500)
model8_table$indepX$n1000 <- sprintf("%.4f(%.4f)", model8_table_mean$indepX$n1000, model8_table_sd$indepX$n1000)
model8_table$corX$n500 <- sprintf("%.4f(%.4f)", model8_table_mean$corX$n500, model8_table_sd$corX$n500)
model8_table$corX$n1000 <- sprintf("%.4f(%.4f)", model8_table_mean$corX$n1000, model8_table_sd$corX$n1000)
# mod8_table <- model8_table
# save(mod8_table, file="./Results/Section5.3/Model8/mod8_table.rda")

mod8_K_table <- list(indepX=list(n500=data.frame(mvFMR=K_unpen_indepX_500, 
                                                 PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_indepX_500,
                                                 ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_indepX_500,
                                                 PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_indepX_500,
                                                 ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_indepX_500,
                                                 PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_indepX_500,
                                                 ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_indepX_500),
                                 n1000=data.frame(mvFMR=K_unpen_indepX_1000, 
                                                  PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_indepX_1000,
                                                  ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_indepX_1000,
                                                  PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_indepX_1000,
                                                  ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_indepX_1000,
                                                  PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_indepX_1000,
                                                  ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_indepX_1000)),
                     corX=list(n500=data.frame(mvFMR=K_unpen_corX_500, 
                                               PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_corX_500,
                                               ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_corX_500,
                                               PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_corX_500,
                                               ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_corX_500,
                                               PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_corX_500,
                                               ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_corX_500),
                               n1000=data.frame(mvFMR=K_unpen_corX_1000, 
                                                PGM_mvFMR_LASSO=K_PGM_mvFMR_LASSO_corX_1000,
                                                ADMM_mvFMR_LASSO=K_ADMM_mvFMR_LASSO_corX_1000,
                                                PGM_mvFMR_SCAD=K_PGM_mvFMR_SCAD_corX_1000,
                                                ADMM_mvFMR_SCAD=K_ADMM_mvFMR_SCAD_corX_1000,
                                                PGM_mvFMR_MCP=K_PGM_mvFMR_MCP_corX_1000,
                                                ADMM_mvFMR_MCP=K_ADMM_mvFMR_MCP_corX_1000)))
# save(mod8_K_table, file="./Results/Section5.3/Model8/mod8_K_table.rda")