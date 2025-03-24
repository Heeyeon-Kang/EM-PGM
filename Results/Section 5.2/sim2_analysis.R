###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####              R-code for reproducing Table 4 in Section 5.2.               ####
###################################################################################

# The following R-code contains functions for calculating TPR, FPR, MSE.

# Each ".rda" file contains the optimal estimators('pi', 'Bk', 'sigma') and 
# the estimators of mixing proportion 'w'.

### Data ###
source("./Data/simulation_data.R")

### Functions ###
source("./Scripts/Functions/simulation_summary_functions.R")

### Model 4 ###
## results ##
load("./Results/Section 5.2/Model 4/sim2_mod4_flexmix_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_oracle_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_mvFMR_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_EM_PGM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_EM_PGM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_EM_PGM_mvFMR_MCP_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_ADMM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_ADMM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 4/sim2_mod4_ADMM_mvFMR_MCP_result.rda")
true <- data_generate_5.2.1(1, 500)$true
K <- 2

## flexmix ##
TPR_flexmix_500 <- sapply(sim2_mod4_flexmix_result$n500$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_500 <- sapply(sim2_mod4_flexmix_result$n500$optimal, function(x) FPR_flexmix(true, x))
TPR_flexmix_1000 <- sapply(sim2_mod4_flexmix_result$n1000$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_1000 <- sapply(sim2_mod4_flexmix_result$n1000$optimal, function(x) FPR_flexmix(true, x))
MSE_pi_flexmix_500 <- sapply(sim2_mod4_flexmix_result$n500$MSE, function(x) x[1])
MSE_Bk_flexmix_500 <- sapply(sim2_mod4_flexmix_result$n500$MSE, function(x) x[2])
MSE_sigma_flexmix_500 <- sapply(sim2_mod4_flexmix_result$n500$MSE, function(x) x[3])
MSE_flexmix_500 <- data.frame(pi=MSE_pi_flexmix_500, Bk=MSE_Bk_flexmix_500, sigma=MSE_sigma_flexmix_500)
MSE_pi_flexmix_1000 <- sapply(sim2_mod4_flexmix_result$n1000$MSE, function(x) x[1])
MSE_Bk_flexmix_1000 <- sapply(sim2_mod4_flexmix_result$n1000$MSE, function(x) x[2])
MSE_sigma_flexmix_1000 <- sapply(sim2_mod4_flexmix_result$n1000$MSE, function(x) x[3])
MSE_flexmix_1000 <- data.frame(pi=MSE_pi_flexmix_1000, Bk=MSE_Bk_flexmix_1000, sigma=MSE_sigma_flexmix_1000)

## oracle ##
TPR_oracle_500 <- sapply(sim2_mod4_oracle_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_500 <- sapply(sim2_mod4_oracle_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_oracle_1000 <- sapply(sim2_mod4_oracle_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_1000 <- sapply(sim2_mod4_oracle_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_oracle_500 <- as.data.frame(t(sapply(sim2_mod4_oracle_result$n500$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_500) <- c("pi", "Bk", "sigma")
MSE_oracle_1000 <- as.data.frame(t(sapply(sim2_mod4_oracle_result$n1000$optimal, 
                                          function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_1000) <- c("pi", "Bk", "sigma")

## mvFMR ##
TPR_mvFMR_500 <- sapply(sim2_mod4_mvFMR_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_500 <- sapply(sim2_mod4_mvFMR_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_mvFMR_1000 <- sapply(sim2_mod4_mvFMR_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_1000 <- sapply(sim2_mod4_mvFMR_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_mvFMR_500 <- as.data.frame(t(sapply(sim2_mod4_mvFMR_result$n500$optimal, 
                                        function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_500) <- c("pi", "Bk", "sigma")
MSE_mvFMR_1000 <- as.data.frame(t(sapply(sim2_mod4_mvFMR_result$n1000$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_LASSO ##
TPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                                      function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_SCAD ##
TPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n500$optimal,
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_MCP ##
TPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod4_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_LASSO ##
TPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_ADMM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_SCAD ##
TPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) FPR(true, x))
TPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) FPR(true, x))
MSE_ADMM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_MCP ##
TPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) FPR(true, x))
TPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) FPR(true, x))
MSE_ADMM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n500$optimal, 
                                                 function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod4_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## Table ##
model4_table <- list(n500=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                               sd=as.data.frame(matrix(nrow=9, ncol=5))),
                     n1000=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                                sd=as.data.frame(matrix(nrow=9, ncol=5))))
colnames(model4_table$n500$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model4_table$n500$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                      "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model4_table$n500$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model4_table$n500$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                    "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model4_table$n1000$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model4_table$n1000$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                       "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model4_table$n1000$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model4_table$n1000$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                     "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")

model4_table$n500$mean$TPR <- c(mean(TPR_flexmix_500), mean(TPR_oracle_500), mean(TPR_mvFMR_500),
                                mean(TPR_EM_PGM_mvFMR_LASSO_500), mean(TPR_ADMM_mvFMR_LASSO_500),
                                mean(TPR_EM_PGM_mvFMR_SCAD_500), mean(TPR_ADMM_mvFMR_SCAD_500),
                                mean(TPR_EM_PGM_mvFMR_MCP_500), mean(TPR_ADMM_mvFMR_MCP_500))
model4_table$n500$mean$FPR <- c(mean(FPR_flexmix_500), mean(FPR_oracle_500), mean(FPR_mvFMR_500),
                                mean(FPR_EM_PGM_mvFMR_LASSO_500), mean(FPR_ADMM_mvFMR_LASSO_500),
                                mean(FPR_EM_PGM_mvFMR_SCAD_500), mean(FPR_ADMM_mvFMR_SCAD_500),
                                mean(FPR_EM_PGM_mvFMR_MCP_500), mean(FPR_ADMM_mvFMR_MCP_500))
model4_table$n500$mean$MSE_pi_k <- c(mean(MSE_flexmix_500$pi), mean(MSE_oracle_500$pi), 
                                     mean(MSE_mvFMR_500$pi), mean(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                     mean(MSE_ADMM_mvFMR_LASSO_500$pi), mean(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                     mean(MSE_ADMM_mvFMR_SCAD_500$pi), mean(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                     mean(MSE_ADMM_mvFMR_MCP_500$pi))
model4_table$n500$mean$MSE_B_k <- c(mean(MSE_flexmix_500$Bk), mean(MSE_oracle_500$Bk), 
                                    mean(MSE_mvFMR_500$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                    mean(MSE_ADMM_mvFMR_LASSO_500$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                    mean(MSE_ADMM_mvFMR_SCAD_500$Bk), mean(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                    mean(MSE_ADMM_mvFMR_MCP_500$Bk))
model4_table$n500$mean$MSE_sigma_k <- c(mean(MSE_flexmix_500$sigma), mean(MSE_oracle_500$sigma), 
                                        mean(MSE_mvFMR_500$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                        mean(MSE_ADMM_mvFMR_LASSO_500$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                        mean(MSE_ADMM_mvFMR_SCAD_500$sigma), mean(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                        mean(MSE_ADMM_mvFMR_MCP_500$sigma))
model4_table$n500$sd$TPR <- c(sd(TPR_flexmix_500), sd(TPR_oracle_500), sd(TPR_mvFMR_500),
                              sd(TPR_EM_PGM_mvFMR_LASSO_500), sd(TPR_ADMM_mvFMR_LASSO_500),
                              sd(TPR_EM_PGM_mvFMR_SCAD_500), sd(TPR_ADMM_mvFMR_SCAD_500),
                              sd(TPR_EM_PGM_mvFMR_MCP_500), sd(TPR_ADMM_mvFMR_MCP_500))
model4_table$n500$sd$FPR <- c(sd(FPR_flexmix_500), sd(FPR_oracle_500), sd(FPR_mvFMR_500),
                              sd(FPR_EM_PGM_mvFMR_LASSO_500), sd(FPR_ADMM_mvFMR_LASSO_500),
                              sd(FPR_EM_PGM_mvFMR_SCAD_500), sd(FPR_ADMM_mvFMR_SCAD_500),
                              sd(FPR_EM_PGM_mvFMR_MCP_500), sd(FPR_ADMM_mvFMR_MCP_500))
model4_table$n500$sd$MSE_pi_k <- c(sd(MSE_flexmix_500$pi), sd(MSE_oracle_500$pi), 
                                   sd(MSE_mvFMR_500$pi), sd(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                   sd(MSE_ADMM_mvFMR_LASSO_500$pi), sd(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                   sd(MSE_ADMM_mvFMR_SCAD_500$pi), sd(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                   sd(MSE_ADMM_mvFMR_MCP_500$pi))
model4_table$n500$sd$MSE_B_k <- c(sd(MSE_flexmix_500$Bk), sd(MSE_oracle_500$Bk), 
                                  sd(MSE_mvFMR_500$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                  sd(MSE_ADMM_mvFMR_LASSO_500$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                  sd(MSE_ADMM_mvFMR_SCAD_500$Bk), sd(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                  sd(MSE_ADMM_mvFMR_MCP_500$Bk))
model4_table$n500$sd$MSE_sigma_k <- c(sd(MSE_flexmix_500$sigma), sd(MSE_oracle_500$sigma), 
                                      sd(MSE_mvFMR_500$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                      sd(MSE_ADMM_mvFMR_LASSO_500$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                      sd(MSE_ADMM_mvFMR_SCAD_500$sigma), sd(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                      sd(MSE_ADMM_mvFMR_MCP_500$sigma))
model4_table$n1000$mean$TPR <- c(mean(TPR_flexmix_1000), mean(TPR_oracle_1000), mean(TPR_mvFMR_1000),
                                 mean(TPR_EM_PGM_mvFMR_LASSO_1000), mean(TPR_ADMM_mvFMR_LASSO_1000),
                                 mean(TPR_EM_PGM_mvFMR_SCAD_1000), mean(TPR_ADMM_mvFMR_SCAD_1000),
                                 mean(TPR_EM_PGM_mvFMR_MCP_1000), mean(TPR_ADMM_mvFMR_MCP_1000))
model4_table$n1000$mean$FPR <- c(mean(FPR_flexmix_1000), mean(FPR_oracle_1000), mean(FPR_mvFMR_1000),
                                 mean(FPR_EM_PGM_mvFMR_LASSO_1000), mean(FPR_ADMM_mvFMR_LASSO_1000),
                                 mean(FPR_EM_PGM_mvFMR_SCAD_1000), mean(FPR_ADMM_mvFMR_SCAD_1000),
                                 mean(FPR_EM_PGM_mvFMR_MCP_1000), mean(FPR_ADMM_mvFMR_MCP_1000))
model4_table$n1000$mean$MSE_pi_k <- c(mean(MSE_flexmix_1000$pi), mean(MSE_oracle_1000$pi), 
                                      mean(MSE_mvFMR_1000$pi), mean(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                      mean(MSE_ADMM_mvFMR_LASSO_1000$pi), mean(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                      mean(MSE_ADMM_mvFMR_SCAD_1000$pi), mean(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                      mean(MSE_ADMM_mvFMR_MCP_1000$pi))
model4_table$n1000$mean$MSE_B_k <- c(mean(MSE_flexmix_1000$Bk), mean(MSE_oracle_1000$Bk), 
                                     mean(MSE_mvFMR_1000$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_LASSO_1000$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_SCAD_1000$Bk), mean(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_MCP_1000$Bk))
model4_table$n1000$mean$MSE_sigma_k <- c(mean(MSE_flexmix_1000$sigma), mean(MSE_oracle_1000$sigma), 
                                         mean(MSE_mvFMR_1000$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_LASSO_1000$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_SCAD_1000$sigma), mean(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_MCP_1000$sigma))
model4_table$n1000$sd$TPR <- c(sd(TPR_flexmix_1000), sd(TPR_oracle_1000), sd(TPR_mvFMR_1000),
                               sd(TPR_EM_PGM_mvFMR_LASSO_1000), sd(TPR_ADMM_mvFMR_LASSO_1000),
                               sd(TPR_EM_PGM_mvFMR_SCAD_1000), sd(TPR_ADMM_mvFMR_SCAD_1000),
                               sd(TPR_EM_PGM_mvFMR_MCP_1000), sd(TPR_ADMM_mvFMR_MCP_1000))
model4_table$n1000$sd$FPR <- c(sd(FPR_flexmix_1000), sd(FPR_oracle_1000), sd(FPR_mvFMR_1000),
                               sd(FPR_EM_PGM_mvFMR_LASSO_1000), sd(FPR_ADMM_mvFMR_LASSO_1000),
                               sd(FPR_EM_PGM_mvFMR_SCAD_1000), sd(FPR_ADMM_mvFMR_SCAD_1000),
                               sd(FPR_EM_PGM_mvFMR_MCP_1000), sd(FPR_ADMM_mvFMR_MCP_1000))
model4_table$n1000$sd$MSE_pi_k <- c(sd(MSE_flexmix_1000$pi), sd(MSE_oracle_1000$pi), 
                                    sd(MSE_mvFMR_1000$pi), sd(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                    sd(MSE_ADMM_mvFMR_LASSO_1000$pi), sd(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                    sd(MSE_ADMM_mvFMR_SCAD_1000$pi), sd(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                    sd(MSE_ADMM_mvFMR_MCP_1000$pi))
model4_table$n1000$sd$MSE_B_k <- c(sd(MSE_flexmix_1000$Bk), sd(MSE_oracle_1000$Bk), 
                                   sd(MSE_mvFMR_1000$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_LASSO_1000$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_SCAD_1000$Bk), sd(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_MCP_1000$Bk))
model4_table$n1000$sd$MSE_sigma_k <- c(sd(MSE_flexmix_1000$sigma), sd(MSE_oracle_1000$sigma), 
                                       sd(MSE_mvFMR_1000$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_LASSO_1000$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_SCAD_1000$sigma), sd(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_MCP_1000$sigma))


### Model 5 ###
## results ##
load("./Results/Section 5.2/Model 5/sim2_mod5_flexmix_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_oracle_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_mvFMR_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_EM_PGM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_EM_PGM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_EM_PGM_mvFMR_MCP_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_ADMM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_ADMM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 5/sim2_mod5_ADMM_mvFMR_MCP_result.rda")
true <- data_generate_5.2.2(1, 500)$true
K <- 2

## flexmix ##
TPR_flexmix_500 <- sapply(sim2_mod5_flexmix_result$n500$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_500 <- sapply(sim2_mod5_flexmix_result$n500$optimal, function(x) FPR_flexmix(true, x))
TPR_flexmix_1000 <- sapply(sim2_mod5_flexmix_result$n1000$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_1000 <- sapply(sim2_mod5_flexmix_result$n1000$optimal, function(x) FPR_flexmix(true, x))
MSE_pi_flexmix_500 <- sapply(sim2_mod5_flexmix_result$n500$MSE, function(x) x[1])
MSE_Bk_flexmix_500 <- sapply(sim2_mod5_flexmix_result$n500$MSE, function(x) x[2])
MSE_sigma_flexmix_500 <- sapply(sim2_mod5_flexmix_result$n500$MSE, function(x) x[3])
MSE_flexmix_500 <- data.frame(pi=MSE_pi_flexmix_500, Bk=MSE_Bk_flexmix_500, sigma=MSE_sigma_flexmix_500)
MSE_pi_flexmix_1000 <- sapply(sim2_mod5_flexmix_result$n1000$MSE, function(x) x[1])
MSE_Bk_flexmix_1000 <- sapply(sim2_mod5_flexmix_result$n1000$MSE, function(x) x[2])
MSE_sigma_flexmix_1000 <- sapply(sim2_mod5_flexmix_result$n1000$MSE, function(x) x[3])
MSE_flexmix_1000 <- data.frame(pi=MSE_pi_flexmix_1000, Bk=MSE_Bk_flexmix_1000, sigma=MSE_sigma_flexmix_1000)

## oracle ##
TPR_oracle_500 <- sapply(sim2_mod5_oracle_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_500 <- sapply(sim2_mod5_oracle_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_oracle_1000 <- sapply(sim2_mod5_oracle_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_1000 <- sapply(sim2_mod5_oracle_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_oracle_500 <- as.data.frame(t(sapply(sim2_mod5_oracle_result$n500$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_500) <- c("pi", "Bk", "sigma")
MSE_oracle_1000 <- as.data.frame(t(sapply(sim2_mod5_oracle_result$n1000$optimal, 
                                          function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_1000) <- c("pi", "Bk", "sigma")

## mvFMR ##
TPR_mvFMR_500 <- sapply(sim2_mod5_mvFMR_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_500 <- sapply(sim2_mod5_mvFMR_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_mvFMR_1000 <- sapply(sim2_mod5_mvFMR_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_1000 <- sapply(sim2_mod5_mvFMR_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_mvFMR_500 <- as.data.frame(t(sapply(sim2_mod5_mvFMR_result$n500$optimal, 
                                        function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_500) <- c("pi", "Bk", "sigma")
MSE_mvFMR_1000 <- as.data.frame(t(sapply(sim2_mod5_mvFMR_result$n1000$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_LASSO ##
TPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                                      function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_SCAD ##
TPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n500$optimal,
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_MCP ##
TPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod5_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_LASSO ##
TPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_ADMM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_SCAD ##
TPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) FPR(true, x))
TPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) FPR(true, x))
MSE_ADMM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_MCP ##
TPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) FPR(true, x))
TPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) FPR(true, x))
MSE_ADMM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n500$optimal, 
                                                 function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod5_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## Table ##
model5_table <- list(n500=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                               sd=as.data.frame(matrix(nrow=9, ncol=5))),
                     n1000=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                                sd=as.data.frame(matrix(nrow=9, ncol=5))))
colnames(model5_table$n500$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model5_table$n500$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                      "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model5_table$n500$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model5_table$n500$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                    "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model5_table$n1000$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model5_table$n1000$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                       "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model5_table$n1000$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model5_table$n1000$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                     "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")

model5_table$n500$mean$TPR <- c(mean(TPR_flexmix_500), mean(TPR_oracle_500), mean(TPR_mvFMR_500),
                                mean(TPR_EM_PGM_mvFMR_LASSO_500), mean(TPR_ADMM_mvFMR_LASSO_500),
                                mean(TPR_EM_PGM_mvFMR_SCAD_500), mean(TPR_ADMM_mvFMR_SCAD_500),
                                mean(TPR_EM_PGM_mvFMR_MCP_500), mean(TPR_ADMM_mvFMR_MCP_500))
model5_table$n500$mean$FPR <- c(mean(FPR_flexmix_500), mean(FPR_oracle_500), mean(FPR_mvFMR_500),
                                mean(FPR_EM_PGM_mvFMR_LASSO_500), mean(FPR_ADMM_mvFMR_LASSO_500),
                                mean(FPR_EM_PGM_mvFMR_SCAD_500), mean(FPR_ADMM_mvFMR_SCAD_500),
                                mean(FPR_EM_PGM_mvFMR_MCP_500), mean(FPR_ADMM_mvFMR_MCP_500))
model5_table$n500$mean$MSE_pi_k <- c(mean(MSE_flexmix_500$pi), mean(MSE_oracle_500$pi), 
                                     mean(MSE_mvFMR_500$pi), mean(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                     mean(MSE_ADMM_mvFMR_LASSO_500$pi), mean(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                     mean(MSE_ADMM_mvFMR_SCAD_500$pi), mean(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                     mean(MSE_ADMM_mvFMR_MCP_500$pi))
model5_table$n500$mean$MSE_B_k <- c(mean(MSE_flexmix_500$Bk), mean(MSE_oracle_500$Bk), 
                                    mean(MSE_mvFMR_500$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                    mean(MSE_ADMM_mvFMR_LASSO_500$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                    mean(MSE_ADMM_mvFMR_SCAD_500$Bk), mean(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                    mean(MSE_ADMM_mvFMR_MCP_500$Bk))
model5_table$n500$mean$MSE_sigma_k <- c(mean(MSE_flexmix_500$sigma), mean(MSE_oracle_500$sigma), 
                                        mean(MSE_mvFMR_500$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                        mean(MSE_ADMM_mvFMR_LASSO_500$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                        mean(MSE_ADMM_mvFMR_SCAD_500$sigma), mean(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                        mean(MSE_ADMM_mvFMR_MCP_500$sigma))
model5_table$n500$sd$TPR <- c(sd(TPR_flexmix_500), sd(TPR_oracle_500), sd(TPR_mvFMR_500),
                              sd(TPR_EM_PGM_mvFMR_LASSO_500), sd(TPR_ADMM_mvFMR_LASSO_500),
                              sd(TPR_EM_PGM_mvFMR_SCAD_500), sd(TPR_ADMM_mvFMR_SCAD_500),
                              sd(TPR_EM_PGM_mvFMR_MCP_500), sd(TPR_ADMM_mvFMR_MCP_500))
model5_table$n500$sd$FPR <- c(sd(FPR_flexmix_500), sd(FPR_oracle_500), sd(FPR_mvFMR_500),
                              sd(FPR_EM_PGM_mvFMR_LASSO_500), sd(FPR_ADMM_mvFMR_LASSO_500),
                              sd(FPR_EM_PGM_mvFMR_SCAD_500), sd(FPR_ADMM_mvFMR_SCAD_500),
                              sd(FPR_EM_PGM_mvFMR_MCP_500), sd(FPR_ADMM_mvFMR_MCP_500))
model5_table$n500$sd$MSE_pi_k <- c(sd(MSE_flexmix_500$pi), sd(MSE_oracle_500$pi), 
                                   sd(MSE_mvFMR_500$pi), sd(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                   sd(MSE_ADMM_mvFMR_LASSO_500$pi), sd(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                   sd(MSE_ADMM_mvFMR_SCAD_500$pi), sd(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                   sd(MSE_ADMM_mvFMR_MCP_500$pi))
model5_table$n500$sd$MSE_B_k <- c(sd(MSE_flexmix_500$Bk), sd(MSE_oracle_500$Bk), 
                                  sd(MSE_mvFMR_500$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                  sd(MSE_ADMM_mvFMR_LASSO_500$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                  sd(MSE_ADMM_mvFMR_SCAD_500$Bk), sd(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                  sd(MSE_ADMM_mvFMR_MCP_500$Bk))
model5_table$n500$sd$MSE_sigma_k <- c(sd(MSE_flexmix_500$sigma), sd(MSE_oracle_500$sigma), 
                                      sd(MSE_mvFMR_500$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                      sd(MSE_ADMM_mvFMR_LASSO_500$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                      sd(MSE_ADMM_mvFMR_SCAD_500$sigma), sd(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                      sd(MSE_ADMM_mvFMR_MCP_500$sigma))
model5_table$n1000$mean$TPR <- c(mean(TPR_flexmix_1000), mean(TPR_oracle_1000), mean(TPR_mvFMR_1000),
                                 mean(TPR_EM_PGM_mvFMR_LASSO_1000), mean(TPR_ADMM_mvFMR_LASSO_1000),
                                 mean(TPR_EM_PGM_mvFMR_SCAD_1000), mean(TPR_ADMM_mvFMR_SCAD_1000),
                                 mean(TPR_EM_PGM_mvFMR_MCP_1000), mean(TPR_ADMM_mvFMR_MCP_1000))
model5_table$n1000$mean$FPR <- c(mean(FPR_flexmix_1000), mean(FPR_oracle_1000), mean(FPR_mvFMR_1000),
                                 mean(FPR_EM_PGM_mvFMR_LASSO_1000), mean(FPR_ADMM_mvFMR_LASSO_1000),
                                 mean(FPR_EM_PGM_mvFMR_SCAD_1000), mean(FPR_ADMM_mvFMR_SCAD_1000),
                                 mean(FPR_EM_PGM_mvFMR_MCP_1000), mean(FPR_ADMM_mvFMR_MCP_1000))
model5_table$n1000$mean$MSE_pi_k <- c(mean(MSE_flexmix_1000$pi), mean(MSE_oracle_1000$pi), 
                                      mean(MSE_mvFMR_1000$pi), mean(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                      mean(MSE_ADMM_mvFMR_LASSO_1000$pi), mean(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                      mean(MSE_ADMM_mvFMR_SCAD_1000$pi), mean(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                      mean(MSE_ADMM_mvFMR_MCP_1000$pi))
model5_table$n1000$mean$MSE_B_k <- c(mean(MSE_flexmix_1000$Bk), mean(MSE_oracle_1000$Bk), 
                                     mean(MSE_mvFMR_1000$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_LASSO_1000$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_SCAD_1000$Bk), mean(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_MCP_1000$Bk))
model5_table$n1000$mean$MSE_sigma_k <- c(mean(MSE_flexmix_1000$sigma), mean(MSE_oracle_1000$sigma), 
                                         mean(MSE_mvFMR_1000$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_LASSO_1000$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_SCAD_1000$sigma), mean(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_MCP_1000$sigma))
model5_table$n1000$sd$TPR <- c(sd(TPR_flexmix_1000), sd(TPR_oracle_1000), sd(TPR_mvFMR_1000),
                               sd(TPR_EM_PGM_mvFMR_LASSO_1000), sd(TPR_ADMM_mvFMR_LASSO_1000),
                               sd(TPR_EM_PGM_mvFMR_SCAD_1000), sd(TPR_ADMM_mvFMR_SCAD_1000),
                               sd(TPR_EM_PGM_mvFMR_MCP_1000), sd(TPR_ADMM_mvFMR_MCP_1000))
model5_table$n1000$sd$FPR <- c(sd(FPR_flexmix_1000), sd(FPR_oracle_1000), sd(FPR_mvFMR_1000),
                               sd(FPR_EM_PGM_mvFMR_LASSO_1000), sd(FPR_ADMM_mvFMR_LASSO_1000),
                               sd(FPR_EM_PGM_mvFMR_SCAD_1000), sd(FPR_ADMM_mvFMR_SCAD_1000),
                               sd(FPR_EM_PGM_mvFMR_MCP_1000), sd(FPR_ADMM_mvFMR_MCP_1000))
model5_table$n1000$sd$MSE_pi_k <- c(sd(MSE_flexmix_1000$pi), sd(MSE_oracle_1000$pi), 
                                    sd(MSE_mvFMR_1000$pi), sd(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                    sd(MSE_ADMM_mvFMR_LASSO_1000$pi), sd(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                    sd(MSE_ADMM_mvFMR_SCAD_1000$pi), sd(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                    sd(MSE_ADMM_mvFMR_MCP_1000$pi))
model5_table$n1000$sd$MSE_B_k <- c(sd(MSE_flexmix_1000$Bk), sd(MSE_oracle_1000$Bk), 
                                   sd(MSE_mvFMR_1000$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_LASSO_1000$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_SCAD_1000$Bk), sd(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_MCP_1000$Bk))
model5_table$n1000$sd$MSE_sigma_k <- c(sd(MSE_flexmix_1000$sigma), sd(MSE_oracle_1000$sigma), 
                                       sd(MSE_mvFMR_1000$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_LASSO_1000$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_SCAD_1000$sigma), sd(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_MCP_1000$sigma))


### Model 6 ###
## results ##
load("./Results/Section 5.2/Model 6/sim2_mod6_flexmix_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_oracle_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_mvFMR_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_EM_PGM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_EM_PGM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_EM_PGM_mvFMR_MCP_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_ADMM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_ADMM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.2/Model 6/sim2_mod6_ADMM_mvFMR_MCP_result.rda")
true <- data_generate_5.2.3(1, 500)$true
K <- 2

## flexmix ##
TPR_flexmix_500 <- sapply(sim2_mod6_flexmix_result$n500$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_500 <- sapply(sim2_mod6_flexmix_result$n500$optimal, function(x) FPR_flexmix(true, x))
TPR_flexmix_1000 <- sapply(sim2_mod6_flexmix_result$n1000$optimal, function(x) TPR_flexmix(true, x))
FPR_flexmix_1000 <- sapply(sim2_mod6_flexmix_result$n1000$optimal, function(x) FPR_flexmix(true, x))
MSE_pi_flexmix_500 <- sapply(sim2_mod6_flexmix_result$n500$MSE, function(x) x[1])
MSE_Bk_flexmix_500 <- sapply(sim2_mod6_flexmix_result$n500$MSE, function(x) x[2])
MSE_sigma_flexmix_500 <- sapply(sim2_mod6_flexmix_result$n500$MSE, function(x) x[3])
MSE_flexmix_500 <- data.frame(pi=MSE_pi_flexmix_500, Bk=MSE_Bk_flexmix_500, sigma=MSE_sigma_flexmix_500)
MSE_pi_flexmix_1000 <- sapply(sim2_mod6_flexmix_result$n1000$MSE, function(x) x[1])
MSE_Bk_flexmix_1000 <- sapply(sim2_mod6_flexmix_result$n1000$MSE, function(x) x[2])
MSE_sigma_flexmix_1000 <- sapply(sim2_mod6_flexmix_result$n1000$MSE, function(x) x[3])
MSE_flexmix_1000 <- data.frame(pi=MSE_pi_flexmix_1000, Bk=MSE_Bk_flexmix_1000, sigma=MSE_sigma_flexmix_1000)

## oracle ##
TPR_oracle_500 <- sapply(sim2_mod6_oracle_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_500 <- sapply(sim2_mod6_oracle_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_oracle_1000 <- sapply(sim2_mod6_oracle_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_oracle_1000 <- sapply(sim2_mod6_oracle_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_oracle_500 <- as.data.frame(t(sapply(sim2_mod6_oracle_result$n500$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_500) <- c("pi", "Bk", "sigma")
MSE_oracle_1000 <- as.data.frame(t(sapply(sim2_mod6_oracle_result$n1000$optimal, 
                                          function(x) MSE_no_penalty(true, x))))
names(MSE_oracle_1000) <- c("pi", "Bk", "sigma")

## mvFMR ##
TPR_mvFMR_500 <- sapply(sim2_mod6_mvFMR_result$n500$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_500 <- sapply(sim2_mod6_mvFMR_result$n500$optimal, function(x) FPR_no_penalty(true, x))
TPR_mvFMR_1000 <- sapply(sim2_mod6_mvFMR_result$n1000$optimal, function(x) TPR_no_penalty(true, x))
FPR_mvFMR_1000 <- sapply(sim2_mod6_mvFMR_result$n1000$optimal, function(x) FPR_no_penalty(true, x))
MSE_mvFMR_500 <- as.data.frame(t(sapply(sim2_mod6_mvFMR_result$n500$optimal, 
                                        function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_500) <- c("pi", "Bk", "sigma")
MSE_mvFMR_1000 <- as.data.frame(t(sapply(sim2_mod6_mvFMR_result$n1000$optimal, 
                                         function(x) MSE_no_penalty(true, x))))
names(MSE_mvFMR_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_LASSO ##
TPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                     function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_LASSO_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                      function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n500$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_LASSO_result$n1000$optimal, 
                                                      function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_SCAD ##
TPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n500$optimal, 
                                    function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_SCAD_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                     function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n500$optimal,
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_SCAD_result$n1000$optimal, 
                                                     function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## EM_PGM_mvFMR_MCP ##
TPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_500 <- sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_EM_PGM_mvFMR_MCP_1000 <- sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_EM_PGM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_EM_PGM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod6_EM_PGM_mvFMR_MCP_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_EM_PGM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_LASSO ##
TPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_500 <- sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                   function(x) FPR(true, x))
TPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) TPR(true, x))
FPR_ADMM_mvFMR_LASSO_1000 <- sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                    function(x) FPR(true, x))
MSE_ADMM_mvFMR_LASSO_500 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n500$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_LASSO_1000 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_LASSO_result$n1000$optimal, 
                                                    function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_LASSO_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_SCAD ##
TPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_500 <- sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                  function(x) FPR(true, x))
TPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) TPR(true, x))
FPR_ADMM_mvFMR_SCAD_1000 <- sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                   function(x) FPR(true, x))
MSE_ADMM_mvFMR_SCAD_500 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n500$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_SCAD_1000 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_SCAD_result$n1000$optimal, 
                                                   function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_SCAD_1000) <- c("pi", "Bk", "sigma")

## ADMM_mvFMR_MCP ##
TPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_500 <- sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n500$optimal, 
                                 function(x) FPR(true, x))
TPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) TPR(true, x))
FPR_ADMM_mvFMR_MCP_1000 <- sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                  function(x) FPR(true, x))
MSE_ADMM_mvFMR_MCP_500 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n500$optimal, 
                                                 function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_500) <- c("pi", "Bk", "sigma")
MSE_ADMM_mvFMR_MCP_1000 <- as.data.frame(t(sapply(sim2_mod6_ADMM_mvFMR_MCP_result$n1000$optimal, 
                                                  function(x) MSE(true, x))))
names(MSE_ADMM_mvFMR_MCP_1000) <- c("pi", "Bk", "sigma")

## Table ##
model6_table <- list(n500=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                               sd=as.data.frame(matrix(nrow=9, ncol=5))),
                     n1000=list(mean=as.data.frame(matrix(nrow=9, ncol=5)), 
                                sd=as.data.frame(matrix(nrow=9, ncol=5))))
colnames(model6_table$n500$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model6_table$n500$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                      "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model6_table$n500$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model6_table$n500$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                    "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model6_table$n1000$mean) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model6_table$n1000$mean) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                       "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model6_table$n1000$sd) <- c("TPR", "FPR", "MSE_pi_k", "MSE_B_k", "MSE_sigma_k")
rownames(model6_table$n1000$sd) <- c("Flexmix", "Oracle", "mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L",
                                     "PGM mvFMR-S", "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")

model6_table$n500$mean$TPR <- c(mean(TPR_flexmix_500), mean(TPR_oracle_500), mean(TPR_mvFMR_500),
                                mean(TPR_EM_PGM_mvFMR_LASSO_500), mean(TPR_ADMM_mvFMR_LASSO_500),
                                mean(TPR_EM_PGM_mvFMR_SCAD_500), mean(TPR_ADMM_mvFMR_SCAD_500),
                                mean(TPR_EM_PGM_mvFMR_MCP_500), mean(TPR_ADMM_mvFMR_MCP_500))
model6_table$n500$mean$FPR <- c(mean(FPR_flexmix_500), mean(FPR_oracle_500), mean(FPR_mvFMR_500),
                                mean(FPR_EM_PGM_mvFMR_LASSO_500), mean(FPR_ADMM_mvFMR_LASSO_500),
                                mean(FPR_EM_PGM_mvFMR_SCAD_500), mean(FPR_ADMM_mvFMR_SCAD_500),
                                mean(FPR_EM_PGM_mvFMR_MCP_500), mean(FPR_ADMM_mvFMR_MCP_500))
model6_table$n500$mean$MSE_pi_k <- c(mean(MSE_flexmix_500$pi), mean(MSE_oracle_500$pi), 
                                     mean(MSE_mvFMR_500$pi), mean(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                     mean(MSE_ADMM_mvFMR_LASSO_500$pi), mean(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                     mean(MSE_ADMM_mvFMR_SCAD_500$pi), mean(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                     mean(MSE_ADMM_mvFMR_MCP_500$pi))
model6_table$n500$mean$MSE_B_k <- c(mean(MSE_flexmix_500$Bk), mean(MSE_oracle_500$Bk), 
                                    mean(MSE_mvFMR_500$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                    mean(MSE_ADMM_mvFMR_LASSO_500$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                    mean(MSE_ADMM_mvFMR_SCAD_500$Bk), mean(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                    mean(MSE_ADMM_mvFMR_MCP_500$Bk))
model6_table$n500$mean$MSE_sigma_k <- c(mean(MSE_flexmix_500$sigma), mean(MSE_oracle_500$sigma), 
                                        mean(MSE_mvFMR_500$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                        mean(MSE_ADMM_mvFMR_LASSO_500$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                        mean(MSE_ADMM_mvFMR_SCAD_500$sigma), mean(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                        mean(MSE_ADMM_mvFMR_MCP_500$sigma))
model6_table$n500$sd$TPR <- c(sd(TPR_flexmix_500), sd(TPR_oracle_500), sd(TPR_mvFMR_500),
                              sd(TPR_EM_PGM_mvFMR_LASSO_500), sd(TPR_ADMM_mvFMR_LASSO_500),
                              sd(TPR_EM_PGM_mvFMR_SCAD_500), sd(TPR_ADMM_mvFMR_SCAD_500),
                              sd(TPR_EM_PGM_mvFMR_MCP_500), sd(TPR_ADMM_mvFMR_MCP_500))
model6_table$n500$sd$FPR <- c(sd(FPR_flexmix_500), sd(FPR_oracle_500), sd(FPR_mvFMR_500),
                              sd(FPR_EM_PGM_mvFMR_LASSO_500), sd(FPR_ADMM_mvFMR_LASSO_500),
                              sd(FPR_EM_PGM_mvFMR_SCAD_500), sd(FPR_ADMM_mvFMR_SCAD_500),
                              sd(FPR_EM_PGM_mvFMR_MCP_500), sd(FPR_ADMM_mvFMR_MCP_500))
model6_table$n500$sd$MSE_pi_k <- c(sd(MSE_flexmix_500$pi), sd(MSE_oracle_500$pi), 
                                   sd(MSE_mvFMR_500$pi), sd(MSE_EM_PGM_mvFMR_LASSO_500$pi),
                                   sd(MSE_ADMM_mvFMR_LASSO_500$pi), sd(MSE_EM_PGM_mvFMR_SCAD_500$pi),
                                   sd(MSE_ADMM_mvFMR_SCAD_500$pi), sd(MSE_EM_PGM_mvFMR_MCP_500$pi),
                                   sd(MSE_ADMM_mvFMR_MCP_500$pi))
model6_table$n500$sd$MSE_B_k <- c(sd(MSE_flexmix_500$Bk), sd(MSE_oracle_500$Bk), 
                                  sd(MSE_mvFMR_500$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_500$Bk),
                                  sd(MSE_ADMM_mvFMR_LASSO_500$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_500$Bk),
                                  sd(MSE_ADMM_mvFMR_SCAD_500$Bk), sd(MSE_EM_PGM_mvFMR_MCP_500$Bk),
                                  sd(MSE_ADMM_mvFMR_MCP_500$Bk))
model6_table$n500$sd$MSE_sigma_k <- c(sd(MSE_flexmix_500$sigma), sd(MSE_oracle_500$sigma), 
                                      sd(MSE_mvFMR_500$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_500$sigma),
                                      sd(MSE_ADMM_mvFMR_LASSO_500$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_500$sigma),
                                      sd(MSE_ADMM_mvFMR_SCAD_500$sigma), sd(MSE_EM_PGM_mvFMR_MCP_500$sigma),
                                      sd(MSE_ADMM_mvFMR_MCP_500$sigma))
model6_table$n1000$mean$TPR <- c(mean(TPR_flexmix_1000), mean(TPR_oracle_1000), mean(TPR_mvFMR_1000),
                                 mean(TPR_EM_PGM_mvFMR_LASSO_1000), mean(TPR_ADMM_mvFMR_LASSO_1000),
                                 mean(TPR_EM_PGM_mvFMR_SCAD_1000), mean(TPR_ADMM_mvFMR_SCAD_1000),
                                 mean(TPR_EM_PGM_mvFMR_MCP_1000), mean(TPR_ADMM_mvFMR_MCP_1000))
model6_table$n1000$mean$FPR <- c(mean(FPR_flexmix_1000), mean(FPR_oracle_1000), mean(FPR_mvFMR_1000),
                                 mean(FPR_EM_PGM_mvFMR_LASSO_1000), mean(FPR_ADMM_mvFMR_LASSO_1000),
                                 mean(FPR_EM_PGM_mvFMR_SCAD_1000), mean(FPR_ADMM_mvFMR_SCAD_1000),
                                 mean(FPR_EM_PGM_mvFMR_MCP_1000), mean(FPR_ADMM_mvFMR_MCP_1000))
model6_table$n1000$mean$MSE_pi_k <- c(mean(MSE_flexmix_1000$pi), mean(MSE_oracle_1000$pi), 
                                      mean(MSE_mvFMR_1000$pi), mean(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                      mean(MSE_ADMM_mvFMR_LASSO_1000$pi), mean(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                      mean(MSE_ADMM_mvFMR_SCAD_1000$pi), mean(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                      mean(MSE_ADMM_mvFMR_MCP_1000$pi))
model6_table$n1000$mean$MSE_B_k <- c(mean(MSE_flexmix_1000$Bk), mean(MSE_oracle_1000$Bk), 
                                     mean(MSE_mvFMR_1000$Bk), mean(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_LASSO_1000$Bk), mean(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_SCAD_1000$Bk), mean(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                     mean(MSE_ADMM_mvFMR_MCP_1000$Bk))
model6_table$n1000$mean$MSE_sigma_k <- c(mean(MSE_flexmix_1000$sigma), mean(MSE_oracle_1000$sigma), 
                                         mean(MSE_mvFMR_1000$sigma), mean(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_LASSO_1000$sigma), mean(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_SCAD_1000$sigma), mean(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                         mean(MSE_ADMM_mvFMR_MCP_1000$sigma))
model6_table$n1000$sd$TPR <- c(sd(TPR_flexmix_1000), sd(TPR_oracle_1000), sd(TPR_mvFMR_1000),
                               sd(TPR_EM_PGM_mvFMR_LASSO_1000), sd(TPR_ADMM_mvFMR_LASSO_1000),
                               sd(TPR_EM_PGM_mvFMR_SCAD_1000), sd(TPR_ADMM_mvFMR_SCAD_1000),
                               sd(TPR_EM_PGM_mvFMR_MCP_1000), sd(TPR_ADMM_mvFMR_MCP_1000))
model6_table$n1000$sd$FPR <- c(sd(FPR_flexmix_1000), sd(FPR_oracle_1000), sd(FPR_mvFMR_1000),
                               sd(FPR_EM_PGM_mvFMR_LASSO_1000), sd(FPR_ADMM_mvFMR_LASSO_1000),
                               sd(FPR_EM_PGM_mvFMR_SCAD_1000), sd(FPR_ADMM_mvFMR_SCAD_1000),
                               sd(FPR_EM_PGM_mvFMR_MCP_1000), sd(FPR_ADMM_mvFMR_MCP_1000))
model6_table$n1000$sd$MSE_pi_k <- c(sd(MSE_flexmix_1000$pi), sd(MSE_oracle_1000$pi), 
                                    sd(MSE_mvFMR_1000$pi), sd(MSE_EM_PGM_mvFMR_LASSO_1000$pi),
                                    sd(MSE_ADMM_mvFMR_LASSO_1000$pi), sd(MSE_EM_PGM_mvFMR_SCAD_1000$pi),
                                    sd(MSE_ADMM_mvFMR_SCAD_1000$pi), sd(MSE_EM_PGM_mvFMR_MCP_1000$pi),
                                    sd(MSE_ADMM_mvFMR_MCP_1000$pi))
model6_table$n1000$sd$MSE_B_k <- c(sd(MSE_flexmix_1000$Bk), sd(MSE_oracle_1000$Bk), 
                                   sd(MSE_mvFMR_1000$Bk), sd(MSE_EM_PGM_mvFMR_LASSO_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_LASSO_1000$Bk), sd(MSE_EM_PGM_mvFMR_SCAD_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_SCAD_1000$Bk), sd(MSE_EM_PGM_mvFMR_MCP_1000$Bk),
                                   sd(MSE_ADMM_mvFMR_MCP_1000$Bk))
model6_table$n1000$sd$MSE_sigma_k <- c(sd(MSE_flexmix_1000$sigma), sd(MSE_oracle_1000$sigma), 
                                       sd(MSE_mvFMR_1000$sigma), sd(MSE_EM_PGM_mvFMR_LASSO_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_LASSO_1000$sigma), sd(MSE_EM_PGM_mvFMR_SCAD_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_SCAD_1000$sigma), sd(MSE_EM_PGM_mvFMR_MCP_1000$sigma),
                                       sd(MSE_ADMM_mvFMR_MCP_1000$sigma))
