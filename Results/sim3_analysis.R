##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####            R-code for reproducing Table 6 in Section 5.3.            ####
##############################################################################

# The following R-code contains functions for calculating predictive log-likelihood,
# which is the twice negative log-likelihood.

# Each ".rda" file contains the optimal estimators('pi', 'Bk', 'sigma'),
# the optimal number of components 'K', and the probability 'density' 
# at each point for the selected optimal K.

### Data ###
source("./Data/simulation_data.R")

### Functions ###
source("./Scripts/Functions/simulation_summary_functions.R")

### Model 7 ###
## results ##
load("./Results/Section 5.3/Model 7/sim3_mod7_mvFMR_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_EM_PGM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_EM_PGM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_EM_PGM_mvFMR_MCP_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_ADMM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_ADMM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.3/Model 7/sim3_mod7_ADMM_mvFMR_MCP_result.rda")

pred_ll_mvFMR_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_LASSO_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_SCAD_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_MCP_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_LASSO_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_SCAD_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_MCP_500 <- numeric(length=500)
pred_ll_mvFMR_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_LASSO_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_SCAD_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_MCP_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_LASSO_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_SCAD_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_MCP_1000 <- numeric(length=500)
for(i in 1:500){
  pred_ll_mvFMR_500[i] <- predictive_ll(sim3_mod7_mvFMR_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_LASSO_500[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_LASSO_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_SCAD_500[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_SCAD_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_MCP_500[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_MCP_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_LASSO_500[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_SCAD_500[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_MCP_500[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_MCP_result$n500, i, 500)
  
  pred_ll_mvFMR_1000[i] <- predictive_ll(sim3_mod7_mvFMR_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_LASSO_1000[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_LASSO_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_SCAD_1000[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_SCAD_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_MCP_1000[i] <- predictive_ll(sim3_mod7_EM_PGM_mvFMR_MCP_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_LASSO_1000[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_LASSO_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_SCAD_1000[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_SCAD_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_MCP_1000[i] <- predictive_ll(sim3_mod7_ADMM_mvFMR_MCP_result$n1000, i, 1000)
}

## Table ##
model7_table <- list(n500=list(mean=as.data.frame(matrix(nrow=7, ncol=1)), 
                               sd=as.data.frame(matrix(nrow=7, ncol=1))),
                     n1000=list(mean=as.data.frame(matrix(nrow=7, ncol=1)), 
                                sd=as.data.frame(matrix(nrow=7, ncol=1))))
colnames(model7_table$n500$mean) <- c("pred_ll")
rownames(model7_table$n500$mean) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                      "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model7_table$n500$sd) <- c("pred_ll")
rownames(model7_table$n500$sd) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                      "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model7_table$n1000$mean) <- c("pred_ll")
rownames(model7_table$n1000$mean) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                      "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model7_table$n1000$sd) <- c("pred_ll")
rownames(model7_table$n1000$sd) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                    "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
model7_table$n500$mean$pred_ll <- c(mean(pred_ll_mvFMR_500), mean(pred_ll_EM_PGM_mvFMR_LASSO_500), 
                                    mean(pred_ll_ADMM_mvFMR_LASSO_500), mean(pred_ll_EM_PGM_mvFMR_SCAD_500),
                                    mean(pred_ll_ADMM_mvFMR_SCAD_500), mean(pred_ll_EM_PGM_mvFMR_MCP_500),
                                    mean(pred_ll_ADMM_mvFMR_MCP_500))
model7_table$n500$sd$pred_ll <- c(sd(pred_ll_mvFMR_500), sd(pred_ll_EM_PGM_mvFMR_LASSO_500), 
                                    sd(pred_ll_ADMM_mvFMR_LASSO_500), sd(pred_ll_EM_PGM_mvFMR_SCAD_500),
                                    sd(pred_ll_ADMM_mvFMR_SCAD_500), sd(pred_ll_EM_PGM_mvFMR_MCP_500),
                                    sd(pred_ll_ADMM_mvFMR_MCP_500))
model7_table$n1000$mean$pred_ll <- c(mean(pred_ll_mvFMR_1000), mean(pred_ll_EM_PGM_mvFMR_LASSO_1000), 
                                     mean(pred_ll_ADMM_mvFMR_LASSO_1000), mean(pred_ll_EM_PGM_mvFMR_SCAD_1000),
                                     mean(pred_ll_ADMM_mvFMR_SCAD_1000), mean(pred_ll_EM_PGM_mvFMR_MCP_1000),
                                     mean(pred_ll_ADMM_mvFMR_MCP_1000))
model7_table$n1000$sd$pred_ll <- c(sd(pred_ll_mvFMR_1000), sd(pred_ll_EM_PGM_mvFMR_LASSO_1000), 
                                   sd(pred_ll_ADMM_mvFMR_LASSO_1000), sd(pred_ll_EM_PGM_mvFMR_SCAD_1000),
                                   sd(pred_ll_ADMM_mvFMR_SCAD_1000), sd(pred_ll_EM_PGM_mvFMR_MCP_1000),
                                   sd(pred_ll_ADMM_mvFMR_MCP_1000))


### Model 8 ###
## results ##
load("./Results/Section 5.3/Model 8/sim3_mod8_mvFMR_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_EM_PGM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_EM_PGM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_EM_PGM_mvFMR_MCP_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_ADMM_mvFMR_LASSO_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_ADMM_mvFMR_SCAD_result.rda")
load("./Results/Section 5.3/Model 8/sim3_mod8_ADMM_mvFMR_MCP_result.rda")

pred_ll_mvFMR_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_LASSO_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_SCAD_500 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_MCP_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_LASSO_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_SCAD_500 <- numeric(length=500)
pred_ll_ADMM_mvFMR_MCP_500 <- numeric(length=500)
pred_ll_mvFMR_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_LASSO_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_SCAD_1000 <- numeric(length=500)
pred_ll_EM_PGM_mvFMR_MCP_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_LASSO_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_SCAD_1000 <- numeric(length=500)
pred_ll_ADMM_mvFMR_MCP_1000 <- numeric(length=500)
for(i in 1:500){
  pred_ll_mvFMR_500[i] <- predictive_ll(sim3_mod8_mvFMR_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_LASSO_500[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_LASSO_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_SCAD_500[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_SCAD_result$n500, i, 500)
  pred_ll_EM_PGM_mvFMR_MCP_500[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_MCP_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_LASSO_500[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_SCAD_500[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD_result$n500, i, 500)
  pred_ll_ADMM_mvFMR_MCP_500[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_MCP_result$n500, i, 500)
  
  pred_ll_mvFMR_1000[i] <- predictive_ll(sim3_mod8_mvFMR_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_LASSO_1000[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_LASSO_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_SCAD_1000[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_SCAD_result$n1000, i, 1000)
  pred_ll_EM_PGM_mvFMR_MCP_1000[i] <- predictive_ll(sim3_mod8_EM_PGM_mvFMR_MCP_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_LASSO_1000[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_LASSO_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_SCAD_1000[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_SCAD_result$n1000, i, 1000)
  pred_ll_ADMM_mvFMR_MCP_1000[i] <- predictive_ll(sim3_mod8_ADMM_mvFMR_MCP_result$n1000, i, 1000)
}

## Table ##
model8_table <- list(n500=list(mean=as.data.frame(matrix(nrow=7, ncol=1)), 
                               sd=as.data.frame(matrix(nrow=7, ncol=1))),
                     n1000=list(mean=as.data.frame(matrix(nrow=7, ncol=1)), 
                                sd=as.data.frame(matrix(nrow=7, ncol=1))))
colnames(model8_table$n500$mean) <- c("pred_ll")
rownames(model8_table$n500$mean) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                      "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model8_table$n500$sd) <- c("pred_ll")
rownames(model8_table$n500$sd) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                    "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model8_table$n1000$mean) <- c("pred_ll")
rownames(model8_table$n1000$mean) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                       "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
colnames(model8_table$n1000$sd) <- c("pred_ll")
rownames(model8_table$n1000$sd) <- c("mvFMR", "PGM mvFMR-L", "ADMM mvFMR-L", "PGM mvFMR-S", 
                                     "ADMM mvFMR-S", "PGM mvFMR-M", "ADMM mvFMR-M")
model8_table$n500$mean$pred_ll <- c(mean(pred_ll_mvFMR_500), mean(pred_ll_EM_PGM_mvFMR_LASSO_500), 
                                    mean(pred_ll_ADMM_mvFMR_LASSO_500), mean(pred_ll_EM_PGM_mvFMR_SCAD_500),
                                    mean(pred_ll_ADMM_mvFMR_SCAD_500), mean(pred_ll_EM_PGM_mvFMR_MCP_500),
                                    mean(pred_ll_ADMM_mvFMR_MCP_500))
model8_table$n500$sd$pred_ll <- c(sd(pred_ll_mvFMR_500), sd(pred_ll_EM_PGM_mvFMR_LASSO_500), 
                                  sd(pred_ll_ADMM_mvFMR_LASSO_500), sd(pred_ll_EM_PGM_mvFMR_SCAD_500),
                                  sd(pred_ll_ADMM_mvFMR_SCAD_500), sd(pred_ll_EM_PGM_mvFMR_MCP_500),
                                  sd(pred_ll_ADMM_mvFMR_MCP_500))
model8_table$n1000$mean$pred_ll <- c(mean(pred_ll_mvFMR_1000), mean(pred_ll_EM_PGM_mvFMR_LASSO_1000), 
                                     mean(pred_ll_ADMM_mvFMR_LASSO_1000), mean(pred_ll_EM_PGM_mvFMR_SCAD_1000),
                                     mean(pred_ll_ADMM_mvFMR_SCAD_1000), mean(pred_ll_EM_PGM_mvFMR_MCP_1000),
                                     mean(pred_ll_ADMM_mvFMR_MCP_1000))
model8_table$n1000$sd$pred_ll <- c(sd(pred_ll_mvFMR_1000), sd(pred_ll_EM_PGM_mvFMR_LASSO_1000), 
                                   sd(pred_ll_ADMM_mvFMR_LASSO_1000), sd(pred_ll_EM_PGM_mvFMR_SCAD_1000),
                                   sd(pred_ll_ADMM_mvFMR_SCAD_1000), sd(pred_ll_EM_PGM_mvFMR_MCP_1000),
                                   sd(pred_ll_ADMM_mvFMR_MCP_1000))
