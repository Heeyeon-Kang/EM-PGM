###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####                      R-code for reproducing Table 12                      ####
####         of the diabetes diagnosis data presented in Section 6.1.          ####
###################################################################################

# The following R-code is designed to reproduce Table 12 from Section 6.1, 
# based on the results of analyzing the diabetes diagnosis data using mvFMR-MCP 
# with EM-PGM algorithm.

### Packages ###
requiredPackages <- c("ggplot2")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

### Data ###
source("./Data/diabetes_diagnosis_data.R")

### Result ###
load("./Results/Section6.1/PGM_diabetes_diagnosis.rda")
load("./Results/Section6.1/ADMM_diabetes_diagnosis.rda")

### Analysis ###
diab <- PGM_diabetes_diagnosis
# diab <- ADMM_diabetes_diagnosis

optimal_K <- which.min(sapply(diab$BIC, min)) + 1

## correlation ##
corr <- vector(length=optimal_K)
for(i in 1:length(corr)){
  corr[i] <- diab$optimal[[optimal_K-1]]$sigma[[i]][1,2] / 
    sqrt(diab$optimal[[optimal_K-1]]$sigma[[i]][1,1] * 
           diab$optimal[[optimal_K-1]]$sigma[[i]][2,2])
}

## coefficients ##
sd_X <- c(intercept = 1, apply(diabetes_X_1, 2, sd), gender = 1, apply(diabetes_X_2, 2, sd))
new_Bk <- list()
for(k in 1:optimal_K){
  new_Bk[[k]] <- apply(diab$optimal[[optimal_K-1]]$Bk[[k]], 2, function(x) x/sd_X)
}