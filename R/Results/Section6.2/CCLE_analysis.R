###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####                      R-code for reproducing Table 13                      ####
####   of Cancer Cell Line Encyclopedia (CCLE) data presented in Section 6.2.  ####
###################################################################################

# The following R-code is designed to reproduce Table 13 from Section 6.2, 
# based on the results of analyzing the CCLE data applied mvFMR-MCP with EM-PGM algorithm.

### Data ###
source("./Data/CCLE_data.R")

### Result ###
load("./Results/Section6.2/PGM_CCLE.rda")
load("./Results/Section6.2/ADMM_CCLE.rda")

### Analysis ###
CCLE <- PGM_CCLE
# CCLE <- ADMM_CCLE

optimal_K <- which.min(sapply(CCLE$BIC, min)) + 1

## correlation ##
corr <- list()
for(k in 1:optimal_K){
  corr[[k]] <- as.data.frame(diag(3))
  colnames(corr[[k]]) <- c("Erlotinib", "AZD6244", "PD-0325901")
  rownames(corr[[k]]) <- c("Erlotinib", "AZD6244", "PD-0325901")
  
  corr[[k]][1,2] <- CCLE$optimal[[optimal_K-1]]$sigma[[k]][1,2] / 
    sqrt(CCLE$optimal[[optimal_K-1]]$sigma[[k]][1,1]*CCLE$optimal[[optimal_K-1]]$sigma[[k]][2,2])
  corr[[k]][2,1] <- corr[[k]][1,2]
  
  corr[[k]][1,3] <- CCLE$optimal[[optimal_K-1]]$sigma[[k]][1,3] / 
    sqrt(CCLE$optimal[[optimal_K-1]]$sigma[[k]][1,1]*CCLE$optimal[[optimal_K-1]]$sigma[[k]][3,3])
  corr[[k]][3,1] <- corr[[k]][1,3]
  
  corr[[k]][2,3] <- CCLE$optimal[[optimal_K-1]]$sigma[[k]][2,3] / 
    sqrt(CCLE$optimal[[optimal_K-1]]$sigma[[k]][2,2]*CCLE$optimal[[optimal_K-1]]$sigma[[k]][3,3])
  corr[[k]][3,2] <- corr[[k]][2,3]
}

## coefficients ##
new_Bk <- CCLE$optimal[[optimal_K-1]]$Bk
for(k in 1:optimal_K){
  for(j in 1:length(X_sd)){
    new_Bk[[k]][j+1,] <- CCLE$optimal[[optimal_K-1]]$Bk[[k]][j+1,] / X_sd[j]
  }
}

for(k in 1:optimal_K){
  rownames(new_Bk[[k]]) <- c("intercept", top_genes)
}

active_Bk <- list()
for(k in 1:optimal_K){
  active_Bk[[k]] <- new_Bk[[k]][rowSums(abs(new_Bk[[k]]) > 1e-6, na.rm = TRUE) > 0, , drop = FALSE]
}