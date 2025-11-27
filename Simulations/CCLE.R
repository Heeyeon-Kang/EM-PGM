###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####             R-code for applying mvFMR-MCP to the CCLE dataset             ####
####           using the EM-PGM algorithm presented in Section 6.2.            ####
###################################################################################

# The following R-code demonstrates how to implement mvFMR-MCP using the EM-PGM algorithm 
# and the ADMM solver to the CCLE dataset presented in Section 6.2.

# But, it can be analyzed by using the other methods (e.g., mvFMR-LASSO, mvFMR-SCAD).

### Data ###
source("./Data/CCLE_data.R")

### Functions ###
source("./Scripts/Analysis/PGM_mvFMR_MCP.R")
source("./Scripts/Analysis/ADMM_mvFMR_MCP.R")

### Fitting ###
log_lambdas <- seq(from=log10(10^-2), to=log10(1), length.out=30)
lambda_s <- 10^log_lambdas

## mvFMR-MCP with the EM-PGM algorithm ##
PGM_CCLE <- PGM_mvFMR_MCP_nonfixedK(X, Y, eta_val=0.01, a=3.7, lambda=lambda_s)
save(PGM_CCLE, file="./Results/Section6.2/PGM_CCLE.rda")

## mvFMR-MCP with the ADMM solver ##
ADMM_CCLE <- ADMM_mvFMR_MCP_nonfixedK(X, Y, rho=0.5, a=3.7, lambda=lambda_s)
save(ADMM_CCLE, file="./Results/Section6.2/ADMM_CCLE.rda")