###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####       R-code for applying mvFMR-MCP to the diabetes diagnosis data        ####
####           using the EM-PGM algorithm presented in Section 6.1.            ####
###################################################################################

# The following R-code demonstrates how to implement mvFMR-MCP using the EM-PGM algorithm 
# and the ADMM solver to the diabetes diagnosis data presented in Section 6.1.

# But, it can be analyzed by using the other methods (e.g., mvFMR-LASSO, mvFMR-SCAD).

### Data ###
source("./Data/diabetes_diagnosis_data.R")

### Functions ###
source("./Scripts/Analysis/PGM_mvFMR_MCP.R")
source("./Scripts/Analysis/ADMM_mvFMR_MCP.R")

### Fitting ###
log_lambdas <- seq(from=log10(10^-4), to=log10(1), length.out=30)
lambda_s <- 10^log_lambdas

## mvFMR-MCP with the EM-PGM algorithm ##
PGM_diabetes_diagnosis <- PGM_mvFMR_MCP_nonfixedK(X, Y, eta_val=0.25, a=3.7, lambda=lambda_s)
save(PGM_diabetes_diagnosis, file="./Results/Section6.1/PGM_diabetes_diagnosis.rda")

## mvFMR-MCP with the ADMM solver ##
ADMM_diabetes_diagnosis <- ADMM_mvFMR_MCP_nonfixedK(X, Y, rho=1, a=3.7, lambda=lambda_s)
save(ADMM_diabetes_diagnosis, file="./Results/Section6.1/ADMM_diabetes_diagnosis.rda")