###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####   R-code for applying mvFMR-MCP to the age-specific mortality rates data  #### 
####            using the EM-PGM algorithm presented in Section 6.2.           ####
###################################################################################

# The following R-code demonstrates how to implement mvFMR-MCP using the EM-PGM algorithm 
# to the age-specific mortality rates data presented in Section 6.2.

# But, it can be analyzed by using the other methods (e.g., mvFMR-LASSO, mvFMR-SCAD).

### Data ###
source("./Data/age_specific_mortality_rates_data.R")

### Functions ###
source("./Scripts/Analysis/PGM_mvFMR_MCP.R")

### Fitting ###
lambda_s <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
              0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1,
              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2,
              2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.5, 5)

age_specific_mortality_rates_result <- PGM_mvFMR_MCP_nonfixedK(X, Y, eta=0.3, a=3.7, lambda=lambda_s, maxiter=50)

save(age_specific_mortality_rates_result, file="./Results/Section 6.2/age_specific_mortality_rates_result.rda")
