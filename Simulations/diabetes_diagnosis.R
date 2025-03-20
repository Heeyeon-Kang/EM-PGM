##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####     R-code for applying mvFMR-MCP to the diabetes diagnosis data     ####
####         using the EM-PGM algorithm presented in Section 6.1.         ####
##############################################################################

# The following R-code demonstrates how to implement mvFMR-MCP using 
# the EM-PGM algorithm to the diabetes diagnosis data presented in Section 6.1.

# But, it can be analyzed by using the other methods (e.g., mvFMR-LASSO, mvFMR-SCAD).

### Data ###
source("./Data/diabetes_diagnosis_data.R")

### Functions ###
source("./Scripts/Analysis/PGM_mvFMR_MCP.R")

### Fitting ###
lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
              30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
              60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
              90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
              35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
              7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 

diabetes_diagnosis_result <- PGM_mvFMR_MCP_nonfixedK(X, Y, eta=0.5, a=3.7, lambda=lambda_s, maxiter=50)

save(diabetes_diagnosis_result, file="./Results/diabetes_diagnosis_result.rda")
