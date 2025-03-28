###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####            R-codes to fit the simulation data in Section 5.2.             ####
###################################################################################

# The following R-code fits the simulation data presented in Section 5.2 using parallel computing.

# We performed 500 simulated trials using parallel computing with 50 cores.
# The following R-codes are designed to automatically adjust and use a smaller number 
# of cores based on the available cores of the user's computer.


### Packages ###
requiredPackages <- c("parallel", "foreach", "doParallel")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}


### Scripts ###
source("./Scripts/Analysis/flexmix_fixedK.R")
source("./Scripts/Analysis/oracle_fixedK.R")
source("./Scripts/Analysis/mvFMR_fixedK.R")
source("./Scripts/Analysis/PGM_mvFMR_LASSO_fixedK.R")
source("./Scripts/Analysis/PGM_mvFMR_SCAD_fixedK.R")
source("./Scripts/Analysis/PGM_mvFMR_MCP_fixedK.R")
source("./Scripts/Analysis/ADMM_mvFMR_LASSO_fixedK.R")
source("./Scripts/Analysis/ADMM_mvFMR_SCAD_fixedK.R")
source("./Scripts/Analysis/ADMM_mvFMR_MCP_fixedK.R")

lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
              30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
              60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
              90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
              35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
              7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 

numCores <- detectCores()
registerDoParallel(numCores)

factors_of_500 <- c(1, 2, 4, 5, 10, 20, 25, 50, 100, 125, 250, 500)
valid_factors <- factors_of_500[factors_of_500 < numCores]
parallel_col <- max(valid_factors)


### Model 4 ###
index_parallel <- matrix(seed_number_5.2.1, ncol=parallel_col, byrow=TRUE)

## Please set n4 by 500, 1000, or values you want. ##
n4 <- 500
# n4 <- 1000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                true <- data_generate_5.2.1(i, n4)$true
                                cluster <- data_generate_5.2.1(i, n4)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim2_mod4_flexmix_", n4, ".rda")
  result_name <- paste0("sim2_mod4_flexmix_", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_flexmix_", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## Oracle ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                true_Bk <- data_generate_5.2.1(i, n4)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim2_mod4_oracle_", n4, ".rda")
  result_name <- paste0("sim2_mod4_oracle_", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_oracle_", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim2_mod4_mvFMR_", n4, ".rda")
  result_name <- paste0("sim2_mod4_mvFMR_", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_mvFMR_", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_LASSO_", n4, ".rda")
  result_name <- paste0("sim2_mod4_EM_PGM_mvFMR_LASSO_", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_LASSO_", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_SCAD", n4, ".rda")
  result_name <- paste0("sim2_mod4_EM_PGM_mvFMR_SCAD", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_SCAD", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_MCP", n4, ".rda")
  result_name <- paste0("sim2_mod4_EM_PGM_mvFMR_MCP", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_EM_PGM_mvFMR_MCP", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_ADMM_mvFMR_LASSO_", n4, ".rda")
  result_name <- paste0("sim2_mod4_ADMM_mvFMR_LASSO_", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_ADMM_mvFMR_LASSO_", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_ADMM_mvFMR_SCAD", n4, ".rda")
  result_name <- paste0("sim2_mod4_ADMM_mvFMR_SCAD", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_ADMM_mvFMR_SCAD", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.1(i, n4)$X
                                Y <- data_generate_5.2.1(i, n4)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod4_ADMM_mvFMR_MCP", n4, ".rda")
  result_name <- paste0("sim2_mod4_ADMM_mvFMR_MCP", n4)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod4_ADMM_mvFMR_MCP", n4, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}


### Model 5 ###
index_parallel <- matrix(seed_number_5.2.2, ncol=parallel_col, byrow=TRUE)

## Please set n5 by 500, 1000, or values you want. ##
n5 <- 500
# n5 <- 1000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                true <- data_generate_5.2.2(i, n5)$true
                                cluster <- data_generate_5.2.2(i, n5)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim2_mod5_flexmix_", n5, ".rda")
  result_name <- paste0("sim2_mod5_flexmix_", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_flexmix_", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## Oracle ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                true_Bk <- data_generate_5.2.2(i, n5)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim2_mod5_oracle_", n5, ".rda")
  result_name <- paste0("sim2_mod5_oracle_", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_oracle_", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim2_mod5_mvFMR_", n5, ".rda")
  result_name <- paste0("sim2_mod5_mvFMR_", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_mvFMR_", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_LASSO_", n5, ".rda")
  result_name <- paste0("sim2_mod5_EM_PGM_mvFMR_LASSO_", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_LASSO_", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_SCAD", n5, ".rda")
  result_name <- paste0("sim2_mod5_EM_PGM_mvFMR_SCAD", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_SCAD", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_MCP", n5, ".rda")
  result_name <- paste0("sim2_mod5_EM_PGM_mvFMR_MCP", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_EM_PGM_mvFMR_MCP", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_ADMM_mvFMR_LASSO_", n5, ".rda")
  result_name <- paste0("sim2_mod5_ADMM_mvFMR_LASSO_", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_ADMM_mvFMR_LASSO_", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_ADMM_mvFMR_SCAD", n5, ".rda")
  result_name <- paste0("sim2_mod5_ADMM_mvFMR_SCAD", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_ADMM_mvFMR_SCAD", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.2(i, n5)$X
                                Y <- data_generate_5.2.2(i, n5)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod5_ADMM_mvFMR_MCP", n5, ".rda")
  result_name <- paste0("sim2_mod5_ADMM_mvFMR_MCP", n5)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod5_ADMM_mvFMR_MCP", n5, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}


### Model 6 ###
index_parallel <- matrix(seed_number_5.2.3, ncol=parallel_col, byrow=TRUE)

## Please set n6 by 500, 1000, or values you want. ##
n6 <- 500
# n6 <- 1000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                true <- data_generate_5.2.3(i, n6)$true
                                cluster <- data_generate_5.2.3(i, n6)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim2_mod6_flexmix_", n6, ".rda")
  result_name <- paste0("sim2_mod6_flexmix_", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_flexmix_", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## Oracle ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                true_Bk <- data_generate_5.2.3(i, n6)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim2_mod6_oracle_", n6, ".rda")
  result_name <- paste0("sim2_mod6_oracle_", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_oracle_", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim2_mod6_mvFMR_", n6, ".rda")
  result_name <- paste0("sim2_mod6_mvFMR_", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_mvFMR_", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_LASSO_", n6, ".rda")
  result_name <- paste0("sim2_mod6_EM_PGM_mvFMR_LASSO_", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_LASSO_", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_SCAD", n6, ".rda")
  result_name <- paste0("sim2_mod6_EM_PGM_mvFMR_SCAD", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_SCAD", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the EM-PGM algorithm ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_MCP", n6, ".rda")
  result_name <- paste0("sim2_mod6_EM_PGM_mvFMR_MCP", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_EM_PGM_mvFMR_MCP", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_ADMM_mvFMR_LASSO_", n6, ".rda")
  result_name <- paste0("sim2_mod6_ADMM_mvFMR_LASSO_", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_ADMM_mvFMR_LASSO_", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_ADMM_mvFMR_SCAD", n6, ".rda")
  result_name <- paste0("sim2_mod6_ADMM_mvFMR_SCAD", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_ADMM_mvFMR_SCAD", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the ADMM solver ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.2.3(i, n6)$X
                                Y <- data_generate_5.2.3(i, n6)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim2_mod6_ADMM_mvFMR_MCP", n6, ".rda")
  result_name <- paste0("sim2_mod6_ADMM_mvFMR_MCP", n6)
  
  if(file_name %in% list.files(path="./Results")){
    load(file=paste0("./Results/", file_name))
    loaded_result <- get(result_name)
    for(p in 1:length(parallel_outputs)){
      loaded_result[[length(loaded_result)+1]] <- parallel_outputs[[p]]
    }
    assign(result_name, loaded_result)
    save(list=result_name, file=file.path("./Results", file_name))
  }else{
    assign(result_name, parallel_outputs)
    file_name <- paste0("sim2_mod6_ADMM_mvFMR_MCP", n6, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}