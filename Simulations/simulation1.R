###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####            R-codes to fit the simulation data in Section 5.1.             ####
###################################################################################

# The following R-code fits the simulation data presented in Section 5.1 using parallel computing.

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


### Model 1 ###
index_parallel <- matrix(seed_number_5.1.1, ncol=parallel_col, byrow=TRUE)

## Please set n1 by 500, 1000, or values you want. ##
n1 <- 500
# n1 <- 1000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                true <- data_generate_5.1.1(i, n1)$true
                                cluster <- data_generate_5.1.1(i, n1)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim1_mod1_flexmix_", n1, ".rda")
  result_name <- paste0("sim1_mod1_flexmix_", n1)
  
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
    file_name <- paste0("sim1_mod1_flexmix_", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                true_Bk <- data_generate_5.1.1(i, n1)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim1_mod1_oracle_", n1, ".rda")
  result_name <- paste0("sim1_mod1_oracle_", n1)
  
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
    file_name <- paste0("sim1_mod1_oracle_", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim1_mod1_mvFMR_", n1, ".rda")
  result_name <- paste0("sim1_mod1_mvFMR_", n1)
  
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
    file_name <- paste0("sim1_mod1_mvFMR_", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_LASSO_", n1, ".rda")
  result_name <- paste0("sim1_mod1_EM_PGM_mvFMR_LASSO_", n1)
  
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
    file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_LASSO_", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_SCAD", n1, ".rda")
  result_name <- paste0("sim1_mod1_EM_PGM_mvFMR_SCAD", n1)
  
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
    file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_SCAD", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_MCP", n1, ".rda")
  result_name <- paste0("sim1_mod1_EM_PGM_mvFMR_MCP", n1)
  
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
    file_name <- paste0("sim1_mod1_EM_PGM_mvFMR_MCP", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_ADMM_mvFMR_LASSO_", n1, ".rda")
  result_name <- paste0("sim1_mod1_ADMM_mvFMR_LASSO_", n1)
  
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
    file_name <- paste0("sim1_mod1_ADMM_mvFMR_LASSO_", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_ADMM_mvFMR_SCAD", n1, ".rda")
  result_name <- paste0("sim1_mod1_ADMM_mvFMR_SCAD", n1)
  
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
    file_name <- paste0("sim1_mod1_ADMM_mvFMR_SCAD", n1, ".rda")
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
                                X <- data_generate_5.1.1(i, n1)$X
                                Y <- data_generate_5.1.1(i, n1)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod1_ADMM_mvFMR_MCP", n1, ".rda")
  result_name <- paste0("sim1_mod1_ADMM_mvFMR_MCP", n1)
  
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
    file_name <- paste0("sim1_mod1_ADMM_mvFMR_MCP", n1, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}



### Model 2 ###
index_parallel <- matrix(seed_number_5.1.2, ncol=parallel_col, byrow=TRUE)

## Please set n2 by 500, 1000, or values you want. ##
n2 <- 500
# n2 <- 1000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                true <- data_generate_5.1.2(i, n2)$true
                                cluster <- data_generate_5.1.2(i, n2)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim1_mod2_flexmix_", n2, ".rda")
  result_name <- paste0("sim1_mod2_flexmix_", n2)
  
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
    file_name <- paste0("sim1_mod2_flexmix_", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                true_Bk <- data_generate_5.1.2(i, n2)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim1_mod2_oracle_", n2, ".rda")
  result_name <- paste0("sim1_mod2_oracle_", n2)
  
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
    file_name <- paste0("sim1_mod2_oracle_", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim1_mod2_mvFMR_", n2, ".rda")
  result_name <- paste0("sim1_mod2_mvFMR_", n2)
  
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
    file_name <- paste0("sim1_mod2_mvFMR_", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_LASSO", n2, ".rda")
  result_name <- paste0("sim1_mod2_EM_PGM_mvFMR_LASSO", n2)
  
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
    file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_LASSO", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_SCAD", n2, ".rda")
  result_name <- paste0("sim1_mod2_EM_PGM_mvFMR_SCAD", n2)
  
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
    file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_SCAD", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_MCP", n2, ".rda")
  result_name <- paste0("sim1_mod2_EM_PGM_mvFMR_MCP", n2)
  
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
    file_name <- paste0("sim1_mod2_EM_PGM_mvFMR_MCP", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_ADMM_mvFMR_LASSO_", n2, ".rda")
  result_name <- paste0("sim1_mod2_ADMM_mvFMR_LASSO_", n2)
  
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
    file_name <- paste0("sim1_mod2_ADMM_mvFMR_LASSO_", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_ADMM_mvFMR_SCAD", n2, ".rda")
  result_name <- paste0("sim1_mod2_ADMM_mvFMR_SCAD", n2)
  
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
    file_name <- paste0("sim1_mod2_ADMM_mvFMR_SCAD", n2, ".rda")
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
                                X <- data_generate_5.1.2(i, n2)$X
                                Y <- data_generate_5.1.2(i, n2)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod2_ADMM_mvFMR_MCP", n2, ".rda")
  result_name <- paste0("sim1_mod2_ADMM_mvFMR_MCP", n2)
  
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
    file_name <- paste0("sim1_mod2_ADMM_mvFMR_MCP", n2, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}


### Model 3 ###
index_parallel <- matrix(seed_number_5.1.3, ncol=parallel_col, byrow=TRUE)

## Please set n3 by 1000, 2000, or values you want. ##
n3 <- 1000
# n3 <- 2000

## Flexmix ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                true <- data_generate_5.1.3(i, n3)$true
                                cluster <- data_generate_5.1.3(i, n3)$cluster
                                flexmix_fixedK(X, Y, true=true, cluster=cluster, K=2)
                              }
  file_name <- paste0("sim1_mod3_flexmix_", n3, ".rda")
  result_name <- paste0("sim1_mod3_flexmix_", n3)
  
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
    file_name <- paste0("sim1_mod3_flexmix_", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                true_Bk <- data_generate_5.1.3(i, n3)$true$Bk
                                oracle_fixedK(X, Y, K=2, true_B=true_Bk)
                              }
  file_name <- paste0("sim1_mod3_oracle_", n3, ".rda")
  result_name <- paste0("sim1_mod3_oracle_", n3)
  
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
    file_name <- paste0("sim1_mod3_oracle_", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                mvFMR_fixedK(X, Y, K=2)
                              }
  file_name <- paste0("sim1_mod3_mvFMR_", n3, ".rda")
  result_name <- paste0("sim1_mod3_mvFMR_", n3)
  
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
    file_name <- paste0("sim1_mod3_mvFMR_", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                PGM_mvFMR_LASSO_fixedK(X, Y, K=2, eta=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_LASSO", n3, ".rda")
  result_name <- paste0("sim1_mod3_EM_PGM_mvFMR_LASSO", n3)
  
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
    file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_LASSO", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                PGM_mvFMR_SCAD_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_SCAD", n3, ".rda")
  result_name <- paste0("sim1_mod3_EM_PGM_mvFMR_SCAD", n3)
  
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
    file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_SCAD", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                PGM_mvFMR_MCP_fixedK(X, Y, K=2, eta=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_MCP", n3, ".rda")
  result_name <- paste0("sim1_mod3_EM_PGM_mvFMR_MCP", n3)
  
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
    file_name <- paste0("sim1_mod3_EM_PGM_mvFMR_MCP", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                ADMM_mvFMR_LASSO_fixedK(X, Y, K=2, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_ADMM_mvFMR_LASSO_", n3, ".rda")
  result_name <- paste0("sim1_mod3_ADMM_mvFMR_LASSO_", n3)
  
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
    file_name <- paste0("sim1_mod3_ADMM_mvFMR_LASSO_", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_ADMM_mvFMR_SCAD", n3, ".rda")
  result_name <- paste0("sim1_mod3_ADMM_mvFMR_SCAD", n3)
  
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
    file_name <- paste0("sim1_mod3_ADMM_mvFMR_SCAD", n3, ".rda")
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
                                X <- data_generate_5.1.3(i, n3)$X
                                Y <- data_generate_5.1.3(i, n3)$Y
                                ADMM_mvFMR_MCP_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim1_mod3_ADMM_mvFMR_MCP", n3, ".rda")
  result_name <- paste0("sim1_mod3_ADMM_mvFMR_MCP", n3)
  
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
    file_name <- paste0("sim1_mod3_ADMM_mvFMR_MCP", n3, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}
