###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####            R-codes to fit the simulation data in Section 5.3.             ####
###################################################################################

# The following R-code fits the simulation data presented in Section 5.3 using parallel computing.

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
source("./Scripts/Analysis/mvFMR.R")
source("./Scripts/Analysis/PGM_mvFMR_LASSO.R")
source("./Scripts/Analysis/PGM_mvFMR_SCAD.R")
source("./Scripts/Analysis/PGM_mvFMR_MCP.R")
source("./Scripts/Analysis/ADMM_mvFMR_LASSO.R")
source("./Scripts/Analysis/ADMM_mvFMR_SCAD.R")
source("./Scripts/Analysis/ADMM_mvFMR_MCP.R")

lambda_s <- c(1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 
              1e-1, 15*1e-2, 2e-1, 25*1e-2, 3e-1, 35*1e-2, 4e-1, 45*1e-2, 5e-1, 
              55*1e-2, 6e-1, 65*1e-2, 7e-1, 75*1e-2, 8e-1, 85*1e-2, 9e-1, 95*1e-2, 1)

numCores <- detectCores()
registerDoParallel(numCores)

factors_of_500 <- c(1, 2, 4, 5, 10, 20, 25, 50, 100, 125, 250, 500)
valid_factors <- factors_of_500[factors_of_500 < numCores]
parallel_col <- max(valid_factors)


### Model 7 ###
index_parallel <- matrix(seed_number_5.3.1, ncol=parallel_col, byrow=TRUE)

## Please set n7 by 500, 1000, or values you want. ##
## To generate the correlated predictors, please choose "corX" instead of "indepX". ##
n7 <- 500
# n7 <- 1000
indepX <- TRUE
# corX <- FALSE

## mvFMR ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indep)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                mvFMR_nonfixedK(X, Y)
                              }
  file_name <- paste0("sim3_mod7_mvFMR_", n7, ".rda")
  result_name <- paste0("sim3_mod7_mvFMR_", n7)
  
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
    file_name <- paste0("sim3_mod7_mvFMR_", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indep)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_LASSO_nonfixedK(X, Y, eta_val=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_PGM_mvFMR_LASSO_", n7, ".rda")
  result_name <- paste0("sim3_mod7_PGM_mvFMR_LASSO_", n7)
  
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
    file_name <- paste0("sim3_mod7_PGM_mvFMR_LASSO_", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_SCAD_nonfixedK(X, Y, eta_val=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_PGM_mvFMR_SCAD", n7, ".rda")
  result_name <- paste0("sim3_mod7_PGM_mvFMR_SCAD", n7)
  
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
    file_name <- paste0("sim3_mod7_PGM_mvFMR_SCAD", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_MCP_nonfixedK(X, Y, eta_val=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_PGM_mvFMR_MCP", n7, ".rda")
  result_name <- paste0("sim3_mod7_PGM_mvFMR_MCP", n7)
  
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
    file_name <- paste0("sim3_mod7_PGM_mvFMR_MCP", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_LASSO_nonfixedK(X, Y, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_ADMM_mvFMR_LASSO_", n7, ".rda")
  result_name <- paste0("sim3_mod7_ADMM_mvFMR_LASSO_", n7)
  
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
    file_name <- paste0("sim3_mod7_ADMM_mvFMR_LASSO_", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_SCAD_nonfixedK(X, Y, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_ADMM_mvFMR_SCAD", n7, ".rda")
  result_name <- paste0("sim3_mod7_ADMM_mvFMR_SCAD", n7)
  
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
    file_name <- paste0("sim3_mod7_ADMM_mvFMR_SCAD", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.1(i, n7, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_MCP_nonfixedK(X, Y, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod7_ADMM_mvFMR_MCP", n7, ".rda")
  result_name <- paste0("sim3_mod7_ADMM_mvFMR_MCP", n7)
  
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
    file_name <- paste0("sim3_mod7_ADMM_mvFMR_MCP", n7, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}


### Model 8 ###
index_parallel <- matrix(seed_number_5.3.2, ncol=parallel_col, byrow=TRUE)

## Please set n8 by 500, 1000, or values you want. ##
## To generate the correlated predictors, please choose "corX" instead of "indepX". ##
n8 <- 500
# n8 <- 1000
indepX <- TRUE
# corX <- FALSE

## mvFMR ##
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indep)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                mvFMR_nonfixedK(X, Y)
                              }
  file_name <- paste0("sim3_mod8_mvFMR_", n8, ".rda")
  result_name <- paste0("sim3_mod8_mvFMR_", n8)
  
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
    file_name <- paste0("sim3_mod8_mvFMR_", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indep)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_LASSO_nonfixedK(X, Y, eta_val=1, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_PGM_mvFMR_LASSO_", n8, ".rda")
  result_name <- paste0("sim3_mod8_PGM_mvFMR_LASSO_", n8)
  
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
    file_name <- paste0("sim3_mod8_PGM_mvFMR_LASSO_", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_SCAD_nonfixedK(X, Y, eta_val=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_PGM_mvFMR_SCAD", n8, ".rda")
  result_name <- paste0("sim3_mod8_PGM_mvFMR_SCAD", n8)
  
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
    file_name <- paste0("sim3_mod8_PGM_mvFMR_SCAD", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val as 1 when the predictors are independent
# and 0.25 when they are correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                PGM_mvFMR_MCP_nonfixedK(X, Y, eta_val=1, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_PGM_mvFMR_MCP", n8, ".rda")
  result_name <- paste0("sim3_mod8_PGM_mvFMR_MCP", n8)
  
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
    file_name <- paste0("sim3_mod8_PGM_mvFMR_MCP", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-LASSO using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_LASSO_nonfixedK(X, Y, rho=1, alpha=1.5, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_ADMM_mvFMR_LASSO_", n8, ".rda")
  result_name <- paste0("sim3_mod8_ADMM_mvFMR_LASSO_", n8)
  
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
    file_name <- paste0("sim3_mod8_ADMM_mvFMR_LASSO_", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-SCAD using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_SCAD_nonfixedK(X, Y, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_ADMM_mvFMR_SCAD", n8, ".rda")
  result_name <- paste0("sim3_mod8_ADMM_mvFMR_SCAD", n8)
  
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
    file_name <- paste0("sim3_mod8_ADMM_mvFMR_SCAD", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}

## mvFMR-MCP using the ADMM solver ##
# We used rho as 1 in both predictor designs; independent and correlated.
for(j in 1:nrow(index_parallel)){
  parallel_outputs <- foreach(i=index_parallel[j,],
                              .combine=list,
                              .multicombine=TRUE) %dopar% {
                                dat <- data_generate_5.3.2(i, n8, indepX)
                                X <- dat$train_X
                                Y <- dat$train_Y
                                ADMM_mvFMR_MCP_nonfixedK(X, Y, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)
                              }
  file_name <- paste0("sim3_mod8_ADMM_mvFMR_MCP", n8, ".rda")
  result_name <- paste0("sim3_mod8_ADMM_mvFMR_MCP", n8)
  
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
    file_name <- paste0("sim3_mod8_ADMM_mvFMR_MCP", n8, ".rda")
    save(list=result_name, file=file.path("./Results", file_name))
  }
  rm(parallel_outputs)
  rm(file_name)
  rm(result_name)
}