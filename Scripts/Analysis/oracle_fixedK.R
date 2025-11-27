###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####      R-code of estimating the oracle estimator to the simulated data      ####
####               with known K presented in Section 5.1 - 5.2.                ####
###################################################################################

# The following R-code demonstrates how to implement finding oracle estimator to the 
# simulated data from Section 5.1 - 5.2.

# As a reference, we assumed that the number of components K was known in Section 5.1 - 5.2.

### Example ###
# source("./Data/simulation_seed_number.R")
# source("./Data/simulation_data.R")
# source("./Scripts/Functions/functions.R")
# X <- data_generate_5.1.1(seed_number_5.1.1[1], 500)$X
# Y <- data_generate_5.1.1(seed_number_5.1.1[1], 500)$Y
# true <- data_generate_5.1.1(seed_number_5.1.1[1], 500)$true
# true_Bk <- true$Bk
# result_oracle_fixedK <- mvFMR_oracle_fixedK(X, Y, K=2, true_B=true_Bk)

source("./Data/simulation_seed_number.R")
source("./Data/simulation_data.R")
source("./Scripts/Functions/functions.R")

mvFMR_oracle_fixedK <- function(X, Y, n = nrow(X), P = (ncol(X)-1), m = ncol(Y), init = TRUE, K, true_Bk, maxiter = 100){
 
  EM_computation_time <- list()
  Bk_computation_time <- list()
  col_names <- c("iter_num", "time", "acc_time")
  
  EM_computation_time <- as.data.frame(matrix(NA, nrow = 1, ncol = length(col_names)))
  colnames(EM_computation_time) <- col_names
  
  try({  
    theta_diff <- vector()
    density <- list(list())
    w <- list(list())
    Bk <- list(list())
    sigma <- list(list())
    inv_sigma <- list(list())
    pi_ <- list(vector(length=K))
    optimal_w <- list()
    optimal_density <- list()
    
    ## Initialization ##
    Bk[[1]] <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
    
    # If the estimation does not change, try init=FALSE #
    if(init == TRUE){
      pi_[[1]] <- rep(1/K, K)
    }else{
      if(K == 1){
        pi_[[1]] <- 1
      }else if(K == 2){
        pi_[[1]] <- c(0.45, 0.55)
      }else if(K == 3){
        pi_[[1]] <- c(0.27, 0.33, 0.4)
      }else if(K == 4){
        pi_[[1]] <- c(0.23, 0.24, 0.26, 0.27)
      }else if(K == 5){
        pi_[[1]] <- c(0.16, 0.18, 0.2, 0.22, 0.24)
      }else if(K == 6){
        pi_[[1]] <- c(0.14, 0.15, 0.16, 0.17, 0.18, 0.2)
      }
    }
    
    sigma[[1]] <- replicate(K, diag(m), simplify = FALSE)
    
    theta_diff[1] <- 0
    
    w[[1]] <- replicate(K, rep(1/K, n), simplify = FALSE)
    
    EM_start_time <- as.numeric(proc.time()["elapsed"])
    
    inv_sigma[[1]] <- lapply(sigma[[1]], function(s) solve(s))
    density[[1]] <- lapply(1:K, function(k) {
      dens_k <- density_f(X, Y, Bk[[1]][[k]], sigma[[1]][[k]], inv_sigma[[1]][[k]])
      dens_k[dens_k == 0] <- 1e-300
      return(dens_k)
    })
    
    ## The first E and M steps ##
    w[[2]] <- e_step_w(pi_, density, 1)
    
    pi_[[2]] <- m_step_pi(w, 1)
    
    new_Bk <- m_step_Bk_oracle(X, Y, Bk, inv_sigma, true_Bk, w, 1)
    Bk[[2]] <- new_Bk$B
    
    sigma[[2]] <- m_step_sigma(X, Y, Bk, w, 1)
    sigma[[2]] <- lapply(1:K, function(k) safe_regularize_sigma(sigma[[2]][[k]]))
    
    EM_end_time <- as.numeric(proc.time()["elapsed"])
    EM_computation_time[1,] <- data.frame(iter_num = 1,
                                          time = EM_end_time-EM_start_time, 
                                          acc_time = EM_end_time-EM_start_time)
    Bk_computation_time[[1]] <- new_Bk$time
    
    theta_diff[2] <- total_diff_norm(pi_, Bk, sigma, 1)
    
    inv_sigma[[2]] <- lapply(sigma[[2]], function(s) solve(s))
    density[[2]] <- lapply(1:K, function(k) {
      dens_k <- density_f(X, Y, Bk[[2]][[k]], sigma[[2]][[k]], inv_sigma[[2]][[k]])
      dens_k[dens_k == 0] <- 1e-300
      return(dens_k)
    })
    
    ## The iteration of E and M steps ##
    t <- 2
    while(theta_diff[t] >= 1e-6){
      EM_start_time <- as.numeric(proc.time()["elapsed"])
      
      w[[t+1]] <- e_step_w(pi_, density, t)
      
      pi_[[t+1]] <- m_step_pi(w, t)
      
      new_Bk <- m_step_Bk_oracle(X, Y, Bk, inv_sigma, true_Bk, w, t)
      Bk[[t+1]] <- new_Bk$B
      
      sigma[[t+1]] <- m_step_sigma(X, Y, Bk, w, t)
      sigma[[t+1]] <- lapply(1:K, function(k) safe_regularize_sigma(sigma[[t+1]][[k]]))
      
      inv_sigma[[t+1]] <- lapply(sigma[[t+1]], function(s) solve(s))
      density[[t+1]] <- lapply(1:K, function(k) {
        dens_k <- density_f(X, Y, Bk[[t+1]][[k]], sigma[[t+1]][[k]], inv_sigma[[t+1]][[k]])
        dens_k[dens_k == 0] <- 1e-300
        return(dens_k)
      })
      
      EM_end_time <- as.numeric(proc.time()["elapsed"])
      EM_computation_time <- rbind(EM_computation_time, 
                                   data.frame(iter_num = t, 
                                              time = EM_end_time-EM_start_time, 
                                              acc_time = EM_computation_time[t-1,3]+EM_end_time-EM_start_time))
      Bk_computation_time[[t]] <- new_Bk$time
      
      theta_diff[t+1] <- total_diff_norm(pi_, Bk, sigma, t)
      
      output <- list(pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
      print(output)
      cat("Iteration :", t+1, "\n")
      cat("K :", K, "\n")
      
      if(theta_diff[t+1] == Inf){
        break
      }else{
        if(t == maxiter){
          break
        }else{
          t <- t+1
        }
      }
    }
  }, silent = TRUE)
  
  optimal_w <- w[[t]]
  optimal_density <- density[[t]]
  
  return(list(total_density = density, total_w = w, 
              EM_time = EM_computation_time, Bk_time = Bk_computation_time,
              optimal_density = optimal_density, optimal = output, optimal_w = optimal_w))
}