###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####          R-code for applying mvFMR-LASSO to the simulated data            ####
####             using the ADMM solver presented in Section 5.3.               ####
###################################################################################

# The following R-code demonstrates how to implement mvFMR-LASSO using the ADMM solver 
# to the simulated data from Section 5.3.

# As a reference, we assumed that the number of components K was unknown in Section 5.3.

### Example ###
# source("./Data/simulation_seed_number.R")
# source("./Data/simulation_data.R")
# source("./Scripts/Functions/functions.R")
# X <- data_generate_5.3.1(seed_number_5.3.1[1], 500)$train_X
# Y <- data_generate_5.3.1(seed_number_5.3.1[1], 500)$train_Y
# admm_result_lasso_nonfixedK <- ADMM_mvFMR_LASSO_nonfixedK(X, Y, rho=1, alpha=1.5, lambda=lambda_s)

source("./Data/simulation_seed_number.R")
source("./Data/simulation_data.R")
source("./Scripts/Functions/functions.R")

lambda_s <- c(1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 
              1e-1, 15*1e-2, 2e-1, 25*1e-2, 3e-1, 35*1e-2, 4e-1, 45*1e-2, 5e-1, 
              55*1e-2, 6e-1, 65*1e-2, 7e-1, 75*1e-2, 8e-1, 85*1e-2, 9e-1, 95*1e-2, 1)

ADMM_mvFMR_LASSO_nonfixedK <- function(X, Y, n = nrow(X), P = (ncol(X)-1), m = ncol(Y), init = TRUE,
                                       e_abs = 1e-6, e_rel = 1e-4, rho, alpha, lambda, maxiter = 100){
  
  BIC_LASSO_s <- list()
  BIC_LASSO <- vector(length = length(lambda))
  total_LASSO <- list(list())
  total_OUTPUTS_LASSO <- list(list())
  optimal_OUTPUTS_LASSO <- list()
  w_s <- list(list())
  optimal_w <- list()
  density_s <- list(list())
  optimal_density <- list()
  optimal_Bk_time_s <- list()
  optimal_EM_time_s <- list()
  
  EM_computation_time <- list()
  Bk_computation_time <- list()
  col_names <- c("iter_num", "time", "acc_time")
  
  for(K in 1:6){
    EM_computation_time[[K]] <- lapply(1:length(lambda), function(i){
      setNames(data.frame(
        integer(), numeric(), numeric(),
        stringsAsFactors = FALSE
      ), col_names)
    })
    Bk_computation_time[[K]] <- lapply(1:length(lambda), function(i){
      list()  
    })
    for(z in 1:length(lambda)){
      try({
        theta_diff <- vector()
        density <- list(list())
        w <- list(list())
        Bk <- list(list())
        sigma <- list(list())
        inv_sigma <- list(list())
        pi_ <- list(vector())
        
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
        
        lamb <- rep(lambda[z], K)
        
        theta_diff[1] <- 0
        
        w[[1]] <- replicate(K, rep(1/K, n), simplify = FALSE)
        
        ## The first E and M steps ##
        # E step #
        EM_start_time <- as.numeric(proc.time()["elapsed"])
        
        inv_sigma[[1]] <- lapply(sigma[[1]], function(s) solve(s))
        density[[1]] <- lapply(1:K, function(k){
          dens_k <- density_f(X, Y, Bk[[1]][[k]], sigma[[1]][[k]], inv_sigma[[1]][[k]])
          dens_k[dens_k == 0] <- 1e-300
          return(dens_k)
        })
        
        w[[2]] <- e_step_w(pi_, density, 1)
        
        # M step #
        pi_[[2]] <- m_step_pi(w, 1)
        
        new_Bk <- m_step_Bk_ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, 1)
        Bk[[2]] <- new_Bk$B
        
        sigma[[2]] <- m_step_sigma(X, Y, Bk, w, 1)
        sigma[[2]] <- lapply(1:K, function(k) safe_regularize_sigma(sigma[[2]][[k]]))
        
        EM_end_time <- as.numeric(proc.time()["elapsed"])
        EM_computation_time[[K]][[z]] <- rbind(EM_computation_time[[K]][[z]], 
                                               data.frame(iter_num = 1,
                                                          time = EM_end_time-EM_start_time, 
                                                          acc_time = EM_end_time-EM_start_time))
        Bk_computation_time[[K]][[z]][[1]] <- new_Bk$time
        
        theta_diff[2] <- total_diff_norm(pi_, Bk, sigma, 1)
        
        inv_sigma[[2]] <- lapply(sigma[[2]], function(s) solve(s))
        density[[2]] <- lapply(1:K, function(k){
          dens_k <- density_f(X, Y, Bk[[2]][[k]], sigma[[2]][[k]], inv_sigma[[2]][[k]])
          dens_k[dens_k == 0] <- 1e-300
          return(dens_k)
        })
        
        ## The iteration of E and M steps ##
        t <- 2
        while(theta_diff[t] >= 1e-6){
          # E step #
          EM_start_time <- as.numeric(proc.time()["elapsed"])
          
          w[[t+1]] <- e_step_w(pi_, density, t)
          
          # M step #
          pi_[[t+1]] <- m_step_pi(w, t)
          
          new_Bk <- m_step_Bk_ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, t)
          Bk[[t+1]] <- new_Bk$B
          
          sigma[[t+1]] <- m_step_sigma(X, Y, Bk, w, t)
          sigma[[t+1]] <- lapply(1:K, function(k) safe_regularize_sigma(sigma[[t+1]][[k]]))
          
          inv_sigma[[t+1]] <- lapply(sigma[[t+1]], function(s) solve(s))
          density[[t+1]] <- lapply(1:K, function(k){
            dens_k <- density_f(X, Y, Bk[[t+1]][[k]], sigma[[t+1]][[k]], inv_sigma[[t+1]][[k]])
            dens_k[dens_k == 0] <- 1e-300
            return(dens_k)
          })
          
          EM_end_time <- as.numeric(proc.time()["elapsed"])
          EM_computation_time[[K]][[z]] <- rbind(EM_computation_time[[K]][[z]], 
                                                 data.frame(iter_num = t,
                                                            time = EM_end_time-EM_start_time, 
                                                            acc_time = EM_computation_time[[K]][[z]]$acc_time[t-1]+EM_end_time-EM_start_time))
          Bk_computation_time[[K]][[z]][[t]] <- new_Bk$time
          
          theta_diff[t+1] <- total_diff_norm(pi_, Bk, sigma, t)
          
          output <- list(lambda=lambda[z], pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
          print(output)
          cat("Iteration :", t+1, "\n")
          cat("Lambda :", lambda[z], "\n")
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
      
      BIC_LASSO[z] <- modified_BIC_6(pi_, Bk, density, w, t)
      
      w_s[[z]] <- w[[t]]
      density_s[[z]] <- density[[t]]
      total_OUTPUTS_LASSO[[z]] <- output
    }
    
    BIC_LASSO_s[[K]] <- replace(BIC_LASSO, which(is.infinite(BIC_LASSO) == T), NA)
    min_LASSO <- which.min(BIC_LASSO_s[[K]])
    total_LASSO[[K]] <- total_OUTPUTS_LASSO
    optimal_OUTPUTS_LASSO[[K]] <- total_OUTPUTS_LASSO[[min_LASSO]]
    optimal_w[[K]] <- w_s[[min_LASSO]]
    optimal_density[[K]] <- density_s[[min_LASSO]]
    optimal_Bk_time_s[[K]] <- Bk_computation_time[[K]][[min_LASSO]]
    optimal_EM_time_s[[K]] <- EM_computation_time[[K]][[min_LASSO]]
  }
  
  optimal_BIC_LASSO_s <- sapply(BIC_LASSO_s, function(x) min(x, na.rm = TRUE))
  optimal_K <- which.min(optimal_BIC_LASSO_s)
  optimal_OUTPUTS_LASSO_K <- optimal_OUTPUTS_LASSO[[optimal_K]]
  optimal_w_K <- optimal_w[[optimal_K]]
  optimal_density_K <- optimal_density[[optimal_K]]
  optimal_Bk_time <- optimal_Bk_time_s[[optimal_K]]
  optimal_EM_time <- optimal_EM_time_s[[optimal_K]]
  
  return(list(total = total_LASSO, total_w = optimal_w, total_density = optimal_density,
              total_EM_time = EM_computation_time, total_Bk_time = Bk_computation_time,
              K = optimal_K, optimal_w = optimal_w_K, BIC = BIC_LASSO_s, 
              optimal = optimal_OUTPUTS_LASSO_K, optimal_density = optimal_density_K, optimal_BIC = optimal_BIC_LASSO_s,
              optimal_EM_time = optimal_EM_time, optimal_Bk_time = optimal_Bk_time))
}