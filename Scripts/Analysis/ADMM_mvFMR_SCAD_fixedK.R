###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####     R-code for applying mvFMR-SCAD to the simulated data with known K     ####
####          using the ADMM solver presented in Section 5.1 - 5.2.            ####
###################################################################################

# The following R-code demonstrates how to implement mvFMR-SCAD using the ADMM solver 
# to the simulated data from Section 5.1 - 5.2.

# As a reference, we assumed that the number of components K was known in Section 5.1 - 5.2.

### Example ###
# source("./Data/simulation_seed_number.R")
# source("./Data/simulation_data.R")
# source("./Functions/functions.R")
# X <- data_generate_5.1.1(seed_number_5.1.1[1], 500)$X
# Y <- data_generate_5.1.1(seed_number_5.1.1[1], 500)$Y
# admm_result_scad_fixedK <- ADMM_mvFMR_SCAD_fixedK(X, Y, K=2, rho=1, alpha=1.5, a=3.7, lambda=lambda_s)

source("./Data/simulation_seed_number.R")
source("./Data/simulation_data.R")
source("./Scripts/Functions/functions.R")

lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3,
              30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3,
              60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3,
              90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1,
              35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2,
              7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1)

ADMM_mvFMR_SCAD_fixedK <- function(X, Y, n=nrow(X), P=(ncol(X)-1), m=ncol(Y), e_abs=1e-6, e_rel=1e-4, K, rho, alpha, a, lambda, maxiter=100){
  
  BIC_SCAD <- vector(length=length(lambda))
  total_SCAD <- list(list())
  total_OUTPUTS_SCAD <- list(list())
  optimal_OUTPUTS_SCAD <- list()
  w_s <- list(list())
  optimal_w <- list()
  density_s <- list(list())
  optimal_density <- list()
  
  for(z in 1:length(lambda)){
    try({
      theta_diff <- vector()
      density <- list(list())
      w <- list(list())
      Bk <- list(list())
      sigma <- list(list())
      inv_sigma <- list(list())
      pi_ <- list(vector())
      
      # Initialization #
      Bk_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
      Bk[[1]] <- lapply(seq(dim(Bk_array)[3]), function(x) Bk_array[ , , x])
      
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
      
      sigma_array <- array(rep(0, m*m*K), c(m, m, K))
      sigma[[1]] <- lapply(lapply(seq(dim(sigma_array)[3]), function(x) sigma_array[ , , x]), function(x) x+diag(m))
      
      lamb <- rep(lambda[z], K)
      
      theta_diff[1] <- 0
      
      w_array <- array(rep(0, n*1*K), c(n, 1, K))
      w[[1]] <- lapply(seq(dim(w_array)[3]), function(x) w_array[ , , x])
      
      density[[1]] <- list()
      inv_sigma[[1]] <- list()
      for(k in 1:K){
        inv_sigma[[1]][[k]] <- solve(sigma[[1]][[k]])
      }
      for(k in 1:K){
        density[[1]][[k]] <- density_f(X, Y, Bk[[1]][[k]], sigma[[1]][[k]], inv_sigma[[1]][[k]])
        for(i in 1:n){
          if(density[[1]][[k]][i] == 0){
            density[[1]][[k]][i] <- 1e-300
          }
        }
      }
      
      # The first E and M steps #
      w[[2]] <- e_step_w(pi_, density, 1)
      
      pi_[[2]] <- m_step_pi(w, 1)
      Bk[[2]] <- m_step_Bk_ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, 1)
      sigma[[2]] <- m_step_sigma(X, Y, Bk, w, 1)
      
      density[[2]] <- list()
      inv_sigma[[2]] <- list()
      for(k in 1:K){
        inv_sigma[[2]][[k]] <- solve(sigma[[2]][[k]])
      }
      for(k in 1:K){
        density[[2]][[k]] <- density_f(X, Y, Bk[[2]][[k]], sigma[[2]][[k]], inv_sigma[[2]][[k]])
        for(i in 1:n){
          if(density[[2]][[k]][i] == 0){
            density[[2]][[k]][i] <- 1e-300
          }
        }
      }
      
      theta_diff[2] <- total_diff_norm(pi_, Bk, sigma, 1)
      
      # The iteration of E and M steps #
      t <- 2
      while(theta_diff[t] >= 1e-6){
        w[[t+1]] <- e_step_w(pi_, density, t)
        
        pi_[[t+1]] <- m_step_pi(w, t)
        Bk[[t+1]] <- m_step_Bk_ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, t)
        sigma[[t+1]] <- m_step_sigma(X, Y, Bk, w, t)
        
        density[[t+1]] <- list()
        inv_sigma[[t+1]] <- list()
        for(k in 1:K){
          inv_sigma[[t+1]][[k]] <- solve(sigma[[t+1]][[k]])
        }
        for(k in 1:K){
          density[[t+1]][[k]] <- density_f(X, Y, Bk[[t+1]][[k]], sigma[[t+1]][[k]], inv_sigma[[t+1]][[k]])
          for(i in 1:n){
            if(density[[t+1]][[k]][i] == 0){
              density[[t+1]][[k]][i] <- 1e-300
            }
          }
        }
        theta_diff[t+1] <- total_diff_norm(pi_, Bk, sigma, t)
        
        output <- list(lambda=lambda[z], pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
        print(output)
        cat("Iteration :", t+1, "\n")
        cat("Lambda :", lambda[z], "\n")
        cat("K :", K, "\n")
        
        if(t == maxiter){
          break
        }else{
          t <- t+1
        }
      }
    }, silent = TRUE)
    
    BIC_SCAD[z] <- modified_BIC_6(pi_, Bk, density, w, t)
    
    w_s[[z]] <- w[[t]]
    density_s[[z]] <- density[[t]]
    total_OUTPUTS_SCAD[[z]] <- output
  }
  
  min_SCAD <- which.min(replace(BIC_SCAD, which(is.infinite(BIC_SCAD)==T), NA))
  optimal_OUTPUTS_SCAD <- total_OUTPUTS_SCAD[[min_SCAD]]
  optimal_w <- w_s[[min_SCAD]]
  optimal_density <- density_s[[min_SCAD]]
  
  return(list(total=total_OUTPUTS_SCAD, optimal=optimal_OUTPUTS_SCAD, total_w=w_s, 
              total_density=density_s, optimal_w=optimal_w, optimal_density=optimal_density, BIC=BIC_SCAD))
}
