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

mvFMR_oracle_fixedK <- function(X, Y, n=nrow(X), P=(ncol(X)-1), m=ncol(Y), K, true_B, maxiter=100){
  
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
    Bk[[2]] <- m_step_Bk_oracle(X, Y, Bk, inv_sigma, true_Bk, w, 1)
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
      Bk[[t+1]] <- m_step_Bk_oracle(X, Y, Bk, inv_sigma, true_Bk, w, t)
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
      
      output <- list(pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
      print(output)
      cat("Iteration :", t+1, "\n")
    
      if(t == maxiter){
        break
      }else{
        t <- t+1
      }
    }
  }, silent = TRUE)
  
  optimal_w <- w[[t]]
  optimal_density <- density[[t]]
  
  return(list(optimal=output, w=optimal_w, density=optimal_density))
}
