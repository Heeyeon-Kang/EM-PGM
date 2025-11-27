###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####          R-code for applying the methods to the simulated data.           ####
###################################################################################

# The following R-code includes all the necessary functions to apply methods 
# (e.g., Flexmix, Oracle, mvFMR, mvFMR-L, mvFMR-S, and mvFMR-M) to the simulated data 
# and reproduce all the results presented in Section 5. 


### Packages ###
requiredPackages <- c("MASS", "pracma", "expm", "causact", "fossil", "gtools", "flexmix")
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


### The density of multivariate Gaussian regression model ###
density_f <- function(X, Y, B, sigma, inv_sigma){
  n <- nrow(X)
  m <- ncol(Y)
  
  a <- (2*pi)^(-m/2)
  b <- determinant(sigma, logarithm = FALSE)$sign
  c <- determinant(sigma, logarithm = FALSE)$modulus[1]
  d <- (b*c)^(-1/2)
  e <- vapply(1:n, function(i){
    res <- Y[i,] - t(B)%*%X[i,]
    as.numeric((-1/2)*t(res)%*%inv_sigma%*%res)
  }, numeric(1))
  f <- a*d*exp(e)
  
  return(f)
}

safe_regularize_sigma <- function(S, alpha = 0.1, epsilon = 1e-3, prior = NULL){
  if (is.null(prior)) prior <- diag(ncol(S))
  if (!all(is.finite(S))) return(prior)
  
  S <- alpha*prior + (1-alpha)*S
  eig <- eigen(S, symmetric = TRUE)
  if (any(!is.finite(eig$values))) return(prior)
  
  D_safe <- pmax(eig$values, epsilon)
  S_psd <- eig$vectors%*%diag(D_safe)%*%t(eig$vectors)
  return(S_psd)
}

### PGM functions by penalty type ###
gradient_obj <- function(X, Y, wk, Bk, inv_sigma_k){
  n <- nrow(X)
  resid <- Y - X%*%Bk

  Z <- resid %*% inv_sigma_k
  Xw <- X*wk
  grad <- -(1/n)*crossprod(Xw, Z)
  
  return(grad)
}

PGM_mvFMR_LASSO <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  inv_sigma_old <- inv_sigma[[l]]
  w_new <- w[[l+1]]
  K <- length(pi_new)
  
  PGM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  
  # Update Bk #
  PGM_start_time <- as.numeric(proc.time()["elapsed"])
  Bk_sol[[2]] <- lapply(1:K, function(k){
    LASSO_thresholding(Bk_sol[[1]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma_old[[k]]), 
                       pi_new[k], lambda[k], eta_val)})
  PGM_end_time <- as.numeric(proc.time()["elapsed"])
  PGM_computation_time <- rbind(PGM_computation_time,
                                data.frame(iter_num = 1,
                                           time = PGM_end_time-PGM_start_time,
                                           acc_time = PGM_end_time-PGM_start_time))
  
  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]]) >= rep(1e-6, K))){
    PGM_start_time <- as.numeric(proc.time()["elapsed"])
    Bk_sol[[t+1]] <- lapply(1:K, function(k){
      LASSO_thresholding(Bk_sol[[t]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma_old[[k]]), 
                         pi_new[k], lambda[k], eta_val)})
    PGM_end_time <- as.numeric(proc.time()["elapsed"])
    PGM_computation_time <- rbind(PGM_computation_time,
                                  data.frame(iter_num = t,
                                             time = PGM_end_time-PGM_start_time,
                                             acc_time = PGM_computation_time[t-1,3]+PGM_end_time-PGM_start_time))

    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = PGM_computation_time))
}

PGM_mvFMR_SCAD <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  inv_sigma_old <- inv_sigma[[l]]
  w_new <- w[[l+1]]
  K <- length(pi_new)
  
  PGM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  
  # Update Bk #
  PGM_start_time <- as.numeric(proc.time()["elapsed"])
  Bk_sol[[2]] <- lapply(1:K, function(k){
    SCAD_thresholding(Bk_sol[[1]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma_old[[k]]), 
                       pi_new[k], lambda[k], eta_val, a)})
  PGM_end_time <- as.numeric(proc.time()["elapsed"])
  PGM_computation_time <- rbind(PGM_computation_time,
                                data.frame(iter_num = 1,
                                           time = PGM_end_time-PGM_start_time,
                                           acc_time = PGM_end_time-PGM_start_time))

  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]]) >= rep(1e-6, K))){
    # Update Bk #
    PGM_start_time <- as.numeric(proc.time()["elapsed"])
    Bk_sol[[t+1]] <- lapply(1:K, function(k){
      SCAD_thresholding(Bk_sol[[t]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma_old[[k]]), 
                         pi_new[k], lambda[k], eta_val, a)})
    PGM_end_time <- as.numeric(proc.time()["elapsed"])
    PGM_computation_time <- rbind(PGM_computation_time,
                                  data.frame(iter_num = t,
                                             time = PGM_end_time-PGM_start_time,
                                             acc_time = PGM_computation_time[t-1,3]+PGM_end_time-PGM_start_time))
    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = PGM_computation_time))
}

PGM_mvFMR_MCP <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  inv_sigma_old <- inv_sigma[[l]]
  w_new <- w[[l+1]]
  K <- length(pi_new)
  
  PGM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  
  # Update Bk #
  PGM_start_time <- as.numeric(proc.time()["elapsed"])
  Bk_sol[[2]] <- lapply(1:K, function(k){
    MCP_thresholding(Bk_sol[[1]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma_old[[k]]), 
                      pi_new[k], lambda[k], eta_val, a)})
  PGM_end_time <- as.numeric(proc.time()["elapsed"])
  PGM_computation_time <- rbind(PGM_computation_time,
                                data.frame(iter_num = 1,
                                           time = PGM_end_time-PGM_start_time,
                                           acc_time = PGM_end_time-PGM_start_time))

  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]]) >= rep(1e-6, K))){
    # Update Bk #
    PGM_start_time <- as.numeric(proc.time()["elapsed"])
    Bk_sol[[t+1]] <- lapply(1:K, function(k){
      MCP_thresholding(Bk_sol[[t]][[k]] - eta_val*gradient_obj(X, Y, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma_old[[k]]), 
                        pi_new[k], lambda[k], eta_val, a)})
    PGM_end_time <- as.numeric(proc.time()["elapsed"])
    PGM_computation_time <- rbind(PGM_computation_time,
                                  data.frame(iter_num = t,
                                             time = PGM_end_time-PGM_start_time,
                                             acc_time = PGM_computation_time[t-1,3]+PGM_end_time-PGM_start_time))

    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = PGM_computation_time))
}


### ADMM functions by penalty type ###
primal_residual <- function(B, D) B-D

dual_residual <- function(U, u, r) -r*(U-u)

epsilon_primal <- function(B, C, e_absolute, e_relative){
  e_abs_term <- sqrt(length(B)) * e_absolute
  e_rel_term <- max(norm(B,"F"), norm(C,"F")) * e_relative
  return(e_abs_term + e_rel_term)
}

epsilon_dual <- function(U, e_absolute, e_relative){
  e_abs_term <- sqrt(length(U)) * e_absolute
  e_rel_term <- norm(U,"F") * e_relative
  return(e_abs_term + e_rel_term)
}

bartels_stewart <- function(A, B, C){
  schur_A <- Schur(A)
  schur_B <- Schur(B)
  
  Q_A <- schur_A$Q
  T_A <- schur_A$T
  Q_B <- schur_B$Q
  T_B <- schur_B$T
  
  C_tilde <- t(Q_A)%*%C%*%Q_B
  
  nrw_A <- nrow(A)
  ncl_B <- ncol(B)
  Y <- matrix(0, nrw_A, ncl_B)
  
  for (j in ncl_B:1){
    rhs <- C_tilde[,j]
    if (j < ncl_B){
      rhs <- rhs - Y[,(j+1):ncl_B, drop=FALSE]%*%matrix(T_B[j,(j+1):ncl_B], nrow=ncl_B-j)
    }
    Y[,j] <- solve(T_A+T_B[j,j]*diag(nrw_A), rhs)
  }
  
  X <- Q_A%*%Y%*%t(Q_B)
  
  return(X)
}

ADMM_mvFMR_LASSO <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  K <- length(pi_new)
  n <- nrow(X)
  P <- ncol(X)-1
  m <- ncol(Y)
  
  Bk_sol <- list(list())
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  
  # R is list of primal residuals and S is list of dual residuals #
  R <- list(list())     
  S <- list(list())    
  
  ADMM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol[[1]] <- Bk[[l]]
  Ck[[1]] <- list()
  U1[[1]] <- list()
  Hk[[1]] <- list()
  for(k in 1:K){
    Ck[[1]][[k]] <- Bk_sol[[1]][[k]]
    U1[[1]][[k]] <- matrix(0, nrow = (P+1), ncol = m)
  }
  
  ADMM_start_time <- as.numeric(proc.time()["elapsed"])
  A <- replicate(K, matrix(0, P+1, P+1), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- lapply(1:K, function(k) rho*sigma_old[[k]])
  
  W <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma_old[[k]])
  
  # Update Bk by solving Sylvester's equation #
  Bk_sol[[2]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
  
  # Update Ck #
  Hk[[2]] <- lapply(1:K, function(k) alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]])
  Ck[[2]] <- lapply(1:K, function(k) LASSO_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_new[k], lambda[k], rho))
  
  # Update U1 #
  U1[[2]] <- lapply(1:K, function(k) U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]])

  ADMM_end_time <- as.numeric(proc.time()["elapsed"])
  ADMM_computation_time <- rbind(ADMM_computation_time,
                                data.frame(iter_num = 1,
                                           time = ADMM_end_time-ADMM_start_time,
                                           acc_time = ADMM_end_time-ADMM_start_time))
  
  # Calculating residuals #
  R[[2]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]]))
  S[[2]] <- lapply(1:K, function(k) dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho))

  epsilon[[2]] <- matrix(nrow = K, ncol = 2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs, e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  t <- 2
  while(stopping_criteria(R[[t]], S[[t]], epsilon[[t]])){
    ADMM_start_time <- as.numeric(proc.time()["elapsed"])
    G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma_old[[k]])

    # Update Bk by solving Sylvester's equation #
    Bk_sol[[t+1]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
    
    # Update Ck #
    Hk[[t+1]] <- lapply(1:K, function(k) alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]])
    Ck[[t+1]] <- lapply(1:K, function(k) LASSO_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_new[k], lambda[k], rho))
    
    # Update U1 #
    U1[[t+1]] <- lapply(1:K, function(k) U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]])
    
    ADMM_end_time <- as.numeric(proc.time()["elapsed"])
    ADMM_computation_time <- rbind(ADMM_computation_time,
                                   data.frame(iter_num = t,
                                              time = ADMM_end_time-ADMM_start_time,
                                              acc_time = ADMM_computation_time[t-1,3]+ADMM_end_time-ADMM_start_time))
    
    # Calculating residuals #
    R[[t+1]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]]))
    S[[t+1]] <- lapply(1:K, function(k) dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho))
    
    epsilon[[t+1]] <- matrix(nrow = K, ncol = 2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = ADMM_computation_time))
}

ADMM_mvFMR_SCAD <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  K <- length(pi_new)
  n <- nrow(X)
  P <- ncol(X)-1
  m <- ncol(Y)
  
  Bk_sol <- list(list())
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  
  # R is list of primal residuals and S is list of dual residuals #
  R <- list(list())     
  S <- list(list())    
  
  ADMM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol[[1]] <- Bk[[l]]
  Hk[[1]] <- list()
  Ck[[1]] <- Bk_sol[[1]]
  U1[[1]] <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
  
  ADMM_start_time <- as.numeric(proc.time()["elapsed"])
  A <- replicate(K, matrix(0, P+1, P+1), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- lapply(1:K, function(k) rho*sigma_old[[k]])
  
  W <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma_old[[k]])
  
  # Update Bk by solving Sylvester's equation #
  Bk_sol[[2]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
  
  # Update Ck #
  Hk[[2]] <- lapply(1:K, function(k) alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]])
  Ck[[2]] <- lapply(1:K, function(k) SCAD_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_new[k], lambda[k], rho, a))
  
  # Update U1 #
  U1[[2]] <- lapply(1:K, function(k) U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]])
  
  ADMM_end_time <- as.numeric(proc.time()["elapsed"])
  ADMM_computation_time <- rbind(ADMM_computation_time,
                                 data.frame(iter_num = 1,
                                            time = ADMM_end_time-ADMM_start_time,
                                            acc_time = ADMM_end_time-ADMM_start_time))
  
  # Calculating residuals #
  R[[2]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]]))
  S[[2]] <- lapply(1:K, function(k) dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho))
  
  epsilon[[2]] <- matrix(nrow = K, ncol = 2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs, e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  t <- 2
  while(stopping_criteria(R[[t]], S[[t]], epsilon[[t]])){
    ADMM_start_time <- as.numeric(proc.time()["elapsed"])
    G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma_old[[k]])
    
    # Update Bk by solving Sylvester's equation #
    Bk_sol[[t+1]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
    
    # Update Ck #
    Hk[[t+1]] <- lapply(1:K, function(k) alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]])
    Ck[[t+1]] <- lapply(1:K, function(k) SCAD_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_new[k], lambda[k], rho, a))
    
    # Update U1 #
    U1[[t+1]] <- lapply(1:K, function(k) U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]])
    
    ADMM_end_time <- as.numeric(proc.time()["elapsed"])
    ADMM_computation_time <- rbind(ADMM_computation_time,
                                   data.frame(iter_num = t,
                                              time = ADMM_end_time-ADMM_start_time,
                                              acc_time = ADMM_computation_time[t-1,3]+ADMM_end_time-ADMM_start_time))
    
    # Calculating residuals #
    R[[t+1]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]]))
    S[[t+1]] <- lapply(1:K, function(k) dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho))
    
    epsilon[[t+1]] <- matrix(nrow = K, ncol = 2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = ADMM_computation_time))
}

ADMM_mvFMR_MCP <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l){
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  K <- length(pi_new)
  n <- nrow(X)
  P <- (ncol(X))-1
  m <- ncol(Y)
  
  Bk_sol <- list(list())
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  
  # R is list of primal residuals and S is list of dual residuals #
  R <- list(list())     
  S <- list(list())    
  
  ADMM_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  # Initialization #
  Bk_sol[[1]] <- Bk[[l]]
  Hk[[1]] <- list()
  Ck[[1]] <- Bk_sol[[1]]
  U1[[1]] <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
  
  ADMM_start_time <- as.numeric(proc.time()["elapsed"])
  A <- replicate(K, matrix(0, P+1, P+1), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- lapply(1:K, function(k) rho*sigma_old[[k]])
  
  W <- replicate(K, matrix(0, P+1, m), simplify = FALSE)
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma_old[[k]])
  
  # Update Bk by solving Sylvester's equation #
  Bk_sol[[2]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
  
  # Update Ck #
  Hk[[2]] <- lapply(1:K, function(k) alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]])
  Ck[[2]] <- lapply(1:K, function(k) MCP_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_new[k], lambda[k], rho, a))

  # Update U1 #
  U1[[2]] <- lapply(1:K, function(k) U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]])
  
  # Calculating residuals #
  R[[2]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]]))
  S[[2]] <- lapply(1:K, function(k) dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho))
  
  epsilon[[2]] <- matrix(nrow = K, ncol = 2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs, e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  t <- 2
  while(stopping_criteria(R[[t]], S[[t]], epsilon[[t]])){
    ADMM_start_time <- as.numeric(proc.time()["elapsed"])
    
    G <- lapply(1:K, function(k) W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma_old[[k]])
    
    # Update Bk by solving Sylvester's equation #
    Bk_sol[[t+1]] <- lapply(1:K, function(k) bartels_stewart(A[[k]], E[[k]], G[[k]]))
    
    # Update Ck #
    Hk[[t+1]] <- lapply(1:K, function(k) alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]])
    Ck[[t+1]] <- lapply(1:K, function(k) MCP_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_new[k], lambda[k], rho, a))
    
    # Update U1 #
    U1[[t+1]] <- lapply(1:K, function(k) U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]])
    
    ADMM_end_time <- as.numeric(proc.time()["elapsed"])
    ADMM_computation_time <- rbind(ADMM_computation_time,
                                   data.frame(iter_num = t,
                                              time = ADMM_end_time-ADMM_start_time,
                                              acc_time = ADMM_computation_time[t-1,3]+ADMM_end_time-ADMM_start_time))
    
    # Calculating residuals #
    R[[t+1]] <- lapply(1:K, function(k) primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]]))
    S[[t+1]] <- lapply(1:K, function(k) dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho))
    
    epsilon[[t+1]] <- matrix(nrow = K, ncol = 2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(list(B = Bk_new, time = ADMM_computation_time))
}

## Stopping Criteria ##
stopping_criteria <- function(R, S, E){    
  # E[[1]] is a vector of primal residuals and E[[2]] is a vector of dual residuals #
  R_norm <- sapply(R, function(x) norm(x, "F"))            
  S_norm <- sapply(S, function(x) norm(x, "F"))
  
  return(any(any(R_norm > E[,1]), any(S_norm > E[,2])))
}


### Thresholding functions ###
LASSO_thresholding <- function(BB, p, lam, step_size){
  threshold <- step_size*p*lam
  Bk_new <- apply(BB, c(1, 2), function(x){
    if(abs(x) <= threshold){
      return(0)
    }else{
      return(x-sign(x)*threshold)
    }
  })
  return(Bk_new)
}

SCAD_thresholding <- function(BB, p, lam, step_size, a){
  threshold1 <- step_size*p*lam
  threshold2 <- (1+step_size*p)*lam
  threshold3 <- a*lam
  
  Bk_new <- apply(BB, c(1, 2), function(x){
    abs_x <- abs(x)
    if(abs_x <= threshold1){
      return(0)
    }else if(abs_x <= threshold2){
      return(x-sign(x)*threshold1)
    }else if (abs_x <= threshold3){
      return((x-sign(x)*(a/(a-1))*threshold1)/(1-step_size*p/(a-1)))
    }else{
      return(x)
    }
  })
  return(Bk_new)
}

MCP_thresholding <- function(BB, p, lam, step_size, a){
  threshold1 <- step_size*p*lam
  threshold2 <- a*lam
  
  Bk_new <- apply(BB, c(1, 2), function(x){
    abs_x <- abs(x)
    if(abs_x > threshold2){
      return(x)
    }else if(abs_x <= threshold1){
      return(0)
    }else{
      return((x-sign(x)*threshold1)/(1-step_size*p/a))
    }
  })
  return(Bk_new)
}

diff_norm <- function(B, B_new){
  K <- length(B)
  aa <- vapply(1:K, function(k) norm(B_new[[k]] - B[[k]], "F"), numeric(1))
  return(aa)
}


### EM algorithm ###
## E step ##
e_step_w <- function(pi_, density, l){
  density_old <- density[[l]]

  n <- length(density_old[[1]])
  K <- length(density_old)

  safe_log <- function(x) log(pmax(x, .Machine$double.xmin))
  
  log_wk <- sapply(1:K, function(k) log(pi_[[l]][k]) + safe_log(density_old[[k]]))
  
  m <- apply(log_wk, 1, max)
  w_mat <- exp(log_wk-m)
  w_mat <- sweep(w_mat, 1, rowSums(w_mat), "/")
  
  w_mat <- pmin(pmax(w_mat, 1e-3), 1-1e-3)
  w_mat <- sweep(w_mat, 1, rowSums(w_mat), "/")
  
  w_new <- lapply(seq_len(K), function(k) w_mat[, k])
  
  return(w_new)
}

## M step ##
# Update pi_k in M step #
m_step_pi <- function(w, l){
  w_new <- w[[l+1]]
  n <- length(w_new[[1]])
  
  pi_new <- (1/n) * sapply(w_new, sum)
  return(pi_new)
}

# Update B_k in M step #
m_step_Bk_mvFMR <- function(X, Y, Bk, w, l){
  Bk_new <- list()
  Bk_old <- Bk[[l]]
  w_new <- w[[l+1]]
  
  mvFMR_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  n <- nrow(X)
  K <- length(Bk_old)
  
  mvFMR_start_time <- as.numeric(proc.time()["elapsed"])
  
  WXX_array <- array(rep(0, nrow(Bk_old[[1]])*nrow(Bk_old[[1]])*K), c(nrow(Bk_old[[1]]), nrow(Bk_old[[1]]), K))
  WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXX[[k]] <- WXX[[k]] + w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }

  WXY_array <- array(rep(0, nrow(Bk_old[[1]])*ncol(Bk_old[[1]])*K), c(nrow(Bk_old[[1]]), ncol(Bk_old[[1]]), K))
  WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXY[[k]] <- WXY[[k]] + w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }

  Bk_new <- lapply(1:K, function(k) solve(WXX[[k]])%*%WXY[[k]])
  
  mvFMR_end_time <- as.numeric(proc.time()["elapsed"])
  mvFMR_computation_time <- rbind(mvFMR_computation_time,
                                   data.frame(iter_num = 1,
                                              time = mvFMR_end_time-mvFMR_start_time,
                                              acc_time = mvFMR_end_time-mvFMR_start_time))

  return(list(B = Bk_new, time = mvFMR_computation_time))
}

m_step_Bk_PGM_mvFMR_LASSO <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, l){
  Bk_new <- PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, l)
  return(Bk_new)
}

m_step_Bk_PGM_mvFMR_SCAD <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l){
  Bk_new <- PGM_mvFMR_SCAD(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l)
  return(Bk_new)
}

m_step_Bk_PGM_mvFMR_MCP <- function(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l){
  Bk_new <- PGM_mvFMR_MCP(X, Y, pi_, Bk, inv_sigma, w, lambda, eta_val, a, l)
  return(Bk_new)
}

m_step_Bk_ADMM_mvFMR_LASSO <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, l){
  Bk_new <- ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, l)
  return(Bk_new)
}

m_step_Bk_ADMM_mvFMR_SCAD <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l){
  Bk_new <- ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l)
  return(Bk_new)
}

m_step_Bk_ADMM_mvFMR_MCP <- function(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l){
  Bk_new <- ADMM_mvFMR_MCP(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lambda, rho, alpha, a, l)
  return(Bk_new)
}

m_step_Bk_oracle <- function(X, Y, Bk, inv_sigma, true_Bk, w, l) {
  oracle_computation_time <- data.frame(iter_num = integer(), time = numeric(), acc_time = numeric())
  
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  Bk_old <- Bk[[l]]
  K <- length(Bk_old)
  
  Bk_new <- vector("list", K)
  
  t0 <- as.numeric(proc.time()["elapsed"])
  
  for (k in seq_len(K)) {
    w_k <- as.numeric(w[[l+1]][[k]])
    if (length(w_k) != n) stop(sprintf("Length of w[[l+1]][[%d]] is not n", k))
    if (any(w_k < 0)) stop("w is negative")
    
    XtWX <- crossprod(X, w_k * X) + diag(1e-6, p)
    XtWY <- crossprod(X, w_k * Y)
    
    S_inv <- inv_sigma[[l]][[k]]
    if (!all(dim(S_inv) == c(q, q))) stop(sprintf("Error in inv_sigma[[l]][[%d]] size", k))
    S_inv <- 0.5 * (S_inv + t(S_inv))
    
    mask <- (true_Bk[[k]] != 0)
    if (!all(dim(mask) == c(p, q))) stop(sprintf("Error in true_Bk[[%d]] size", k))
    free_idx <- which(as.vector(mask))
    
    if (length(free_idx) == 0L) {
      Bk_new[[k]] <- matrix(0, nrow = p, ncol = q)
      next
    }
    
    A <- kronecker(S_inv, XtWX)
    rhs <- as.vector(XtWY%*%S_inv)
    
    A_ff   <- A[free_idx, free_idx, drop = FALSE]
    rhs_f  <- rhs[free_idx]
    
    b_f <- tryCatch({
      R <- chol(A_ff)
      backsolve(R, forwardsolve(t(R), rhs_f))
    }, error = function(e) {
      A_ff2 <- A_ff + diag(1e-8, nrow(A_ff))
      solve(A_ff2, rhs_f)
    })
    
    b_vec <- numeric(p * q)
    b_vec[free_idx] <- b_f
    B_hat <- matrix(b_vec, nrow = p, ncol = q)
    
    Bk_new[[k]] <- B_hat
  }
  
  t1 <- as.numeric(proc.time()["elapsed"])
  elapsed <- t1 - t0
  oracle_computation_time <- rbind(
    oracle_computation_time,
    data.frame(iter_num = 1, time = elapsed, acc_time = elapsed)
  )
  
  return(list(B = Bk_new, time = oracle_computation_time))
}

# Update Sigma_k in M step #
m_step_sigma <- function(X, Y, Bk, w, l){
  n <- nrow(X)
  m <- ncol(Y)

  Bk_new <- Bk[[l+1]]
  w_new <- w[[l+1]]
  K <- length(Bk_new)

  sigma_new <- list()
  sigma_sum <- replicate(K, matrix(0, m, m), simplify = FALSE)

  for (k in 1:K) {
    resid_mat <- sapply(1:n, function(i) {
      res <- Y[i, ] - t(Bk_new[[k]])%*%X[i, ]
      w_new[[k]][i]*(res%*%t(res))
    }, simplify = "array")

    sigma_sum[[k]] <- apply(resid_mat, c(1, 2), sum)
    sigma_new[[k]] <- sigma_sum[[k]]/sum(w_new[[k]])
  }

  return(sigma_new)
}

### MSE ###
total_diff_norm <- function(pi_, Bk, sigma, l){
  pi_new <- pi_[[l+1]]
  Bk_new <- Bk[[l+1]]
  sigma_new <- sigma[[l+1]]
  pi_old <- pi_[[l]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  
  K <- length(pi_new)
  
  pi_norm <- (pi_new - pi_old)^2
  Bk_norm <- sapply(1:K, function(k){
    sum((Bk_new[[k]] - Bk_old[[k]])^2)
  })
  sigma_norm <- sapply(1:K, function(k){
    sum((sigma_new[[k]] - sigma_old[[k]])^2)
  })
  
  return(sqrt(sum(pi_norm, Bk_norm, sigma_norm)))
}

### The modified BIC ###
modified_BIC <- function(pi_, Bk, density, w, l){
  pi_old <- pi_[[l]]
  w_old <- w[[l]]
  density_old <- density[[l]]
  Bk_old <- Bk[[l]]
  
  n <- length(w_old[[1]])
  K <- length(pi_old)
  
  logpi_term <- rep(0, K)
  for(k in 1:K){
    logpi_term[k] <- sum(w_old[[k]])*log(pi_old[k])
  }
  logdensity_term <- rep(0, K)
  for(k in 1:K){
    logdensity_term[k] <- sum(w_old[[k]] * lapply(density_old,log)[[k]])
  }
  negative_ll <- (-2) * (sum(logpi_term+logdensity_term))
  d_e <- K + (K-1) + sum(sapply(Bk_old, function(x) colSums(x != 0)))
  
  return(negative_ll + log(n)*d_e)
}

modified_BIC_6 <- function(pi_, Bk, density, w, l){
  floor_6 <- function(x) as.data.frame(trunc(x*10^5)/10^5)
  pi_old <- pi_[[l]]
  w_old <- w[[l]]
  density_old <- density[[l]]
  Bk_old <- Bk[[l]]
  
  n <- length(w_old[[1]])
  K <- length(pi_old)
  
  logpi_term <- sum(sapply(w_old, sum) * log(pi_old))
  logdensity_term <- 0
  for(k in 1:K){
    logdensity_term <- logdensity_term + sum(w_old[[k]] * lapply(density_old, log)[[k]])
  }
  negative_ll <- (-2) * (logpi_term + logdensity_term)
  d_e <- K + (K-1) + sum(sapply(lapply(Bk_old, floor_6), function(x) colSums(x != 0)))
  
  return(negative_ll + log(n)*d_e)
}

