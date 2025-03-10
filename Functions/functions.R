##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####           R-code of functions used for EM-PGM algorithm.             ####
##############################################################################

# The following R-code contains all functions for running EM-PGM algorithm.

### Packages install ###
requiredPackages <- c("MASS", "pracma", "expm", "corrcoverage", 
                      "causact", "fossil", "flexmix", "gtools")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

### The density of multivariate Gaussian regression model ###
density_f <- function(X, Y, B, sigma_k, inv_sigma){
  
  n <- nrow(X)
  m <- ncol(Y)
  
  a <- (2*pi)^(-m/2)
  b <- determinant(sigma_k, logarithm=FALSE)$sign
  c <- determinant(sigma_k, logarithm=FALSE)$modulus[1]
  d <- (b*c)^(-1/2)
  e <- vector(length=n)
  f <- vector(length=n)
  for(i in 1:n){
    e[i] <- (-1/2)*(t(Y[i,]-t(B)%*%X[i,])%*%inv_sigma%*%(Y[i,]-t(B)%*%X[i,]))
    f[i] <- a*d*exp(e[i])
  }
  return(f)
}


### PGM functions by penalty type ###
gradient_obj <- function(n, wk, B, inv_sigma){
  
  a <- list()
  for(i in 1:n){
    a[[i]] <- wk[i]*(X[i,]%*%t(Y[i,]-t(B)%*%X[i,])%*%inv_sigma)
  }
  b <- (-1/n)*Reduce("+",a)
  
  return(b)
}

PGM_mvFMR_LASSO <- function(pi_, Bk, sigma, w, lambda, eta, n, K, l){
  
  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  # Initialization
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  Bk_sol[[2]] <- list()
  inv_sigma <- list()
  for(k in 1:K){
    inv_sigma[[k]] <- solve(sigma_old[[k]])  
  }
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- LASSO_thresholding(Bk_sol[[1]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma[[k]]), 
                                           pi_new[k], lambda[k], eta)
  }
  
  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]], K) >= rep(1e-6, K))){
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- LASSO_thresholding(Bk_sol[[t]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma[[k]]), 
                                               pi_new[k], lambda[k], eta)
    }
    
    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(Bk_new)
}

PGM_mvFMR_SCAD <- function(pi_, Bk, sigma, w, lambda, eta, a, n, K, l){

  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  # Initialization
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  Bk_sol[[2]] <- list()
  inv_sigma <- list()
  for(k in 1:K){
    inv_sigma[[k]] <- solve(sigma_old[[k]])
  }
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- SCAD_thresholding(Bk_sol[[1]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma[[k]]),
                                          pi_new[k], lambda[k], eta, a)
  }

  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]], K) >= rep(1e-6, K))){
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- SCAD_thresholding(Bk_sol[[t]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma[[k]]),
                                              pi_new[k], lambda[k], eta, a)
    }

    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(Bk_new)
}

PGM_mvFMR_MCP <- function(pi_, Bk, sigma, w, lambda, eta, a, n, K, l){

  pi_new <- pi_[[l+1]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  w_new <- w[[l+1]]
  
  # Initialization
  Bk_sol <- list()
  Bk_sol[[1]] <- Bk_old
  Bk_sol[[2]] <- list()
  inv_sigma <- list()
  for(k in 1:K){
    inv_sigma[[k]] <- solve(sigma_old[[k]])
  }
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- MCP_thresholding(Bk_sol[[1]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[1]][[k]], inv_sigma[[k]]),
                                         pi_new[k], lambda[k], eta, a)
  }

  t <- 2
  while(all(diff_norm(Bk_sol[[t-1]], Bk_sol[[t]], K) >= rep(1e-6,K))){
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- MCP_thresholding(Bk_sol[[t]][[k]]-eta*gradient_obj(n, w_new[[k]], Bk_sol[[t]][[k]], inv_sigma[[k]]),
                                             pi_new[k], lambda[k], eta, a)
    }

    if(t == 500){
      break
    }else{
      t <- t+1
    }
  }
  
  Bk_new <- Bk_sol[[t]]
  
  return(Bk_new)
}

diff_norm <- function(B, B_new, K){
  aa <- rep(0,K)
  for(k in 1:K){
    aa[k] <- norm(B_new[[k]]-B[[k]], "F")
  }
  return(aa)
}


### Thresholding functions ###
LASSO_thresholding <- function(BB, p, lam, eta){

  Bk_new <- matrix(0, nrow=nrow(BB), ncol=ncol(BB))
  for(r in 1:nrow(Bk_new)){
    for(s in 1:ncol(Bk_new)){
      if(abs(BB[r,s]) <= eta*p*lam){
        Bk_new[r,s] <- 0
      }else{
        Bk_new[r,s] <- BB[r,s] - sign(BB[r,s])*eta*p*lam
      }
    }
  }
  return(Bk_new)
}

SCAD_thresholding <- function(BB, p, lam, eta, a){
  
  Bk_new <- matrix(0, nrow=nrow(BB), ncol=ncol(BB))
  for(r in 1:nrow(Bk_new)){
    for(s in 1:ncol(Bk_new)){
      if(abs(BB[r,s]) <= eta*p*lam){
        Bk_new[r,s] <- 0
      }else if(abs(BB[r,s]) <= (1+eta*p)*lam){
        Bk_new[r,s] <- BB[r,s] - sign(BB[r,s])*eta*p*lam
      }else if(abs(BB[r,s]) <= a*lam){
        Bk_new[r,s] <- (BB[r,s]-sign(BB[r,s])*(a/(a-1))*eta*p*lam) / (1-eta*p/(a-1))
      }else{
        Bk_new[r,s] <- BB[r,s]
      }
    }
  }
  return(Bk_new)
}

MCP_thresholding <- function(BB, p, lam, eta, a){
  Bk_new <- matrix(0, nrow=nrow(BB), ncol=ncol(BB))
  for(r in 1:nrow(Bk_new)){
    for(s in 1:ncol(Bk_new)){
      if(abs(BB[r,s]) > a*lam){
        Bk_new[r,s] <- BB[r,s]
      }else{
        if(abs(BB[r,s]) <= eta*p*lam){
          Bk_new[r,s] <- 0
        }else{
          Bk_new[r,s] <- (BB[r,s]-sign(BB[r,s])*eta*p*lam) / (1-eta*p/a)
        }
      }
    }
  }
  return(Bk_new)
}


### EM algorithm ###
## E step ##
e_step_w <- function(pi_, density, n, K, l){
  density_old <- density[[l]]
  w_new <- list()
  wk <- matrix(nrow=n, ncol=K)
  for(k in 1:K){
    wk[,k] <- pi_[[l]][k] * density_old[[k]]
  }
  for(k in 1:K){
    w_new[[k]] <- prop.table(wk, margin=1)[,k]
  }
  return(w_new)
}

## M step ##
m_step_pi <- function(w, n, l){
  w_new <- w[[l+1]]
  pi_new <- vector()
  sum_w <- sapply(w_new, sum)
  
  pi_new <- (1/n) * sum_w
  return(pi_new)
}

## Update B_k in M step ##
m_step_Bk_mvFMR <- function(X, Y, Bk, w, K, l){
  Bk_new <- list()
  Bk_old <- Bk[[l]]
  w_new <- w[[l+1]]
  
  WXX_array <- array(rep(0, nrow(Bk_old)*nrow(Bk_old)*K), c(nrow(Bk_old), nrow(Bk_old), K))
  WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXX[[k]] <- WXX[[k]] + w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }

  WXY_array <- array(rep(0, nrow(Bk_old)*ncol(Bk_old)*K), c(nrow(Bk_old), ncol(Bk_old), K))
  WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXY[[k]] <- WXY[[k]] + w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }

  for(k in 1:K){
    Bk_new[[k]] <- solve(WXX[[k]])%*%WXY[[k]]
  }

  return(Bk_new)
}

m_step_Bk_mvFMR_LASSO <- function(pi_, Bk, sigma, w, lambda, eta, n, K, l){
  Bk_new <- list()
  Bk_new <- PGM_mvFMR_LASSO(pi_, Bk, sigma, w, lambda, eta, n, K, l)
  
  return(Bk_new)
}

m_step_Bk_mvFMR_SCAD <- function(pi_, Bk, sigma, w, lambda, eta, a, n, K, l){
  Bk_new <- list()
  Bk_new <- PGM_mvFMR_SCAD(pi_, Bk, sigma, w, lambda, eta, a, n, K, l)

  return(Bk_new)
}

m_step_Bk_mvFMR_MCP <- function(pi_, Bk, sigma, w, lambda, eta, a, n, K, l){
  Bk_new <- list()
  Bk_new <- PGM_mvFMR_MCP(pi_, Bk, sigma, w, lambda, eta, a, n, K, l)

  return(Bk_new)
}

m_step_Bk_oracle <- function(X, Y, Bk, true_Bk, w, n, K, l){
  
  Bk_new <- list()
  Bk_old <- Bk[[l]]
  w_new <- w[[l+1]]
  
  WXX_array <- array(rep(0, nrow(Bk_old)*nrow(Bk_old)*K), c(nrow(Bk_old), nrow(Bk_old), K))
  WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXX[[k]] <- WXX[[k]] + w_new[[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  WXY_array <- array(rep(0, nrow(Bk_old)*ncol(Bk_old)*K), c(nrow(Bk_old), ncol(Bk_old), K))
  WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXY[[k]] <- WXY[[k]] + w_new[[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  Bk_oracle <- list()
  for(k in 1:K){
    Bk_new[[k]] <- solve(WXX[[k]])%*%WXY[[k]]
    Bk_oracle[[k]] <- Bk_new[[k]]
    Bk_oracle[[k]][which(true_Bk[[k]] == 0)] <- 0
    Bk_new[[k]] <- Bk_oracle[[k]]
  }
  
  return(Bk_new)
}

# m_step_Bk_oracle <- function(l){
#   Bk[[l+1]] <- list()
# 
#   WXX_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
#   WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
#   for(k in 1:K){
#     for(i in 1:n){
#       WXX[[k]] <- WXX[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
#     }
#   }
# 
#   WXY_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
#   WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
#   for(k in 1:K){
#     for(i in 1:n){
#       WXY[[k]] <- WXY[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
#     }
#   }
# 
#   Bk_oracle <- list()
#   for(k in 1:K){
#     Bk[[l+1]][[k]] <- ginv(WXX[[k]])%*%WXY[[k]]
#     Bk_oracle[[k]] <- Bk[[l+1]][[k]]
#     Bk_oracle[[k]][which(true_Bk[[k]] == 0)] <- 0
#     Bk[[l+1]][[k]] <- Bk_oracle[[k]]
#   }
# 
#   return(Bk[[l+1]])
# }

## Update Sigma_k in M step ##
m_step_sigma <- function(X, Y, Bk, sigma, w, K, l){
  
  n <- nrow(X)
  m <- ncol(Y)
  
  Bk_new <- Bk[[l+1]]
  w_new <- w[[l+1]]
  
  sigma_new <- list()
  sigma_sum_array <- array(rep(0, m*m*K), c(m, m, K))
  sigma_sum <- lapply(seq(dim(sigma_sum_array)[3]), function(x) sigma_sum_array[ , , x])
  
  for(k in 1:K){
    for(i in 1:n){
      sigma_sum[[k]] <- sigma_sum[[k]] + w_new[[k]][i] * 
        ((Y[i,]-t(Bk_new[[k]])%*%X[i,]) %*% t(Y[i,]-t(Bk_new[[k]])%*%X[i,]))
    }
    sigma_new[[k]] <- (1/sum(w_new[[k]])) * sigma_sum[[k]]
  }
  return(sigma_new)
}

## MSE ##
total_diff_norm <- function(pi_, Bk, sigma, K, l){
  
  pi_new <- pi_[[l+1]]
  Bk_new <- Bk[[l+1]]
  sigma_new <- sigma[[l+1]]
  pi_old <- pi_[[l]]
  Bk_old <- Bk[[l]]
  sigma_old <- sigma[[l]]
  
  pi_norm <- vector()
  Bk_norm <- vector()
  sigma_norm <- vector()
  
  for(k in 1:K){
    pi_norm[k] <- (pi_new[k]-pi_old[k])^2
    Bk_norm[k] <- sum((Bk_new[[k]]-Bk_old[[k]])^2) #same as Frobenius norm
    sigma_norm[k] <- sum((sigma_new[[k]]-sigma_old[[k]])^2) #same as Frobenius norm
  }
  
  return(sqrt(sum(pi_norm,Bk_norm,sigma_norm)))
}


## The modified BIC rounding down to five decimal places ##
modified_BIC <- function(pi_, Bk, density, w, n, K, l){
  
  pi_old <- pi_[[l]]
  w_old <- w[[l]]
  density_old <- density[[l]]
  Bk_old <- Bk[[l]]
  
  logpi_term <- rep(0,K)
  for(k in 1:K){
    logpi_term[k] <- sum(w_old[[k]])*log(pi_old[k])
  }
  logdensity_term <- rep(0,K)
  for(k in 1:K){
    logdensity_term[k] <- sum(w_old[[k]]*lapply(density_old,log)[[k]])
  }
  negative_ll <- (-2)*(sum(logpi_term+logdensity_term))
  d_e <- K + (K-1) + sum(sapply(Bk_old, function(x) colSums(x!=0)))
  
  return(negative_ll + log(n)*d_e)
}

