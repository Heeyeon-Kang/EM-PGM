###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####         R-code for numerical presentation of simulation results.          ####
###################################################################################

# The following R-code includes functions to numerically present the results of the 
# simulations in Section 5.1 - 5.3.


### Packages ###
requiredPackages <- c("MASS", "pracma", "expm", "corrcoverage", 
                      "causact", "fossil", "flexmix", "gtools")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

TPR <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]] - perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]

  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
 
  new_pred_Bk <- pred$Bk
  Number_of_nonzeros <- sapply(true_Bk, function(t) length(which(t != 0)))

  tpr <- vector(length = K)
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    tpr[k] <- length(which(true_Bk[[k]] != 0 & selected)) / Number_of_nonzeros[k]
  }
  
  return(mean(tpr))
}

FPR <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]] - perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  new_pred_Bk <- pred$Bk
  Number_of_zeros <- sapply(true_Bk, function(t) length(which(t == 0)))
  
  fpr <- vector(length = K)
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    fpr[k] <- length(which(true_Bk[[k]] == 0 & selected)) / Number_of_zeros[k]
  }
  
  return(mean(fpr))
}

FDR <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]]-perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  new_pred_Bk <- lapply(1:K, function(x) ifelse(abs(pred$Bk[[x]]) < 5e-6, 0, pred$Bk[[x]]))
  
  fdr <- vector(length = K)
  R <- sapply(new_pred_Bk, function(t) length(which(t != 0)))
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    fdr[k] <- ifelse(R[k] == 0, 0, length(which(true_Bk[[k]] == 0 & selected)) / R[k])
  }
  
  return(mean(fdr))
}

TPR_no_penalty <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]]-perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  new_pred_Bk <- pred$Bk
  Number_of_nonzeros <- sapply(true_Bk, function(t) length(which(t != 0)))
  
  tpr <- vector(length = K)
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    tpr[k] <- length(which(true_Bk[[k]] != 0 & selected)) / Number_of_nonzeros[k]
  }
  
  return(mean(tpr))
}


FPR_no_penalty <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]]-perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  new_pred_Bk <- pred$Bk
  Number_of_zeros <- sapply(true_Bk, function(t) length(which(t == 0)))
  
  fpr <- vector(length = K)
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    fpr[k] <- length(which(true_Bk[[k]] == 0 & selected)) / Number_of_zeros[k]
  }
  
  return(mean(fpr))
}

FDR_no_penalty <- function(true, pred){
  true_Bk <- true$Bk
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]]-perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  new_pred_Bk <- lapply(1:K, function(x) ifelse(abs(pred$Bk[[x]]) < 5e-6, 0, pred$Bk[[x]]))
  
  fdr <- vector(length = K)
  R <- sapply(new_pred_Bk, function(t) length(which(t != 0)))
  for(k in 1:K){
    selected <- (ifelse(abs(new_pred_Bk[[k]]) < 5e-6, 0, new_pred_Bk[[k]]) != 0)
    fdr[k] <- ifelse(R[k] == 0, 0, length(which(true_Bk[[k]] == 0 & selected)) / R[k])
  }
  
  return(mean(fdr))
}

MSE <- function(true, pred){
  true_pi <- true$pi
  true_Bk <- true$Bk
  true_sigma <- true$sigma
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]]-perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  new_pred_pi <- pred$pi
  new_pred_Bk <- pred$Bk
  new_pred_sigma <- pred$sigma
  
  mse <- list(pi = vector(length = K), Bk = vector(length = K), sigma = vector(length = K))
  for(k in 1:K){
    mse$pi[k] <- norm(true_pi[k] - new_pred_pi[k], "2")^2
    mse$Bk[k] <- norm(true_Bk[[k]] - new_pred_Bk[[k]], "F")^2 / length(new_pred_Bk[[k]])
    mse$sigma[k] <- norm(true_sigma[[k]] - new_pred_sigma[[k]], "F")^2 / length(new_pred_sigma[[k]])
  }
  
  return(unlist(lapply(mse, mean)))
}


MSE_no_penalty <- function(true, pred){
  true_pi <- true$pi
  true_Bk <- true$Bk
  true_sigma <- true$sigma
  pred_Bk <- pred$Bk
  
  P <- permutations(K, K)
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(P)){
    perm_Bk <- pred_Bk[P[i,]]
    norm_of_true_and_pred_Bk[i] <- sum(sapply(1:K, function(t) norm(true_Bk[[t]] - perm_Bk[[t]], "F")))
  }
  switch_of_pred <- P[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  new_pred_pi <- pred$pi
  new_pred_Bk <- pred$Bk
  new_pred_sigma <- pred$sigma
  
  mse <- list(pi = vector(length = K), Bk = vector(length = K), sigma = vector(length = K))
  for(k in 1:K){
    mse$pi[k] <- norm(true_pi[k] - new_pred_pi[k], "2")^2
    mse$Bk[k] <- norm(true_Bk[[k]] - new_pred_Bk[[k]], "F")^2 / length(new_pred_Bk[[k]])
    mse$sigma[k] <- norm(true_sigma[[k]] - new_pred_sigma[[k]], "F")^2 / length(new_pred_sigma[[k]])
  }
  
  return(unlist(lapply(mse, mean)))
}

floor_6 <- function(x) as.data.frame(trunc(x*10^5)/10^5)

pred_ll_dens <- function(dat, opt){
  true_B <- dat$true$Bk
  opt_B <- opt$Bk
  opt_sigma <- opt$sigma
  opt_inv_sigma <- lapply(opt$sigma, solve)
  test_x <- dat$test_X
  test_y <- dat$test_Y
  
  dens_list <- lapply(1:length(opt$Bk),
                      function(x) density_f(test_x, test_y,
                                            opt_B[[x]], opt_sigma[[x]], opt_inv_sigma[[x]]))
  return(dens_list)
}

predictive_ll <- function(opt, dens, l, n){
  
  pi_ <- opt[[l]]$pi
  dens <- dens[[l]]
  opt_K <- length(pi_)
  
  in_log <- vector(length = n)
  D <- do.call(cbind, dens)      
  in_log <- as.vector(D %*% pi_)
  
  pred_ll <- (-2)*(1/n)*sum(log(in_log))
  
  return(pred_ll)
}
