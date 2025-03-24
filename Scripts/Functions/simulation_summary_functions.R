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
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_6(total_pred_Bk) != 0))/Number_of_nonzeros
  
  return(True_Positive_Rate)
}

FPR <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)

  False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_6(total_pred_Bk) != 0))/Number_of_zeros
 
  return(False_Positive_Rate)
}

TPR_no_penalty <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  True_Positive_Rate <- length(which(total_true_Bk != 0 & total_pred_Bk != 0))/Number_of_nonzeros
  
  return(True_Positive_Rate)
}

FPR_no_penalty <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  False_Positive_Rate <- length(which(total_true_Bk == 0 & total_pred_Bk != 0))/Number_of_zeros
  
  return(False_Positive_Rate)
}

TPR_flexmix <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]

  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)

  True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_6(total_pred_Bk) != 0))/Number_of_nonzeros
  
  return(True_Positive_Rate)
}

FPR_flexmix <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_6(total_pred_Bk) != 0))/Number_of_zeros
  
  return(False_Positive_Rate)
}

MSE <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K, K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  total_pred_sigma <- sapply(pred$sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  
  MSE_pi <- 1/K*(norm(true$pi-pred$pi, "2"))^2
  MSE_Bk <- 1/K*(norm(total_true_Bk-total_pred_Bk, "F"))^2
  MSE_sigma <- 1/K*(norm(total_true_sigma-total_pred_sigma, "F"))^2
  
  return(c(MSE_pi, MSE_Bk, MSE_sigma))
}

MSE_no_penalty <- function(true, pred){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K, K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk-total_pred_Bk[, permutations(K, K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  total_pred_sigma <- sapply(pred$sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  
  MSE_pi <- 1/K*(norm(true$pi-pred$pi, "2"))^2
  MSE_Bk <- 1/K*(norm(total_true_Bk-total_pred_Bk, "F"))^2
  MSE_sigma <- 1/K*(norm(total_true_sigma-total_pred_sigma, "F"))^2
  
  return(c(MSE_pi, MSE_Bk, MSE_sigma))
}

floor_6 <- function(x) as.data.frame(trunc(x*10^5)/10^5)

predictive_ll <- function(result, l, n){
  optimal <- result$optimal
  density <- result$density
  optimal_K <- result$K

  pi_ <- optimal[[l]]$pi
  dens <- density[[l]]
  opt_K <- optimal_K[[l]]
  
  in_log <- vector(length=n)
  for(i in 1:n){
    in_log[i] <- 0
    for(k in 1:opt_K){
      in_log[i] <- in_log[i]+pi_[k]*dens[[k]][i]
    }
  }
  pred_ll <- (-2)*(1/n)*sum(log(in_log))
  
  return(pred_ll)
}