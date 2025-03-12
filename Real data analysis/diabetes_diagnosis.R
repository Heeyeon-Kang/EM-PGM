##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####            R-code of fitting the diabetes diagnosis data             ####
####                      presented in Section 6.1.                       ####
##############################################################################

# The following R-code fits the diabetes diagnosis data presented in Section 6.1.

# The data is analyzed using mvFMR-MCP.
# But, it can be analyzed by using the other methods.
# (e.g., mvFMR-LASSO, mvFMR-SCAD)

## Packages ##
requiredPackages <- c("ggplot2")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

## Data ##
source("./Data/diabetes_diagnosis_data.R")

## Functions ##
source("./Functions/mvFMR_MCP.R")

## Fitting ##
lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
              30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
              60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
              90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
              35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
              7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 

diabetes_diagnosis_result <- mvFMR_MCP(X, Y, eta=0.5, a=3.7, lambda=lambda_s, maxiter=50)

## Analysis ##
optimal_K <- diabetes_diagnosis_result$K

# correlation #
corr <- vector(length=optimal_K)
for(i in 1:length(corr)){
  corr[i] <- diabetes_diagnosis_result$optimal$sigma[[i]][1,2] / 
    sqrt(diabetes_diagnosis_result$optimal$sigma[[i]][1,1] * 
           diabetes_diagnosis_result$optimal$sigma[[i]][2,2])
}

# coefficients #
sd_X <- c(intercept=1, apply(diabetes_X_1, 2, sd), gender=1, apply(diabetes_X_2, 2, sd))
new_Bk <- list()
for(k in 1:optimal_K){
  new_Bk[[k]] <- apply(diabetes_diagnosis_result$optimal$Bk[[k]], 2, function(x) x/sd_X)
}

# clustering #
cluster <- vector(length=nrow(X))
for(i in 1:length(cluster)){
  cluster[i] <- which.max(lapply(diabetes_diagnosis_result$total_w[[optimal_K]], function(x) x[i]))
}

dat_Y <- as.data.frame(cbind(Y, cluster=cluster))
for(i in 1:length(cluster)){
  if(dat_Y$cluster[i] == 1){
    dat_Y$cluster[i] <- "Group 2"
  }else if(dat_Y$cluster[i] == 2){
    dat_Y$cluster[i] <- "Group 1"
  }
}

# scatter plot #
mvFMR_MCP_plot <- ggplot(data=dat_Y, aes(x=glucose, y=HbA1c, group=cluster)) +
  geom_point(aes(shape=cluster, color=cluster, size=cluster)) + 
  ggthemes::theme_few() + 
  theme(legend.title = element_blank())  +
  labs(x="glucose", y="HbA1c", title="mvFMR-MCP") +
  scale_color_manual(values=c("black", "indianred2")) +
  scale_size_manual(values=c(2, 2)) +
  scale_shape_manual(values=c(15, 16)) +
  theme(axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13)) + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme(legend.text=element_text(size=13))

mvFMR_MCP_plot2 <- mvFMR_MCP_plot + geom_hline(yintercept=6.5, linetype="dashed", size=0.6, col="grey50") + 
  geom_vline(xintercept=126, linetype="dashed", size=0.6, col="grey50")
