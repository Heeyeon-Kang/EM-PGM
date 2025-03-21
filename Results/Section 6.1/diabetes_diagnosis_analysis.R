##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####             R-code for reproducing Table 7 and Figure 1              ####
####       of the diabetes diagnosis data presented in Section 6.1.       ####
##############################################################################

# The following R-code is designed to reproduce Table 7 and Figure 1 
# from Section 6.1, based on the results of analyzing the diabetes diagnosis data 
# using mvFMR-MCP of EM-PGM algorithm.

### Packages ###
requiredPackages <- c("ggplot2")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

### Data ###
source("./Data/diabetes_diagnosis_data.R")

### Result ###
load("./Results/diabetes_diagnosis_result.rda")

### Analysis ###
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
