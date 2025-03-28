###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####               R-code for applying k-means clustering method               ####
####         to the diabetes diagnosis data presented in Section 6.1.          ####
###################################################################################

# The following R-code demonstrates how to implement k-means clustering method
# to the diabetes diagnosis data presented in Section 6.1.

# Since mvFMR-MCP of the EM-PGM algorithm estimated the optimal K to be 2, 
# we applied k-means clustering method when K = 2.

### Packages ###
requiredPackages <- c("ggplot2")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

### Data ###
source("./Data/diabetes_diagnosis_data.R")

### Analysis ###
# k-means clustering method with K=2 #
set.seed(2868)
group <- kmeans(apply(Y, 2, scale), centers=2)$cluster

dat_Y <- as.data.frame(cbind(Y, group=group))
for(i in 1:nrow(Y)){
  if(dat_Y$group[i] == 1){
    dat_Y$group[i] <- "Group 1"
  }else if(dat_Y$group[i] == 2){
    dat_Y$group[i] <- "Group 2"
  }
}

# scatter plot #
k_means_plot <- ggplot(data=dat_Y, aes(x=glucose, y=HbA1c, group=group)) +
  geom_point(aes(shape=group, color=group, size=group)) + 
  ggthemes::theme_few() +
  theme(legend.title=element_blank())  +
  labs(x="glucose", y="HbA1c", title="k-means clustering with K=2") +
  scale_color_manual(values=c("black", "indianred2")) +
  scale_size_manual(values=c(2, 2)) +
  scale_shape_manual(values=c(15, 16)) +
  theme(axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13)) + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme(legend.text=element_text(size=13))

k_means_plot2 <- k_means_plot + geom_hline(yintercept=6.5, linetype="dashed", size=0.6, col="grey50") + 
  geom_vline(xintercept=126, linetype="dashed", size=0.6, col="grey50")
