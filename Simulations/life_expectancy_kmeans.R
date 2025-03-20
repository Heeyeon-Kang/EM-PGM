##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####             R-code for applying k-means clustering method            ####
####         to the life expectancy data presented in Section 6.2.        ####
##############################################################################

# The following R-code demonstrates how to implement k-means clustering method
# to the life expectancy data presented in Section 6.2.

# Since mvFMR-MCP of the EM-PGM algorithm estimated the optimal K to be 2, 
# we applied k-means clustering method when K = 2.

### Packages ###
requiredPackages <- c("ggplot2", "plotly")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

### Data ###
source("./Data/life_expectancy_data.R")

### Analysis ###
# k-means clustering method with K=2 #
set.seed(2868)
group <- kmeans(Y, centers=2)$cluster

dat_Y <- as.data.frame(cbind(Y, group=group))
for(i in 1:nrow(Y)){
  if(dat_Y$group[i] == 1){
    dat_Y$group[i] <- "Group 2"
  }else if(dat_Y$group[i] == 2){
    dat_Y$group[i] <- "Group 1"
  }
}

# scatter plot #
k_means_plot <- plot_ly(dat_Y, x=~adult_mortality, y=~infant_deaths, z=~under_5_deaths, 
                    color=~group, symbol=~group, symbols=c("circle", "square"),
                    colors=c("black", "indianred2"), type="scatter3d", mode="markers",
                    marker=list(size=2, sizemode="diameter", sizeref=0.05))

k_means_plot2 <- k_means_plot %>% 
  layout(title=list(x=0.55, y=0.8, text="k-means clustering with K=2", font=list(size=18)),
  scene=list(aspectmode="manual", aspectratio=list(x=2, y=0.6, z=0.8), 
             xaxis=list(title="Adult Mortality", titlefont=list(size=14), tickfont=list(size=10)),
             yaxis=list(title="Infant Deaths", titlefont=list(size=14), tickfont=list(size=10)),
             zaxis=list(title="Children Under 5 Deaths", titlefont=list(size=14), tickfont=list(size=10)),
             camera=list(eye=list(x=1, y=-2, z=1))),
  legend=list(x=0.9, y=0.5, xanchor="left", yanchor="middle", font=list(size=14))) %>% 
  layout(margin=list(l=0, r=0, t=0, b=0))
