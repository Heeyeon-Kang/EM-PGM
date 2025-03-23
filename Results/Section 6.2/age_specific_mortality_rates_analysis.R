##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####             R-code for reproducing Table 8 and Figure 2              ####
####  of the age-specific mortality rates data presented in Section 6.2.  ####
##############################################################################

# The following R-code is designed to reproduce Table 8 and Figure 2 from Section 6.2, 
# based on the results of analyzing the age-specific mortality rates data 
# applied mvFMR-MCP with EM-PGM algorithm.

### Packages ###
requiredPackages <- c("ggplot2", "plotly")
for(p in requiredPackages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

### Data ###
source("./Data/age_specific_mortality_rates_data.R")

### Result ###
load("./Results/Section 6.2/age_specific_mortality_rates_result.rda")

### Analysis ###
optimal_K <- age_specific_mortality_rates_result$K

# correlation #
corr <- list()
for(k in 1:optimal_K){
  corr[[k]] <- as.data.frame(diag(ncol(Y)))
  colnames(corr[[k]]) <- c("adult mortality", "infant deaths", "under 5 deaths")
  rownames(corr[[k]]) <- c("adult mortality", "infant deaths", "under 5 deaths")
  
  corr[[k]][1,2] <- age_specific_mortality_rates_result$optimal$sigma[[k]][1,2] / 
    sqrt(age_specific_mortality_rates_result$optimal$sigma[[k]][1,1]*age_specific_mortality_rates_result$optimal$sigma[[k]][2,2])
  corr[[k]][2,1] <- corr[[k]][1,2]
  
  corr[[k]][1,3] <- age_specific_mortality_rates_result$optimal$sigma[[k]][1,3] / 
    sqrt(age_specific_mortality_rates_result$optimal$sigma[[k]][1,1]*age_specific_mortality_rates_result$optimal$sigma[[k]][3,3])
  corr[[k]][3,1] <- corr[[k]][1,3]
  
  corr[[k]][2,3] <- age_specific_mortality_rates_result$optimal$sigma[[k]][2,3] / 
    sqrt(age_specific_mortality_rates_result$optimal$sigma[[k]][2,2]*age_specific_mortality_rates_result$optimal$sigma[[k]][3,3])
  corr[[k]][3,2] <- corr[[k]][2,3]
}

# coefficients #
sd_X <- c(intercept=1, apply(age_specific_mortality_rates_dat2[,c(6,5,7,8,9,10,11,12,13,14)], 2, sd))
new_Bk <- list()
for(k in 1:optimal_K){
  new_Bk[[k]] <- apply(age_specific_mortality_rates_result$optimal$Bk[[k]],2,function(x) x/sd_X)
}

# clustering #
cluster <- vector(length=nrow(X))
for(i in 1:nrow(X)){
  cluster[i] <- which.max(lapply(age_specific_mortality_rates_result$optimal_w, function(x) x[i]))
}


dat_Y <- as.data.frame(cbind(Y, cluster=cluster))
for(i in 1:nrow(X)){
  if(dat_Y$cluster[i] == 1){
    dat_Y$cluster[i] <- "Group 1"
  }else if(dat_Y$cluster[i] == 2){
    dat_Y$cluster[i] <- "Group 2"
  }
}

total_dat <- cbind(dat_Y, country=age_specific_mortality_rates_dat2[,16], 
                   year=age_specific_mortality_rates_dat2[,18], 
                   status=age_specific_mortality_rates_dat2[,15],
                   region=age_specific_mortality_rates_dat2[,17])

# scatter plot #
mvFMR_MCP_plot <- plot_ly(total_dat, x=~adult_mortality, y=~infant_deaths, z=~under_5_deaths, 
                          color=~cluster, symbol=~cluster, symbols=c("circle", "square"),
                          colors=c("black", "indianred2"), type="scatter3d", mode="markers",
                          marker=list(size=2, sizemode="diameter", sizeref=0.05))

mvFMR_MCP_plot2 <- mvFMR_MCP_plot %>% 
  layout(title=list(x=0.55, y=0.8, text="mvFMR-MCP", font=list(size=18)),
         scene=list(aspectmode="manual", aspectratio=list(x=2, y=0.6, z=0.8),
                    xaxis=list(title="Adult Mortality", titlefont=list(size=14), tickfont=list(size=10)),
                    yaxis=list(title="Infant Deaths", titlefont=list(size=14), tickfont=list(size=10)),
                    zaxis=list(title="Children Under 5 Deaths", titlefont=list(size=14), tickfont=list(size=10)),
                    camera=list(eye=list(x=1, y=-2, z=1))),
         legend=list(x=0.9, y=0.5, xanchor="left", yanchor="middle", font=list(size=14))) %>% 
  layout(margin=list(l=0, r=0, t=0, b=0))
