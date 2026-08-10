###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####    R-code of the data of real data analysis presented in Section 6.2.     ####
###################################################################################

# The raw data can be accessed from https://sites.broadinstitute.org/ccle/datasets.
# We focus on three chemical compounds, Erlotinib, AZD6244, and PD-0325901.
# Data is cleaned by removing unnecessary variables and missing observations.
# We use only the lines tested for all three, yielding n=489 cell lines. 
# For genes, we keep P=300 genes that are highly correlated with three responses 
# by calculating the sum of its absolute correlations with the three responses.
# We center and scale the predictors.

X <- read.table("./Data/CCLE_genes.csv", sep = ",", header = TRUE)
Y <- read.table("./Data/CCLE_drugs.csv", sep = ",", header = TRUE)

C <- cor(X, Y, use = "pairwise.complete.obs")
score <- rowSums(abs(C))

ord <- order(score, decreasing = TRUE)
Xr <- X[, ord[1:300], drop = FALSE]
top_genes <- colnames(X)[ord[1:300]]

X <- cbind(intercept = rep(1, nrow(Xr)), scale(Xr))
X_sd <- apply(Xr, 2, sd)
