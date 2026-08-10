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

Data <- read.table("./Data/CCLE_expression.csv", sep = ",")
response <- read.table("./Data/CCLE_Drugdata.csv", header = T, sep = ",")

## Cell line names ##
w <- as.character(Data[1, -c(1,2)])

## Gene names ##
Gene <- as.character(Data[2:18927, 2])
Xgenes <- as.matrix(Data[2:18927, 3:1039])    
storage.mode(Xgenes) <- "numeric"              
colnames(Xgenes) <- w
rownames(Xgenes) <- Gene

z1 <- response[response[,3] == "Erlotinib", 1]
z2 <- response[response[,3] == "AZD6244", 1]
z3 <- response[response[,3] == "PD-0325901", 1]

area1 <- response[response[,3] == "Erlotinib", 13]
area2 <- response[response[,3] == "AZD6244", 13]
area3 <- response[response[,3] == "PD-0325901", 13]

## Common cell line for three drugs ##
z <- Reduce(intersect, list(z1, z2, z3))

i1 <- match(z, z1)
j2 <- match(z, z2)
k3 <- match(z, z3)
area <- cbind(Erlotinib = area1[i1], AZD6244 = area2[j2], PD_0325901 = area3[k3])
rownames(area) <- z

w_key <- toupper(trimws(colnames(Xgenes)))
z_key <- toupper(trimws(rownames(area)))

J_all <- match(z_key, w_key)  
I <- which(!is.na(J_all))
J <- J_all[I]        

Y <- area[I, , drop=FALSE]
X <- t(Xgenes[, J, drop=FALSE])

C <- cor(X, Y, use = "pairwise.complete.obs")
score <- rowSums(abs(C))

ord <- order(score, decreasing = TRUE)
Xr <- X[, ord[1:300], drop = FALSE]
top_genes <- colnames(X)[ord[1:300]]

X <- cbind(intercept = rep(1, nrow(Xr)), scale(Xr))
X_sd <- apply(Xr, 2, sd)
