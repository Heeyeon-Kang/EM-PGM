###################################################################################
#### Penalized estimation in finite mixtures of multivariate regression models ####
####                         via the EM-PGM algorithm                          ####
###################################################################################
####                           Author: Heeyeon Kang                            ####
####                         Supervisor: Sunyoung Shin                         ####
###################################################################################
####  R-code of the seed numbers generating the data of simulations presented  ####
####                               in Section 5.                               ####
###################################################################################

# The following R-code contains the seed numbers generating the data of simulations 
# presented in Section 5.1, Section 5.2, and Section 5.3.

# We set the seed numbers to fix all the dataset used in our simulations.
# The seed numbers of Section 5.1, Section 5.2, and Section 5.3 are 
# 2,5,10 / 210,220,230 / 302,304, respectively, for each model.

# Since we generated 500 independent simulations for each model, the following
# R-code contains 500 seed numbers for each model.


### Section 5.1 : model 1 ###
set.seed(511)
seed_number_5.1.1 <- sample(1:1e5, 500, replace = FALSE)

### Section 5.1 : model 2 ###
set.seed(512)
seed_number_5.1.2 <- sample(1:1e5, 500, replace = FALSE)

### Section 5.1 : model 3 ###
set.seed(513)
seed_number_5.1.3 <- sample(1:1e5, 500, replace = FALSE)


### Section 5.2 : model 1 ###
set.seed(521)
seed_number_5.2.1 <- sample(1:1e5, 500, replace = FALSE)

### Section 5.2 : model 2 ###
set.seed(522)
seed_number_5.2.2 <- sample(1:1e5, 500, replace = FALSE)

### Section 5.2 : model 3 ###
set.seed(523)
seed_number_5.2.3 <- sample(1:1e5, 500, replace = FALSE)


### Section 5.3 : model 1 ###
set.seed(531)
seed_number_5.3.1 <- sample(1:1e5, 500, replace = FALSE)

### Section 5.3 : model 2 ###
set.seed(532)
seed_number_5.3.2 <- sample(1:1e5, 500, replace = FALSE)
