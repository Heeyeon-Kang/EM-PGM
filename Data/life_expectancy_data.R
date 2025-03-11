##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####  R-code of the data of real data analysis presented in Section 6.2.  ####
##############################################################################

# The raw data can be accessed from 
# https://www.kaggle.com/datasets/lashagoch/life-expectancy-who-updated/data,
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/current-health-expenditure-(che)-as-percentage-of-gross-domestic-product-(gdp)-(-).
# Data is cleaned by removing unnecessary variables and missing observations 
# and using dummy variables.
# Also, the units are different for each variable, 
# so we so we center and scale the predictors.

life_expectancy_dat <- read.csv(file="./Data/life_expectancy_data.csv")
col <- c(8,6,7,23, 5,9,10,11,12,13,14,15,16,20, 21,2,3,4)
life_expectancy_dat2 <- na.omit(life_expectancy_dat[,col])

X <- as.matrix(data.frame(intercept=rep(1, nrow(life_expectancy_dat2)),
                          alcohol=scale(life_expectancy_dat2[,6]),
                          perc_expenditure=scale(life_expectancy_dat2[,5]),
                          Hep_B=scale(life_expectancy_dat2[,7]),
                          Measles=scale(life_expectancy_dat2[,8]),
                          BMI=scale(life_expectancy_dat2[,9]),
                          Polio=scale(life_expectancy_dat2[,10]),
                          Diphtheria=scale(life_expectancy_dat2[,11]),
                          HIV=scale(life_expectancy_dat2[,12]),
                          GDP=scale(life_expectancy_dat2[,13]),
                          schooling=scale(life_expectancy_dat2[,14])))
Y <- as.matrix(data.frame(adult_mortality=life_expectancy_dat2[,1],
                          infant_deaths=life_expectancy_dat2[,2],
                          under_5_deaths=life_expectancy_dat2[,3]))

country <- life_expectancy_dat2[,16]
status <- life_expectancy_dat2[,15]
region <- life_expectancy_dat2[,17]
year <- life_expectancy_dat2[,18]
