# EM-PGM

### Electronic supplement to "Penalized estimation for a finite mixture of regression models"


This repository has the R-code to do the analyses and reproduce the figures and tables in Sections 5 and Section 6 of the manuscript and the supplementary material.

All code was written in a laptop with R version 4.4.2 and run on a Linux Cluster with R version 4.1.2.
To run the R-code, it is recommended to load the R project ***EM-PGM.Rproj***, as all paths are set relative to this directory.

The repository consists of the following folders:

* Data: R-code for generating or refining the data used in Section 5 and Section 6;
  * ***simulation_seed_number.R*** contains the seed numbers generating the data of simulations in Section 5.
  * ***simulation_data.R*** contains the functions generating the dataset using in Section 5.
  * ***diabetes_diagnosis.csv*** can be accessed from **Diabetes data** section in <http://hbiostat.org/data/>.
  * ***life_expectancy.csv*** can be accessed from <https://www.kaggle.com/datasets/lashagoch/life-expectancy-who-updated/data>, while additional relevant data is available from <https://www.who.int/data/gho/data/indicators/indicator-details/GHO/current-health-expenditure-(che)-as-percentage-of-gross-domestic-product-(gdp)-(-)>.
  * ***diabetes_diagnosis_data.R*** and ***life_expectancy_data.R*** are the R-codes of the process of refining the raw data, ***diabetes_diagnosis.csv*** and ***life_expectancy.csv***, respectively.
 
* Functions: R-code of all functions for running the EM-PGM algorithm and R-code for fitting simulation data using each method;
  * ***functions.R*** is the R-code of all functions for running EM-PGM algorithm.
  * ***mvFMR.R***, ***mvFMR_LASSO.R***, ***mvFMR_SCAD.R***, and ***mvFMR_MCP.R*** contain the R-codes for fitting simulation data using mvFMR with their respective penalty functions presented in Section 5.3.
  * ***mvFMR_fixedK.R***, ***mvFMR_LASSO_fixedK.R***, ***mvFMR_SCAD_fixedK.R***, and ***mvFMR_MCP_fixedK.R*** contain the R-codes for fitting the data using mvFMR with their respective penalty functions presented in Section 5.1, Section 5.2, and Section 5.3, assuming K is known.
  * ***flexmix_fixedK.R*** contains the R-code for fitting the data using the R-package "Flexmix" presented in Section 5.1 and Section 5.2.
  * ***oracle_fixedK.R*** contains the R-code to estimate the oracle estimator presented in Section 5.1 and Section 5.2 by fitting the data.

* Simulations: The results of simulation studies and real data analyses to reproduce the figures and tables presented in Section 5 and Section 6;
  * ***simulation_table.R*** contains the functions for calculating TPR, FPR, MSE, and predictive log-likelihood loss.
  * ***diabetes_diagnosis_table.R*** and ***life_expectancy_table.R*** contain the R-codes for analyzing the results of the data and generating Table 7 and Table 8 presented in Section 6.1 and Section 6.2.
  * ***diabetes_diagnosis_kmeans.R*** and ***life_expectancy_kmeans.R*** contains the R-codes for applying k-means clustering method in diabetes dianosis data and life expectancy data.
  * ***diabetes_diagnosis.R*** and ***life_expectancy.R*** contain the R-codes for fitting the data using mvFMR-MCP presented in Section 6.1 and Section 6.2.
  * ***s1_table.rda***, ***s2_table.rda***, and ***s3_table.rda*** are the rda files that contain the contents of Table 2, Table 4, and Table 6 in Section 5.1, Section 5.2, and Section 5.3, respectively.
  * ***diabetes_diagnosis.rda*** and ***life_expectancy.rda*** are rda files containing the results of the analysis using mvFMR-MCP.

    
