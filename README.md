# EM-PGM

### Electronic supplement to "Penalized estimation for a finite mixture of regression models"


This repository has the R-code to do the analyses and reproduce the figures and tables in Sections 5 and Section 6 of the manuscript and the supplementary material.

All code was written in a laptop with R version 4.4.2 and run on a Linux Cluster with R version 4.1.2.
To run the R-code, it is recommended to load the R project ***EM-PGM.Rproj***, as all paths are set relative to this directory.

The repository consists of the following folders:

* Data: R-codes for generating or refining the data used in Section 5 and Section 6;
  * ***simulation_seed_number.R*** contains the seed numbers generating the data of simulations in Section 5.
  * ***simulation_data.R*** contains the functions generating the dataset using in Section 5.
  * ***diabetes_diagnosis.csv*** can be accessed from **Diabetes data** section in <http://hbiostat.org/data/>.
  * ***life_expectancy.csv*** can be accessed from <https://www.kaggle.com/datasets/lashagoch/life-expectancy-who-updated/data>, while additional relevant data is available from <https://www.who.int/data/gho/data/indicators/indicator-details/GHO/current-health-expenditure-(che)-as-percentage-of-gross-domestic-product-(gdp)-(-)>.
  * ***diabetes_diagnosis_data.R*** and ***life_expectancy_data.R*** are the R-codes of the process of refining the raw data, ***diabetes_diagnosis_data.csv*** and ***life_expectancy_data.csv***, respectively.
 
* Scripts;
  * Analysis: R-codes for implementing 'flexmix', 'oracle', 'mvFMR', 'mvFMR-LASSO', 'mvFMR-SCAD', and 'mvFMR-MCP' to simulation data;
    * ***flexmix_fixedK.R*** contains R-code to fit the R-package "flexmix" to the simulated data presented in Section 5.1 - 5.2.
    * ***oracle_fixedK.R*** contains R-code to estimate the oracle estimator presented in Section 5.1 - 5.2.
    * ***mvFMR_fixedK.R*** contains R-code for implementing 'mvFMR' to simulated data presented in Section 5.1 - 5.2.
    * ***PGM_mvFMR_LASSO_fixedK.R***, ***PGM_mvFMR_SCAD_fixedK.R***, ***PGM_mvFMR_MCP_fixedK.R***, ***ADMM_mvFMR_LASSO_fixedK.R***, ***ADMM_mvFMR_SCAD_fixedK.R***, and ***ADMM_mvFMR_MCP_fixedK.R*** are R-codes for implementing 'mvFMR' with their respective penalty functions(LASSO, SCAD, and MCP) and respective algorithms(the EM-PGM algorithm and the ADMM solver) assuming that the number of components K is known, which corresponds to the setting in Section 5.1 - 5.2.
    * ***mvFMR.R*** contains R-code for implementing 'mvFMR' to simulated data presented in Section 5.3.
    * ***PGM_mvFMR_LASSO.R***, ***PGM_mvFMR_SCAD.R***, ***PGM_mvFMR_MCP.R***, ***ADMM_mvFMR_LASSO.R***, ***ADMM_mvFMR_SCAD.R***, and ***ADMM_mvFMR_MCP.R*** are R-codes for implementing 'mvFMR' with their respective penalty functions(LASSO, SCAD, and MCP) and respective algorithms(the EM-PGM algorithm and the ADMM solver) assuming that the number of components K is unknown, which corresponds to the setting in Section 5.3.
  * Functions: R-codes of functions to implement all methods;
    * ***functions.R*** contains R-codes composed of functions to perform all methods in 'Analysis'.
    * ***simulation_summary_functions.R*** contains R-codes composed of functions for numerical presentation of simulation results.
 
* Simulations: The results of simulation studies and real data analyses to reproduce the figures and tables presented in Section 5 and Section 6;
  * ***simulation_table.R*** contains the functions for calculating TPR, FPR, MSE, and predictive log-likelihood loss.
  * ***s1_table.rda***, ***s2_table.rda***, and ***s3_table.rda*** are the rda files that contain the contents of Table 2, Table 4, and Table 6 in Section 5.1, Section 5.2, and Section 5.3, respectively.
  * ***diabetes_diagnosis.R*** and ***life_expectancy.R*** contain the R-codes for fitting and analyzing the data using mvFMR-MCP and generating Table 7 and Table 8 presented in Section 6.1 and Section 6.2.
  * ***diabetes_diagnosis_result.rda*** and ***life_expectancy_result.rda*** are rda files containing the results after applying mvFMR-MCP.
  * ***diabetes_diagnosis_kmeans.R*** and ***life_expectancy_kmeans.R*** contains the R-codes for applying k-means clustering method in diabetes dianosis data and life expectancy data.


    
