# EM-PGM

### Electronic supplement to Penalized estimation in finite mixtures of multivariate regression models via the EM-PGM algorithm"


This repository has the R-code to do the analyses and reproduce the figures and tables in Sections 5 and Section 6 of the manuscript and the supplementary material.

All code was written in a laptop with R version 4.4.2 and run on a Linux Cluster with R version 4.1.2.
To run the R-code, it is recommended to load the R project ***EM-PGM.Rproj***, as all paths are set relative to this directory.

The repository consists of the following folders:

* Data: R-codes for generating or refining the data used in Section 5 and Section 6;
  * ***simulation_seed_number.R*** contains the seed numbers generating the data of simulations in Section 5.
  * ***simulation_data.R*** contains the functions generating the dataset using in Section 5.
  * ***diabetes_diagnosis_data.csv*** can be accessed from **Diabetes data** section in <http://hbiostat.org/data/>.
  * ***age_specific_mortality_rates_data.csv*** can be accessed from <https://www.kaggle.com/datasets/lashagoch/life-expectancy-who-updated/data>, while additional relevant data is available from <https://www.who.int/data/gho/data/indicators/indicator-details/GHO/current-health-expenditure-(che)-as-percentage-of-gross-domestic-product-(gdp)-(-)>.
  * ***diabetes_diagnosis_data.R*** and ***age_specific_mortality_rates_data.R*** are the R-codes of the process of refining the raw data, ***diabetes_diagnosis_data.csv*** and ***age_specific_mortality_rates_data.csv***, respectively.
 
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

* Results;
  * Section 5.1: The '.rda' files containing the summarized results collected in Section 5.1 and R-code for reproducing Table 2;
    * ***Model 1***, ***Model 2***, and ***Model 3*** are folders that contain the summarized results of Model 1, Model 2, and Model 3, respectively.
    * ***sim1_analysis.R*** contains R-code for reproducing Table 2 in Section 5.1 and ***Table 2.rda*** is the corresponding file.
  * Section 5.2: The '.rda' files containing the summarized results collected in Section 5.2 and R-code for reproducing Table 4;
    * ***Model 4***, ***Model 5***, and ***Model 6*** are folders that contain the summarized results of Model 4, Model 5, and Model 6, respectively.
    * ***sim2_analysis.R*** contains R-code for reproducing Table 4 in Section 5.2 and ***Table 4.rda*** is the corresponding file.
  * Section 5.3: The '.rda' files containing the summarized results collected in Section 5.3 and R-code for reproducing Table 5;
    * ***Model 7*** and ***Model 8*** are folders that contain the summarized results of Model 7 and Model 8, respectively.
    * ***sim3_analysis.R*** contains R-code for reproducing Table 6 in Section 5.3 and ***Table 6.rda*** is the corresponding file.
  * Section 6.1;
    * ***diabetes_diagnosis_result.rda*** contains all results after fitting the data on 'mvFMR-MCP' with the EM-PGM algorithm.
    * ***diabetes_diagnosis_analysis.R*** contains R-code for reproducing Table 7 and Figure 1 of the diabetes diagnosis data presented in Section 6.1.
  * Section 6.2;
    * ***age_specific_mortality_rates_result.rda*** contains all results after fitting the data on 'mvFMR-MCP' with the EM-PGM algorithm.
    * ***age_specific_mortality_rates_analysis.R*** contains R-code for reproducing Table 8 and Figure 2 of the age-specific mortality rates data presented in Section 6.2.

* Simulations: The results of simulation studies and real data analyses to reproduce the figures and tables presented in Section 5 and Section 6;
  * ***diabetes_diagnosis.R*** and ***age_specific_mortality_rates.R*** contain R-codes for applying 'mvFMR-MCP' and 'the EM-PGM algorithm' to the diabetes diagnosis data and age-specific mortality rates data presented in Section 6.1 and 6.2, respectively.
  * ***diabetes_diagnosis_kmeans.R*** and ***age_specific_mortality_rates_kmeans.R*** contain R-codes for applying 'k-means clustering method' to the diabetes diagnosis data and age-specific mortality rates data presented in Section 6.1 and 6.2, respectively.
  * ***simulation1.R***, ***simulation2.R***, and ***simulation3.R*** contains R-codes for implementing all methods on simulated data presented in Section 5.1 - 5.3, using parallel computing.
    
