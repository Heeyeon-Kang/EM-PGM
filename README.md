# EM-PGM

### Electronic supplement to "Penalized estimation in finite mixtures of multivariate regression models via the EM-PGM algorithm"

This repository contains the R and Python implementations used to reproduce the analyses, tables, and results presented in Sections 5 and 6 of the manuscript and in the Supplementary Material.

The code was developed using R version 4.4.2 and Python version 3.13.9 on a laptop.

The R code was tested on a Linux cluster with R version 4.1.2. The R code should be run from the R project ***EM-PGM.Rproj***, as all file paths are defined relative to the project directory.

The Python code should be run from the project root directory, as all paths are defined relative to this directory. To keep the repository lightweight, the precomputed simulation results for Sections 5.1–5.3 are not included. These results can be reproduced by running the corresponding simulation scripts. The scripts and results for the real data analyses in Sections 6.1 and 6.2 are included.

The repository consists of the following folders:

-   R: R code for conducting the analyses.
    -   Data: R scripts for generating or preprocessing the data used in Sections 5 and 6.
        -   ***simulation_seed_number.R*** contains the seed numbers used to generate the simulation results in Section 5.
        -   ***simulation_data.R*** contains functions for generating the datasets used in Section 5.
        -   ***diabetes_diagnosis_data.csv*** can be accessed from **Diabetes data** section in <http://hbiostat.org/data/>.
        -   ***CCLE_genes.csv*** and ***CCLE_drugs.csv*** represent the predictor and response data used in the CCLE data analysis. The complete CCLE dataset is accessible at <https://sites.broadinstitute.org/ccle/datasets>. Because the original files are very large and this repository is intended to contain only the data necessary to reproduce the analysis.
        -   ***diabetes_diagnosis_data.R*** and ***CCLE_data.R*** contain the R scripts for preprocessing ***diabetes_diagnosis_data.csv*** and ***CCLE_genes.csv, CCLE_drugs.csv***, respectively.
    -   Scripts
        -   Analysis: R scripts for applying existing and proposed methods to the simulated data. ***fixedK*** indicates that the number of components K is fixed.
            -   ***flexmix_fixedK.R*** contains the R script used to fit the `flexmix` model to the simulated data presented in Sections 5.1 and 5.2.
            -   ***oracle_fixedK.R*** contains the R script used to implement the oracle method for the simulated data presented in Sections 5.1 and 5.2.
            -   ***mvFMR_fixedK.R*** contains the implementation of the proposed method(mvFMR) for the simulated data presented in Sections 5.1 and 5.2.
            -   ***PGM_mvFMR_LASSO_fixedK.R***, ***PGM_mvFMR_SCAD_fixedK.R***, and ***PGM_mvFMR_MCP_fixedK.R*** contain implementations of the proposed method under different penalty functions(LASSO, SCAD, and MCP) using the EM-PGM algorithm for the simulated data presented in Sections 5.1 and 5.2.
            -   ***ADMM_mvFMR_LASSO_fixedK.R***, ***ADMM_mvFMR_SCAD_fixedK.R***, and ***ADMM_mvFMR_MCP_fixedK.R*** contain implementations of the proposed method under different penalty functions using the EM-ADMM algorithm for the simulated data presented in Sections 5.1 and 5.2.
            -   ***mvFMR.R*** contains the implementation of the proposed method for the simulated data presented in Section 5.3.
            -   ***PGM_mvFMR_LASSO.R***, ***PGM_mvFMR_SCAD.R***, and ***PGM_mvFMR_MCP.R*** contain implementations of the proposed method under different penalty functions using the EM-PGM algorithm for the simulated data presented in Section 5.3.
            -   ***ADMM_mvFMR_LASSO.R***, ***ADMM_mvFMR_SCAD.R***, and ***ADMM_mvFMR_MCP.R*** contain implementations of the proposed method under different penalty functions using the EM-ADMM algorithm for the simulated data presented in Section 5.3.
        -   Functions: R functions used by the analysis scripts.
            -   ***functions.R*** contains R functions used to implement all methods in the 'Analysis' scripts.
            -   ***simulation_summary_functions.R*** contains R functions used to generate numerical summaries of the simulation results.
    -   Results
        -   Section5.1: Contains the summarized simulation results and R scripts for reproducing Tables 2--4.
            -   ***Model1***, ***Model2***, and ***Model3*** contain the summarized results for each simulation model.
            -   ***sim1_analysis.R*** reproduces Tables 2--4 presented in Section 5.1.
            -   ***mod1_table.rda***, ***mod2_table.rda***, and ***mod3_table.rda*** contain the summarized results used to generate the corresponding tables.
        -   Section5.2: Contains the summarized simulation results and R scripts for reproducing Tables 6--8.
            -   ***Model4***, ***Model5***, and ***Model6*** contain the summarized results for each simulation model.
            -   ***sim2_analysis.R*** reproduces Tables 6--8 presented in Section 5.2.
            -   ***mod4_table.rda***, ***mod5_table.rda***, and ***mod6_table.rda*** contain the summarized results used to generate the corresponding tables.
        -   Section5.3: Contains the summarized simulation results and R scripts for reproducing Tables 10.
            -   ***Model7*** and ***Model8*** contain the summarized results for each simulation model.
            -   ***sim3_analysis.R*** reproduces Tables 10 presented in Section 5.3.
            -   ***mod7_table.rda, mod7_K_table.rda***, ***mod8_table.rda***, and ***mod8_K_table.rda*** contain the summarized results used to generate the corresponding tables.
        -   Section6.1: Contains the results from fitting the models and the R script for reproducing Table 11.
            -   ***PGM_diabetes_diagnosis.rda*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-PGM algorithm.
            -   ***ADMM_diabetes_diagnosis.rda*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-ADMM algorithm.
            -   ***diabetes_diagnosis_analysis.R*** contains the R script for reproducing Table 11 presented in Section 6.1.
        -   Section6.2: Contains the results from fitting the models and the R script for reproducing Table 12.
            -   ***PGM_CCLE.rda*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-PGM algorithm.
            -   ***ADMM_CCLE.rda*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-ADMM algorithm.
            -   ***CCLE_analysis.R*** contains the R script for reproducing Table 12 presented in Section 6.2.
    -   Simulations: R scripts for conducting simulation studies and real data analyses, and for reproducing the figures and tables presented in Sections 5 and 6.
        -   ***diabetes_diagnosis.R*** and ***CCLE.R*** contain R scripts for applying 'mvFMR-MCP' using the EM-PGM and the EM-ADMM algorithm to the datasets presented in Sections 6.1 and 6.2, respectively.
        -   ***simulation1.R***, ***simulation2.R***, and ***simulation3.R*** contain R scripts for implementing all methods on simulated data presented in Sections 5.1--5.3 using parallel computing.
-   Python: Python code for conducting the analyses.
    -   Data: Python scripts for generating or preprocessing the data used in Sections 5 and 6.
        -   ***simulation_seed_number.py*** contains the seed numbers used to generate the simulation results in Section 5.
        -   ***simulation_data.py*** contains functions for generating the datasets used in Section 5.
        -   ***diabetes_diagnosis_data.csv*** can be accessed from **Diabetes data** section in <http://hbiostat.org/data/>.
        -   ***CCLE_genes.csv*** and ***CCLE_drugs.csv*** represent the predictor and response data used in the CCLE data analysis. The complete CCLE dataset is accessible at <https://sites.broadinstitute.org/ccle/datasets>. Because the original files are very large and this repository is intended to contain only the data necessary to reproduce the analysis.
        -   ***diabetes_diagnosis_data.py*** and ***CCLE_data.py*** contain the Python scripts for preprocessing ***diabetes_diagnosis_data.csv*** and ***CCLE_genes.csv, CCLE_drugs.csv***, respectively.
    -   Scripts
        -   Analysis: Python scripts for applying existing and proposed methods to the simulated data. ***fixedK*** indicates that the number of components K is fixed.
            -   ***oracle_fixedK.py*** contains the Python script used to implement the oracle method for the simulated data presented in Sections 5.1 and 5.2.
            -   ***mvFMR_fixedK.py*** contains the implementation of the proposed method(mvFMR) for the simulated data presented in Sections 5.1 and 5.2.
            -   ***PGM_mvFMR_LASSO_fixedK.py***, ***PGM_mvFMR_SCAD_fixedK.py***, and ***PGM_mvFMR_MCP_fixedK.py*** contain implementations of the proposed method under different penalty functions(LASSO, SCAD, and MCP) using the EM-PGM algorithm for the simulated data presented in Sections 5.1 and 5.2.
            -   ***ADMM_mvFMR_LASSO_fixedK.py***, ***ADMM_mvFMR_SCAD_fixedK.py***, and ***ADMM_mvFMR_MCP_fixedK.py*** contain implementations of the proposed method under different penalty functions using the EM-ADMM algorithm for the simulated data presented in Sections 5.1 and 5.2.
            -   ***mvFMR.py*** contains the implementation of the proposed method for the simulated data presented in Section 5.3.
            -   ***PGM_mvFMR_LASSO.py***, ***PGM_mvFMR_SCAD.py***, and ***PGM_mvFMR_MCP.py*** contain implementations of the proposed method under different penalty functions using the EM-PGM algorithm for the simulated data presented in Section 5.3.
            -   ***ADMM_mvFMR_LASSO.py***, ***ADMM_mvFMR_SCAD.py***, and ***ADMM_mvFMR_MCP.py*** contain implementations of the proposed method under different penalty functions using the EM-ADMM algorithm for the simulated data presented in Section 5.3.
        -   Functions: Python functions used by the analysis scripts.
            -   ***functions.py*** contains Python functions used to implement all methods in the 'Analysis' scripts.
            -   ***simulation_summary_functions.py*** contains Python functions used to generate numerical summaries of the simulation results.
    -   Results
        -   Section5.1
            -   ***sim1_analysis.py*** demonstrates how to compute TPR, FPR, FDR, and MSE for a single run of the simulation study presented in Section 5.1.
        -   Section5.2
            -   ***sim2_analysis.py*** demonstrates how to compute TPR, FPR, FDR, and MSE for a single run of the simulation study presented in Section 5.2.
        -   Section5.3
            -   ***sim3_analysis.py*** demonstrates how to compute predictive log-likelihood for a single run of the simulation study presented in Section 5.3.
        -   Section6.1: Contains the results from fitting the models and the Python script for reproducing Table 11.
            -   ***PGM_diabetes_diagnosis.pkl*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-PGM algorithm.
            -   ***ADMM_diabetes_diagnosis.pkl*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-ADMM algorithm.
            -   ***diabetes_diagnosis_analysis.py*** contains the Python script for reproducing Table 11 presented in Section 6.1.
        -   Section6.2: Contains the results from fitting the models and the Python script for reproducing Table 12.
            -   ***PGM_CCLE.pkl*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-PGM algorithm.
            -   ***ADMM_CCLE.pkl*** contains the fitting results obtained using 'mvFMR-MCP' with the EM-ADMM algorithm.
            -   ***CCLE_analysis.py*** contains the Python script for reproducing Table 12 presented in Section 6.2.
    -   Simulations: Python scripts for conducting simulation studies and real data analyses, and for reproducing the figures and tables presented in Sections 5 and 6.
        -   ***diabetes_diagnosis.py*** and ***CCLE.py*** contain Python scripts for applying 'mvFMR-MCP' using the EM-PGM and the EM-ADMM algorithm to the datasets presented in Sections 6.1 and 6.2, respectively.
