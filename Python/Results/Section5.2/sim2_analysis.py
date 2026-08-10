#!/usr/bin/env python
# coding: utf-8

# In[9]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


# In[10]:


### Data ###
from Data.simulation_seed_number import *
from Data.simulation_data import *


# In[11]:


### Methods ###
from Scripts.Analysis.oracle_fixedK import *
from Scripts.Analysis.mvFMR_fixedK import *
from Scripts.Analysis.PGM_mvFMR_LASSO_fixedK import *
from Scripts.Analysis.PGM_mvFMR_SCAD_fixedK import *
from Scripts.Analysis.PGM_mvFMR_MCP_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_LASSO_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_SCAD_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_MCP_fixedK import *


# In[12]:


### Analysis function ###
from Scripts.Functions.simulation_summary_functions import *


# In[13]:


def meas_tab(true, pred, idx_name):
    return pd.DataFrame(
        [[TPR(true, pred), FPR(true, pred), FDR(true, pred)] + MSE(true, pred).tolist()],
        index = [idx_name],
        columns = ["TPR", "FPR", "FDR", "MSE_pi", "MSE_Bk", "MSE_sigma"]
    )


# In[14]:


# Assume the number of components (k) is known in Simulation 2.
k = 2


# In[15]:


lambda_s = [0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
            0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
            0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]


# In[16]:


### Model 4 ###
n = 500
# n = 1000

X_521 = data_generate_521(seed_number_521[0], n, indep = True)["X"]
Y_521 = data_generate_521(seed_number_521[0], n, indep = True)["Y"]
true_521 = data_generate_521(seed_number_521[0], n, indep = True)["true"]


# In[ ]:


## Oracle ##
oracle_521 = mvFMR_oracle_fixedK(X_521, Y_521, k=k, true_Bk=true_521["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_521 = mvFMR_fixedK(X_521, Y_521, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_521 = PGM_mvFMR_LASSO_fixedK(X_521, Y_521, k=k, eta_val=1, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_521 = PGM_mvFMR_SCAD_fixedK(X_521, Y_521, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_521 = PGM_mvFMR_MCP_fixedK(X_521, Y_521, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_521 = ADMM_mvFMR_LASSO_fixedK(X_521, Y_521, k=k, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_521 = ADMM_mvFMR_SCAD_fixedK(X_521, Y_521, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[17]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_521 = ADMM_mvFMR_MCP_fixedK(X_521, Y_521, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[18]:


## TPR, FPR, FDR, and MSE ##
tab_521 = pd.concat([meas_tab(true_521, oracle_521["optimal"], "oracle"),
                     meas_tab(true_521, mvFMR_521["optimal"], "mvFMR"),
                     meas_tab(true_521, PGM_mvFMR_LASSO_521["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_521, ADMM_mvFMR_LASSO_521["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_521, PGM_mvFMR_SCAD_521["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_521, ADMM_mvFMR_SCAD_521["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_521, PGM_mvFMR_MCP_521["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_521, ADMM_mvFMR_MCP_521["optimal"], "ADMM_mvFMR_MCP")
                    ])


# In[ ]:


### Model 5 ###
n = 500
# n = 1000

X_522 = data_generate_522(seed_number_522[0], n, indep = True)["X"]
Y_522 = data_generate_522(seed_number_522[0], n, indep = True)["Y"]
true_522 = data_generate_522(seed_number_522[0], n, indep = True)["true"]


# In[ ]:


## Oracle ##
oracle_522 = mvFMR_oracle_fixedK(X_522, Y_522, k=k, true_Bk=true_522["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_522 = mvFMR_fixedK(X_522, Y_522, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.05 when they are correlated.
PGM_mvFMR_LASSO_522 = PGM_mvFMR_LASSO_fixedK(X_522, Y_522, k=k, eta_val=0.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.05 when they are correlated.
PGM_mvFMR_SCAD_522 = PGM_mvFMR_SCAD_fixedK(X_522, Y_522, k=k, eta_val=0.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.05 when they are correlated.
PGM_mvFMR_MCP_522 = PGM_mvFMR_MCP_fixedK(X_522, Y_522, k=k, eta_val=0.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_522 = ADMM_mvFMR_LASSO_fixedK(X_522, Y_522, k=k, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_522 = ADMM_mvFMR_SCAD_fixedK(X_522, Y_522, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_522 = ADMM_mvFMR_MCP_fixedK(X_522, Y_522, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## TPR, FPR, FDR, and MSE ##
tab_522 = pd.concat([meas_tab(true_522, oracle_522["optimal"], "oracle"),
                     meas_tab(true_522, mvFMR_522["optimal"], "mvFMR"),
                     meas_tab(true_522, PGM_mvFMR_LASSO_522["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_522, ADMM_mvFMR_LASSO_522["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_522, PGM_mvFMR_SCAD_522["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_522, ADMM_mvFMR_SCAD_522["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_522, PGM_mvFMR_MCP_522["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_522, ADMM_mvFMR_MCP_522["optimal"], "ADMM_mvFMR_MCP")
                    ])


# In[ ]:


### Model 6 ###
n = 500
# n = 1000

X_523 = data_generate_523(seed_number_523[0], n, indep = True)["X"]
Y_523 = data_generate_523(seed_number_523[0], n, indep = True)["Y"]
true_523 = data_generate_523(seed_number_523[0], n, indep = True)["true"]


# In[ ]:


## Oracle ##
oracle_523 = mvFMR_oracle_fixedK(X_523, Y_523, k=k, true_Bk=true_523["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_523 = mvFMR_fixedK(X_523, Y_523, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.01 when they are correlated.
PGM_mvFMR_LASSO_523 = PGM_mvFMR_LASSO_fixedK(X_523, Y_523, k=k, eta_val=0.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.01 when they are correlated.
PGM_mvFMR_SCAD_523 = PGM_mvFMR_SCAD_fixedK(X_523, Y_523, k=k, eta_val=0.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 0.5 when the predictors are independent and 0.01 when they are correlated.
PGM_mvFMR_MCP_523 = PGM_mvFMR_MCP_fixedK(X_523, Y_523, k=k, eta_val=0.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 0.75 when the predictors are independent and 0.5 when they are correlated.
ADMM_mvFMR_LASSO_523 = ADMM_mvFMR_LASSO_fixedK(X_523, Y_523, k=k, rho=0.75, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 0.75 when the predictors are independent and 0.5 when they are correlated.
ADMM_mvFMR_SCAD_523 = ADMM_mvFMR_SCAD_fixedK(X_523, Y_523, k=k, rho=0.75, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 0.75 when the predictors are independent and 0.5 when they are correlated.
ADMM_mvFMR_MCP_523 = ADMM_mvFMR_MCP_fixedK(X_523, Y_523, k=k, rho=0.75, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## TPR, FPR, FDR, and MSE ##
tab_523 = pd.concat([meas_tab(true_523, oracle_523["optimal"], "oracle"),
                     meas_tab(true_523, mvFMR_523["optimal"], "mvFMR"),
                     meas_tab(true_523, PGM_mvFMR_LASSO_523["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_523, ADMM_mvFMR_LASSO_523["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_523, PGM_mvFMR_SCAD_523["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_523, ADMM_mvFMR_SCAD_523["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_523, PGM_mvFMR_MCP_523["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_523, ADMM_mvFMR_MCP_523["optimal"], "ADMM_mvFMR_MCP")
                    ])

