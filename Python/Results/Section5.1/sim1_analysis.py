#!/usr/bin/env python
# coding: utf-8

# In[31]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


# In[32]:


### Data ###
from Data.simulation_seed_number import *
from Data.simulation_data import *


# In[33]:


### Methods ###
from Scripts.Analysis.oracle_fixedK import *
from Scripts.Analysis.mvFMR_fixedK import *
from Scripts.Analysis.PGM_mvFMR_LASSO_fixedK import *
from Scripts.Analysis.PGM_mvFMR_SCAD_fixedK import *
from Scripts.Analysis.PGM_mvFMR_MCP_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_LASSO_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_SCAD_fixedK import *
from Scripts.Analysis.ADMM_mvFMR_MCP_fixedK import *


# In[34]:


### Analysis function ###
from Scripts.Functions.simulation_summary_functions import *


# In[35]:


def meas_tab(true, pred, idx_name):
    return pd.DataFrame(
        [[TPR(true, pred), FPR(true, pred), FDR(true, pred)] + MSE(true, pred).tolist()],
        index = [idx_name],
        columns = ["TPR", "FPR", "FDR", "MSE_pi", "MSE_Bk", "MSE_sigma"]
    )


# In[36]:


# Assume the number of components (k) is known in Simulation 1.
k = 2


# In[37]:


lambda_s = [0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
            0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
            0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]


# In[38]:


### Model 1 ###
n = 500
# n = 1000

X_511 = data_generate_511(seed_number_511[0], n, indep=True)["X"]
Y_511 = data_generate_511(seed_number_511[0], n, indep=True)["Y"]
true_511 = data_generate_511(seed_number_511[0], n, indep=True)["true"]


# In[ ]:


## Oracle ##
oracle_511 = mvFMR_oracle_fixedK(X_511, Y_511, k=k, true_Bk=true_511["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_511 = mvFMR_fixedK(X_511, Y_511, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_511 = PGM_mvFMR_LASSO_fixedK(X_511, Y_511, k=k, eta_val=1, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_511 = PGM_mvFMR_SCAD_fixedK(X_511, Y_511, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_511 = PGM_mvFMR_MCP_fixedK(X_511, Y_511, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_511 = ADMM_mvFMR_LASSO_fixedK(X_511, Y_511, k=k, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_511 = ADMM_mvFMR_SCAD_fixedK(X_511, Y_511, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[42]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_511 = ADMM_mvFMR_MCP_fixedK(X_511, Y_511, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[43]:


## TPR, FPR, FDR, and MSE ##
tab_511 = pd.concat([meas_tab(true_511, oracle_511["optimal"], "oracle"),
                     meas_tab(true_511, mvFMR_511["optimal"], "mvFMR"),
                     meas_tab(true_511, PGM_mvFMR_LASSO_511["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_511, ADMM_mvFMR_LASSO_511["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_511, PGM_mvFMR_SCAD_511["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_511, ADMM_mvFMR_SCAD_511["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_511, PGM_mvFMR_MCP_511["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_511, ADMM_mvFMR_MCP_511["optimal"], "ADMM_mvFMR_MCP")
                    ])


# In[45]:


### Model 2 ###
n = 500
# n = 1000

X_512 = data_generate_512(seed_number_512[0], n, indep = True)["X"]
Y_512 = data_generate_512(seed_number_512[0], n, indep = True)["Y"]
true_512 = data_generate_512(seed_number_512[0], n, indep = True)["true"]


# In[ ]:


## Oracle ##
oracle_512 = mvFMR_oracle_fixedK(X_512, Y_512, k=k, true_Bk=true_512["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_512 = mvFMR_fixedK(X_512, Y_512, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_512 = PGM_mvFMR_LASSO_fixedK(X_512, Y_512, k=k, eta_val=1, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_512 = PGM_mvFMR_SCAD_fixedK(X_512, Y_512, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_512 = PGM_mvFMR_MCP_fixedK(X_512, Y_512, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_512 = ADMM_mvFMR_LASSO_fixedK(X_512, Y_512, k=k, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_512 = ADMM_mvFMR_SCAD_fixedK(X_512, Y_512, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[46]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_512 = ADMM_mvFMR_MCP_fixedK(X_512, Y_512, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[47]:


## TPR, FPR, FDR, and MSE ##
tab_512 = pd.concat([meas_tab(true_512, oracle_512["optimal"], "oracle"),
                     meas_tab(true_512, mvFMR_512["optimal"], "mvFMR"),
                     meas_tab(true_512, PGM_mvFMR_LASSO_512["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_512, ADMM_mvFMR_LASSO_512["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_512, PGM_mvFMR_SCAD_512["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_512, ADMM_mvFMR_SCAD_512["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_512, PGM_mvFMR_MCP_512["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_512, ADMM_mvFMR_MCP_512["optimal"], "ADMM_mvFMR_MCP")
                    ])


# In[49]:


### Model 3 ###
n = 500
# n = 1000

X_513 = data_generate_513(seed_number_513[0], n, indep = True)["X"]
Y_513 = data_generate_513(seed_number_513[0], n, indep = True)["Y"]
true_513 = data_generate_513(seed_number_513[0], n, indep = True)["true"]


# In[ ]:


## Oracle ##
oracle_513 = mvFMR_oracle_fixedK(X_513, Y_513, k=k, true_Bk=true_513["Bk"])


# In[ ]:


## mvFMR ##
mvFMR_513 = mvFMR_fixedK(X_513, Y_513, k=k)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_513 = PGM_mvFMR_LASSO_fixedK(X_513, Y_513, k=k, eta_val=1, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_513 = PGM_mvFMR_SCAD_fixedK(X_513, Y_513, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_513 = PGM_mvFMR_MCP_fixedK(X_513, Y_513, k=k, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_513 = ADMM_mvFMR_LASSO_fixedK(X_513, Y_513, k=k, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_513 = ADMM_mvFMR_SCAD_fixedK(X_513, Y_513, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[50]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_513 = ADMM_mvFMR_MCP_fixedK(X_513, Y_513, k=k, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[51]:


## TPR, FPR, FDR, and MSE ##
tab_513 = pd.concat([meas_tab(true_513, oracle_513["optimal"], "oracle"),
                     meas_tab(true_513, mvFMR_513["optimal"], "mvFMR"),
                     meas_tab(true_513, PGM_mvFMR_LASSO_513["optimal"], "PGM_mvFMR_LASSO"),
                     meas_tab(true_513, ADMM_mvFMR_LASSO_513["optimal"], "ADMM_mvFMR_LASSO"),
                     meas_tab(true_513, PGM_mvFMR_SCAD_513["optimal"], "PGM_mvFMR_SCAD"),
                     meas_tab(true_513, ADMM_mvFMR_SCAD_513["optimal"], "ADMM_mvFMR_SCAD"),
                     meas_tab(true_513, PGM_mvFMR_MCP_513["optimal"], "PGM_mvFMR_MCP"),
                     meas_tab(true_513, ADMM_mvFMR_MCP_513["optimal"], "ADMM_mvFMR_MCP")
                    ])

