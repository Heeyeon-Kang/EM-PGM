#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


# In[2]:


### Data ###
from Data.simulation_seed_number import *
from Data.simulation_data import *


# In[3]:


### Methods ###
from Scripts.Analysis.mvFMR import *
from Scripts.Analysis.PGM_mvFMR_LASSO import *
from Scripts.Analysis.PGM_mvFMR_SCAD import *
from Scripts.Analysis.PGM_mvFMR_MCP import *
from Scripts.Analysis.ADMM_mvFMR_LASSO import *
from Scripts.Analysis.ADMM_mvFMR_SCAD import *
from Scripts.Analysis.ADMM_mvFMR_MCP import *


# In[4]:


### Analysis function ###
from Scripts.Functions.simulation_summary_functions import *


# In[6]:


def meas_pred_ll(dat, pred, n, idx_name):
    dens = pred_ll_dens(dat, pred)
    return pd.DataFrame(
        [predictive_ll(pred, dens, n)],
        index = [idx_name],
        columns = ["predictive_ll"]
    )


# In[7]:


lambda_s = [0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
            0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
            0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]


# In[8]:


### Model 7 ###
n = 500
# n = 1000

dat_531 = data_generate_531(seed_number_531[0], n, indep=True)
X_531 = dat_531["train_X"]
Y_531 = dat_531["train_Y"]
true_531 = dat_531["true"]


# In[9]:


## mvFMR ##
mvFMR_531 = mvFMR_nonfixedK(X_531, Y_531)


# In[11]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_531 = PGM_mvFMR_LASSO_nonfixedK(X_531, Y_531, eta_val=1, lambdas=lambda_s)


# In[13]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_531 = PGM_mvFMR_SCAD_nonfixedK(X_531, Y_531, eta_val=1, a=3.7, lambdas=lambda_s)


# In[15]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_531 = PGM_mvFMR_MCP_nonfixedK(X_531, Y_531, eta_val=1, a=3.7, lambdas=lambda_s)


# In[17]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_531 = ADMM_mvFMR_LASSO_nonfixedK(X_531, Y_531, rho=1, alpha=1.5, lambdas=lambda_s)


# In[19]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_531 = ADMM_mvFMR_SCAD_nonfixedK(X_531, Y_531, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[22]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_531 = ADMM_mvFMR_MCP_nonfixedK(X_531, Y_531, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[24]:


## TPR, FPR, FDR, and MSE ##
tab_531 = pd.concat([meas_pred_ll(dat_531, mvFMR_531["optimal"], n, "mvFMR"),
                     meas_pred_ll(dat_531, PGM_mvFMR_LASSO_531["optimal"], n, "PGM_mvFMR_LASSO"),
                     meas_pred_ll(dat_531, ADMM_mvFMR_LASSO_531["optimal"], n, "ADMM_mvFMR_LASSO"),
                     meas_pred_ll(dat_531, PGM_mvFMR_SCAD_531["optimal"], n, "PGM_mvFMR_SCAD"),
                     meas_pred_ll(dat_531, ADMM_mvFMR_SCAD_531["optimal"], n, "ADMM_mvFMR_SCAD"),
                     meas_pred_ll(dat_531, PGM_mvFMR_MCP_531["optimal"], n, "PGM_mvFMR_MCP"),
                     meas_pred_ll(dat_531, ADMM_mvFMR_MCP_531["optimal"], n, "ADMM_mvFMR_MCP"),
                    ])


# In[ ]:


### Model 8 ###
n = 500
# n = 1000

dat_532 = data_generate_532(seed_number_532[0], n, indep=True)
X_532 = dat_532["train_X"]
Y_532 = dat_532["train_Y"]
true_532 = dat_532["true"]


# In[ ]:


## mvFMR ##
mvFMR_532 = mvFMR_nonfixedK(X_532, Y_532)


# In[ ]:


## mvFMR-LASSO using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_LASSO_532 = PGM_mvFMR_LASSO_nonfixedK(X_532, Y_532, eta_val=1, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_SCAD_532 = PGM_mvFMR_SCAD_nonfixedK(X_532, Y_532, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-PGM algorithm ##
# We used eta_val = 1 when the predictors are independent and 0.25 when they are correlated.
PGM_mvFMR_MCP_532 = PGM_mvFMR_MCP_nonfixedK(X_532, Y_532, eta_val=1, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-LASSO using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_LASSO_532 = ADMM_mvFMR_LASSO_nonfixedK(X_532, Y_532, rho=1, alpha=1.5, lambdas=lambda_s)


# In[ ]:


## mvFMR-SCAD using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_SCAD_532 = ADMM_mvFMR_SCAD_nonfixedK(X_532, Y_532, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## mvFMR-MCP using the EM-ADMM algorithm ##
# We used rho = 1 in both predictor designs; independent and correlated.
ADMM_mvFMR_MCP_532 = ADMM_mvFMR_MCP_nonfixedK(X_532, Y_532, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


## TPR, FPR, FDR, and MSE ##
tab_532 = pd.concat([meas_pred_ll(dat_532, mvFMR_532["optimal"], n, "mvFMR"),
                     meas_pred_ll(dat_532, PGM_mvFMR_LASSO_532["optimal"], n, "PGM_mvFMR_LASSO"),
                     meas_pred_ll(dat_532, ADMM_mvFMR_LASSO_532["optimal"], n, "ADMM_mvFMR_LASSO"),
                     meas_pred_ll(dat_532, PGM_mvFMR_SCAD_532["optimal"], n, "PGM_mvFMR_SCAD"),
                     meas_pred_ll(dat_532, ADMM_mvFMR_SCAD_532["optimal"], n, "ADMM_mvFMR_SCAD"),
                     meas_pred_ll(dat_532, PGM_mvFMR_MCP_532["optimal"], n, "PGM_mvFMR_MCP"),
                     meas_pred_ll(dat_532, ADMM_mvFMR_MCP_532["optimal"], n, "ADMM_mvFMR_MCP"),
                    ])

