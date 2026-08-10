#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

import Data.diabetes_diagnosis_data as dat
import importlib
importlib.reload(dat)



# In[10]:


import numpy as np
import pickle

log_lambdas = np.linspace(np.log10(10**-4), np.log10(1), 30)
lambda_s = 10 ** log_lambdas


# In[28]:


## mvFMR-MCP with the EM-PGM algorithm ##

import Scripts.Analysis.PGM_mvFMR_MCP as pgm
import importlib
importlib.reload(pgm)


# In[12]:


PGM_diabetes_diagnosis = pgm.PGM_mvFMR_MCP_nonfixedK(dat.X, dat.Y, eta_val=0.25, a=3.7, lambdas=lambda_s)


# In[ ]:


with open(ROOT / "Results/Section6.1/PGM_diabetes_diagnosis.pkl", "wb") as f:
    pickle.dump(PGM_diabetes_diagnosis, f)



# In[22]:


## 불러올 때 ##
# with open(ROOT / "Results/Section6.1/PGM_diabetes_diagnosis.pkl", "rb") as f:
#     PGM_diabetes_diagnosis = pickle.load(f)


# In[24]:


## mvFMR-MCP with the EM-ADMM algorithm ##

import Scripts.Analysis.ADMM_mvFMR_MCP as admm
import importlib
importlib.reload(admm)


# In[26]:


ADMM_diabetes_diagnosis = admm.ADMM_mvFMR_MCP_nonfixedK(dat.X, dat.Y, rho=1, alpha=1.5, a=3.7, lambdas=lambda_s)


# In[ ]:


with open(ROOT / "Results/Section6.1/ADMM_diabetes_diagnosis.pkl", "wb") as f:
    pickle.dump(ADMM_diabetes_diagnosis, f)




# In[ ]:


## 불러올 때 ##
# with open(ROOT / "Results/Section6.1/ADMM_diabetes_diagnosis.pkl", "rb") as f:
#     ADMM_diabetes_diagnosis = pickle.load(f)



