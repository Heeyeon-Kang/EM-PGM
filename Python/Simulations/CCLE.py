#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

import Data.CCLE_data as dat
import importlib
importlib.reload(dat)


# In[3]:


import numpy as np
import pickle

log_lambdas = np.linspace(np.log10(10**-2), np.log10(1), 30)
lambda_s = 10 ** log_lambdas


# In[10]:


## mvFMR-MCP with the EM-PGM algorithm ##

import Scripts.Analysis.PGM_mvFMR_MCP as pgm
import importlib
importlib.reload(pgm)


# In[ ]:


PGM_CCLE = pgm.PGM_mvFMR_MCP_nonfixedK(dat.X, dat.Y, eta_val=0.01, a=3.7, lambdas=lambda_s)



# In[ ]:


with open(ROOT / "Results/Section6.2/PGM_CCLE.pkl", "wb") as f:
    pickle.dump(PGM_CCLE, f)


# In[ ]:


## 불러올 때 ##
# with open(ROOT / "Result/Section6.2/PGM_CCLE.pkl", "rb") as f:
# 	PGM_CCLE = pickle.load(f)


# In[ ]:


## mvFMR-MCP with the EM-ADMM algorithm ##

import Scripts.Analysis.ADMM_mvFMR_MCP as admm
import importlib
importlib.reload(admm)


# In[ ]:


ADMM_CCLE = admm.ADMM_mvFMR_MCP_nonfixedK(dat.X, dat.Y, rho=0.5, alpha=1.5, a=3.7, lambdas=lambda_s)



# In[ ]:


with open(ROOT / "Results/Section6.2/ADMM_CCLE.pkl", "wb") as f:
    pickle.dump(ADMM_CCLE, f)


# In[ ]:


## 불러올 때 ##
# with open(ROOT / "Result/Section6.2/ADMM_CCLE.pkl", "rb") as f:
# 	ADMM_CCLE = pickle.load(f)
	
	