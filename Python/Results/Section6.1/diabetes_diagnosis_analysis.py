#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import numpy as np
import pandas as pd


# In[3]:


### Data ###
import Data.diabetes_diagnosis_data as dat
import importlib
importlib.reload(dat)


# In[2]:


### Results ###
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

with open(ROOT / "Results/Section6.1/PGM_diabetes_diagnosis.pkl", "rb") as f:
    PGM_diabetes_diagnosis = pickle.load(f)
with open(ROOT / "Results/Section6.1/ADMM_diabetes_diagnosis.pkl", "rb") as f:
    ADMM_diabetes_diagnosis = pickle.load(f)


# In[24]:


### PGM_diabetes_diagnosis ###
PGM_optimal_k = np.nanargmin([np.min(PGM_diabetes_diagnosis["BIC"][k] for k in range(len(PGM_diabetes_diagnosis["BIC"])))]) + 2

## optimal results ##
PGM_optimal = PGM_diabetes_diagnosis["optimal"][PGM_optimal_k-2]

## correlation ##
PGM_corr = []
for i in range(PGM_optimal_k):
    PGM_corr.append(PGM_optimal["sigma"][i][0,1] / np.sqrt(PGM_optimal["sigma"][i][0,0] * PGM_optimal["sigma"][i][1,1]))

## coefficients ##
intercept = pd.DataFrame({"intercept": [1]})
gender = pd.DataFrame({"gender": [1]})
sd_X = pd.concat([intercept, pd.DataFrame(dat.diabetes_X1.std()).T.reset_index(drop=True), gender, pd.DataFrame(dat.diabetes_X2.std()).T.reset_index(drop=True)], axis=1)
PGM_new_Bk = [np.zeros_like(x) for x in PGM_optimal["Bk"]]
sd = sd_X.iloc[0].to_numpy()
for k in range(PGM_optimal_k):
    PGM_new_Bk[k] = PGM_optimal["Bk"][k] / sd[:, None]


# In[91]:


### ADMM_diabetes_diagnosis ###
ADMM_optimal_k = np.nanargmin([np.min(ADMM_diabetes_diagnosis["BIC"][k] for k in range(len(ADMM_diabetes_diagnosis["BIC"])))]) + 2

## optimal results ##
ADMM_optimal = ADMM_diabetes_diagnosis["optimal"][ADMM_optimal_k-2]

## correlation ##
ADMM_corr = []
for i in range(ADMM_optimal_k):
    ADMM_corr.append(ADMM_optimal["sigma"][i][0,1] / np.sqrt(ADMM_optimal["sigma"][i][0,0] * ADMM_optimal["sigma"][i][1,1]))

## coefficients ##
intercept = pd.DataFrame({"intercept": [1]})
gender = pd.DataFrame({"gender": [1]})
sd_X = pd.concat([intercept, pd.DataFrame(dat.diabetes_X1.std()).T.reset_index(drop=True), gender, pd.DataFrame(dat.diabetes_X2.std()).T.reset_index(drop=True)], axis=1)
ADMM_new_Bk = [np.zeros_like(x) for x in ADMM_optimal["Bk"]]
sd = sd_X.iloc[0].to_numpy()
for k in range(ADMM_optimal_k):
    ADMM_new_Bk[k] = ADMM_optimal["Bk"][k] / sd[:, None]

