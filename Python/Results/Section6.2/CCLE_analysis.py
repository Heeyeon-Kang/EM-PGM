#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pickle
import numpy as np
import pandas as pd


# In[ ]:


from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))


# In[4]:


### Data ###
import Data.CCLE_data as dat
import importlib
importlib.reload(dat)


# In[ ]:


### Results ###
with open(ROOT / "Results/Section6.2/PGM_CCLE.pkl", "rb") as f:
    PGM_CCLE = pickle.load(f)
with open(ROOT / "Results/Section6.2/ADMM_CCLE.pkl", "rb") as f:
    ADMM_CCLE = pickle.load(f)


# In[5]:


### PGM_CCLE ###
PGM_optimal_k = np.nanargmin([np.min(PGM_CCLE["BIC"][k] for k in range(len(PGM_CCLE["BIC"])))]) + 2

## optimal results ##
PGM_optimal = PGM_CCLE["optimal"][PGM_optimal_k-2]

## correlation ##
row_names = ["intercept"] + dat.top_genes.tolist()
col_names = ["Erlotinib", "AZD6244", "PD-0325901"]
PGM_corr = [
    pd.DataFrame(
        np.eye(3),
        index = col_names,
        columns = col_names
    )
    for _ in range(2)
]
for k in range(PGM_optimal_k):
    for i in range(3):
        for j in range(3):
            if i != j:
                PGM_corr[k].iloc[i,j] = PGM_optimal["sigma"][k][i,j] / np.sqrt(PGM_optimal["sigma"][k][i,i] * PGM_optimal["sigma"][k][j,j])

## coefficients ##
intercept = pd.DataFrame({"intercept": [1]})
sd_X = pd.concat([intercept, pd.DataFrame(dat.X_sd).T.reset_index(drop=True)], axis=1)
new_Bk = [
    pd.DataFrame(
        np.zeros_like(PGM_optimal["Bk"][k]),
        index = row_names,
        columns = col_names
    )
    for k in range(PGM_optimal_k)
]

sd = sd_X.iloc[0].to_numpy()
PGM_new_Bk = new_Bk.copy()
for k in range(PGM_optimal_k):
    PGM_new_Bk[k] = pd.DataFrame(
        PGM_optimal["Bk"][k] / sd[:, None],
        index = row_names,
        columns = col_names
    )

## active coefficients ##
PGM_active_Bk = []
for k in range(PGM_optimal_k):
    mask = (PGM_new_Bk[k].abs() > 1e-6).sum(axis=1) > 0
    PGM_active_Bk.append(PGM_new_Bk[k].loc[mask])


# In[30]:


### ADMM_CCLE ###
ADMM_optimal_k = np.nanargmin([np.min(ADMM_CCLE["BIC"][k] for k in range(len(ADMM_CCLE["BIC"])))]) + 2

## optimal results ##
ADMM_optimal = ADMM_CCLE["optimal"][ADMM_optimal_k-2]

## correlation ##
row_names = ["intercept"] + dat.top_genes.tolist()
col_names = ["Erlotinib", "AZD6244", "PD-0325901"]
ADMM_corr = [
    pd.DataFrame(
        np.eye(3),
        index = col_names,
        columns = col_names
    )
    for _ in range(2)
]
for k in range(ADMM_optimal_k):
    for i in range(3):
        for j in range(3):
            if i != j:
                ADMM_corr[k].iloc[i,j] = ADMM_optimal["sigma"][k][i,j] / np.sqrt(ADMM_optimal["sigma"][k][i,i] * ADMM_optimal["sigma"][k][j,j])

## coefficients ##
intercept = pd.DataFrame({"intercept": [1]})
sd_X = pd.concat([intercept, pd.DataFrame(dat.X_sd).T.reset_index(drop=True)], axis=1)
new_Bk = [
    pd.DataFrame(
        np.zeros_like(ADMM_optimal["Bk"][k]),
        index = row_names,
        columns = col_names
    )
    for k in range(ADMM_optimal_k)
]

sd = sd_X.iloc[0].to_numpy()
ADMM_new_Bk = new_Bk.copy()
for k in range(ADMM_optimal_k):
    ADMM_new_Bk[k] = pd.DataFrame(
        ADMM_optimal["Bk"][k] / sd[:, None],
        index = row_names,
        columns = col_names
    )

## active coefficients ##
ADMM_active_Bk = []
for k in range(ADMM_optimal_k):
    Bk = ADMM_new_Bk[k].copy()
    Bk[Bk.abs() < 1e-6] = 0
    mask = (Bk != 0).any(axis=1)
    ADMM_active_Bk.append(Bk.loc[mask])

