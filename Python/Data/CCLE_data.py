#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

X = pd.read_csv(ROOT / "Data/CCLE_genes.csv", encoding="utf-8")
Y = pd.read_csv(ROOT / "Data/CCLE_drugs.csv", encoding="utf-8")


# In[3]:


C = np.corrcoef(X.T, Y.T)[:X.shape[1], X.shape[1]:]


# In[4]:


score = np.sum(np.abs(C), axis=1)


# In[5]:


ord = np.argsort(-score)


# In[6]:


Xr = X.iloc[:, ord[:300]]
top_genes = X.columns[ord[:300]].tolist()


# In[7]:


X_scaled = (Xr - Xr.mean()) / Xr.std(ddof=1)
X = pd.concat(
    [pd.Series(1, index=Xr.index, name="intercept"), X_scaled],
    axis=1
)


# In[8]:


X_sd = Xr.std(ddof=1)

