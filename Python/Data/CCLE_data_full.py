#!/usr/bin/env python
# coding: utf-8

# In[145]:


import pandas as pd
import numpy as np


# In[146]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

data = pd.read_csv(ROOT / "Data/CCLE_expression.csv", encoding="utf-8")
response = pd.read_csv(ROOT / "Data/CCLE_Drugdata.csv", encoding="utf-8")



# In[134]:


w = data.columns[2:]
gene = data.loc[0:]["Description"]
Xgenes = data.loc[0:].drop(['Name', 'Description'], axis=1)


# In[135]:


Xgenes.columns = w
Xgenes.index = gene


# In[136]:


z1 = response[response["Compound"] == "Erlotinib"].loc[:, ["CCLE Cell Line Name"]]
z2 = response[response["Compound"] == "AZD6244"].loc[:, ["CCLE Cell Line Name"]]
z3 = response[response["Compound"] == "PD-0325901"].loc[:, ["CCLE Cell Line Name"]]


# In[137]:


area1 = response[response["Compound"] == "Erlotinib"].loc[:, ["ActArea"]]
area2 = response[response["Compound"] == "AZD6244"].loc[:, ["ActArea"]]
area3 = response[response["Compound"] == "PD-0325901"].loc[:, ["ActArea"]]


# In[138]:


## Common cell line for three drugs ##
z = list(
    set(z1["CCLE Cell Line Name"])
    & set(z2["CCLE Cell Line Name"])
    & set(z3["CCLE Cell Line Name"])
)


# In[139]:


tmp1 = pd.concat(
    [z1.reset_index(drop=True),
     area1.reset_index(drop=True)],
    axis=1
)
tmp1.columns = ["CCLE Cell Line Name", "Erlotinib"]

tmp2 = pd.concat(
    [z2.reset_index(drop=True),
     area2.reset_index(drop=True)],
    axis=1
)
tmp2.columns = ["CCLE Cell Line Name", "AZD6244"]

tmp3 = pd.concat(
    [z3.reset_index(drop=True),
     area3.reset_index(drop=True)],
    axis=1
)
tmp3.columns = ["CCLE Cell Line Name", "PD-0325901"]


# In[140]:


area = pd.DataFrame({"CCLE Cell Line Name": z})

area = (
    area
    .merge(tmp1, on="CCLE Cell Line Name", how="left")
    .merge(tmp2, on="CCLE Cell Line Name", how="left")
    .merge(tmp3, on="CCLE Cell Line Name", how="left")
)

area = area.set_index("CCLE Cell Line Name")


# In[141]:


w_key = Xgenes.columns.astype(str).str.upper().str.strip()
z_key = area.index.astype(str).str.upper().str.strip()

J_all = w_key.get_indexer(z_key)

I = np.where(J_all != -1)[0]
J = J_all[I]


# In[142]:


## pandas.core.frame.DataFrame ##
df_Y = area.iloc[I, :].copy()
df_X = Xgenes.iloc[:, J].T.copy()


# In[143]:


C_full = np.corrcoef(
    df_X.to_numpy().T,
    df_Y.to_numpy().T
)

p = df_X.shape[1]
m = df_Y.shape[1]

C = C_full[:p, p:p+m]

score = np.sum(np.abs(C), axis=1)

ord = np.argsort(score)[::-1]

## top 300 genes ##
top_idx = ord[:300]

Xr = df_X.iloc[:, top_idx].copy()
top_genes = df_X.columns[top_idx]

X_sd = Xr.std()
X_scaled = (Xr - Xr.mean()) / X_sd

## add intercept ##
df_X = pd.concat(
    [
        pd.DataFrame(
            {"intercept": np.ones(len(Xr))},
            index=Xr.index
        ),
        X_scaled
    ],
    axis=1
)


# In[144]:


## numpy.ndarray ##
X = np.array(df_X)
Y = np.array(df_Y)

