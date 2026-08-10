#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

import pandas as pd
import numpy as np


# In[ ]:


diabetes_dat = pd.read_csv(ROOT / "Data/diabetes_diagnosis_data.csv", encoding="utf-8")

delete_col = ['id', 'ratio', 'location', 'frame', 'bp.2s', 'bp.2d']
diabetes_dat1 = diabetes_dat.drop(delete_col, axis=1).dropna()
gender = diabetes_dat1['gender'].map({'male': 0, 'female': 1})


# In[138]:


Y_col = ['stab.glu', 'glyhb']
diabetes_Y = diabetes_dat1.loc[:, Y_col]
diabetes_Y.columns = ['glucose', 'HbA1c']

# df_Y: data frame #
df_Y = diabetes_Y

# Y: array #
Y = np.array(df_Y)


# In[139]:


X_col1 = ['chol', 'hdl', 'age']
X_col2 = ['height', 'weight', 'bp.1s', 'bp.1d', 'waist', 'hip']
diabetes_X1 = diabetes_dat1.loc[:, X_col1]
diabetes_X2 = diabetes_dat1.loc[:, X_col2]
diabetes_X1_scaled = (diabetes_X1 - diabetes_X1.mean()) / diabetes_X1.std(ddof=1)
diabetes_X2_scaled = (diabetes_X2 - diabetes_X2.mean()) / diabetes_X2.std(ddof=1)

# df_X: data frame #
df_X = pd.concat(
    [
        pd.Series(1, index=gender.index, name="intercept"),
        diabetes_X1_scaled,
        gender,
        diabetes_X2_scaled
    ],
    axis=1
)

# X: array #
X = np.array(df_X)

