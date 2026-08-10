#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np


# In[15]:


### Section 5.1 : model 1 ###
seed_511 = np.random.default_rng(511)
seed_number_511 = seed_511.choice(np.arange(1, 100001), size = 500, replace = False)


# In[5]:


### Section 5.1 : model 2 ###
seed_512 = np.random.default_rng(512)
seed_number_512 = seed_512.choice(np.arange(1, 100001), size = 500, replace = False)


# In[6]:


### Section 5.1 : model 3 ###
seed_513 = np.random.default_rng(513)
seed_number_513 = seed_513.choice(np.arange(1, 100001), size = 500, replace = False)


# In[9]:


### Section 5.2 : model 1 ###
seed_521 = np.random.default_rng(521)
seed_number_521 = seed_521.choice(np.arange(1, 100001), size = 500, replace = False)


# In[10]:


### Section 5.2 : model 2 ###
seed_522 = np.random.default_rng(522)
seed_number_522 = seed_522.choice(np.arange(1, 100001), size = 500, replace = False)


# In[11]:


### Section 5.2 : model 3 ###
seed_523 = np.random.default_rng(523)
seed_number_523 = seed_523.choice(np.arange(1, 100001), size = 500, replace = False)


# In[12]:


### Section 5.3 : model 1 ###
seed_531 = np.random.default_rng(531)
seed_number_531 = seed_531.choice(np.arange(1, 100001), size = 500, replace = False)


# In[13]:


### Section 5.3 : model 2 ###
seed_532 = np.random.default_rng(532)
seed_number_532 = seed_532.choice(np.arange(1, 100001), size = 500, replace = False)

