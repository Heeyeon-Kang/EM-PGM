#!/usr/bin/env python
# coding: utf-8

# In[206]:


import numpy as np
import math


# In[105]:


### Section 5.1 : model 1 ###
### m = 2, p = 10, K = 2  ###
def data_generate_511(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    true_Bk = [
        np.array([[ 1, -1], 
                  [ 0,  2],
                  [ 0,  0], 
                  [ 3,  0], 
                  [-2,  3],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0]]),
        np.array([[-1,  1],
                  [ 3,  0],
                  [ 0,  2],
                  [ 1, -1],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0]])
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[104]:


### Section 5.1 : model 2 ###
### m = 5, p = 10, K = 2  ###
def data_generate_512(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    sigma = np.eye(5)
    sigma[sigma == 0] = 0.5
    true_sigma = [sigma, sigma]

    ## The coefficients matrices of each components ##
    true_Bk = [
        np.array([[ 1, -1,  1,  2, -1], 
                  [ 0,  2,  0,  1,  0],
                  [ 0,  0, -2, -1,  3], 
                  [ 3,  0,  0,  0,  1], 
                  [-2,  3,  0,  1,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0]]),
        np.array([[-1,  1,  2, -2,  1],
                  [ 3,  0,  0,  1,  0],
                  [ 0,  2,  1,  0,  1],
                  [ 1, -1,  0,  0,  0],
                  [ 0,  0,  3,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0]])
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[106]:


### Section 5.1 : model 3 ###
### m = 10, p = 10, K = 2 ###
def data_generate_513(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    sigma = np.eye(10)
    sigma[sigma == 0] = 0.5
    true_sigma = [sigma, sigma]

    ## The coefficients matrices of each components ##
    true_Bk = [
        np.array([[ 1, -1,  1,  2, -1, -2,  1,  2, -3, -1], 
                  [ 0,  2,  0,  1,  0,  0, -3,  0,  1,  0],
                  [ 0,  0, -2, -1,  3,  1, -1,  0, -1,  0], 
                  [ 3,  0,  0,  0,  1, -1,  2,  0,  0,  1], 
                  [-2,  3,  0,  1,  0,  0,  0, -1,  2, -3],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0]]),
        np.array([[-1,  1,  2, -2,  1,  3, -3,  2, -1, -2],
                  [ 3,  0,  0,  1,  0,  0,  2,  0,  2,  0],
                  [ 0,  2,  1,  0,  1,  0, -1,  0,  3,  0],
                  [ 1, -1,  0,  0,  0, -1,  1, -1,  0,  0],
                  [ 0,  0,  3,  0,  0,  0,  0,  3,  0,  1],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0]])
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[134]:


### Section 5.2 : model 4 ###
### m = 2, p = 30, K = 2 ###
def data_generate_521(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    zero_mat = np.zeros((20, 2))
    true_Bk = [
        np.vstack((
            np.array([[   1, -0.5], 
                      [-0.5,    0],
                      [   0,  1.5], 
                      [ 1.5,  0.5], 
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [-0.5,    0],
                      [   0, -1.5]]), zero_mat)),
        np.vstack((
            np.array([[  -1,  0.5],
                      [   0,  1.5],
                      [ 0.5, -0.5],
                      [-1.5,    0],
                      [   0,    1],
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [ 1.5,    0],
                      [-0.5,    1]]), zero_mat))
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[140]:


### Section 5.2 : model 5 ###
### m = 2, p = 100, K = 2 ###
def data_generate_522(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    zero_mat = np.zeros((90, 2))
    true_Bk = [
        np.vstack((
            np.array([[   1, -0.5], 
                      [-0.5,    0],
                      [   0,  1.5], 
                      [ 1.5,  0.5], 
                      [   0,    0], 
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [-0.5,    0],
                      [   0, -1.5]]), zero_mat)),
        np.vstack((
            np.array([[  -1,  0.5],
                      [   0,  1.5],
                      [ 0.5, -0.5],
                      [-1.5,    0],
                      [   0,    1],
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [ 1.5,    0],
                      [-0.5,    1]]), zero_mat))
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[141]:


### Section 5.2 : model 6 ###
### m = 2, p = 300, K = 2 ###
def data_generate_523(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    zero_mat = np.zeros((290, 2))
    true_Bk = [
        np.vstack((
            np.array([[   1, -0.5], 
                      [-0.5,    0],
                      [   0,  1.5], 
                      [ 1.5,  0.5], 
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [-0.5,    0],
                      [   0, -1.5]]), zero_mat)),
        np.vstack((
            np.array([[  -1,  0.5],
                      [   0,  1.5],
                      [ 0.5, -0.5],
                      [-1.5,    0],
                      [   0,    1],
                      [   0,    0],
                      [   0,    0],
                      [   0,    0],
                      [ 1.5,    0],
                      [-0.5,    1]]), zero_mat))
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result
    else:
        p = true_Bk[0].shape[0] - 1
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = size
        )
        X_data = np.column_stack((np.ones(size), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = size),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = size)
        ]

        cluster = rng.choice(2, size = size, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(size) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

    return {
        "X": X_data,
        "Y": Y_data,
        "cluster": cluster,
        "true": true
    }


# In[195]:


### Section 5.3 : model 7 ###
### m = 2, p = 10, K = 2 ###
def data_generate_531(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.5, 0.5]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    true_Bk = [
        np.array([[-5,  2], 
                  [ 2, -4],
                  [ 4,  3], 
                  [ 0,  0], 
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0]]),
        np.array([[ 3,  4],
                  [ 0,  0],
                  [ 0,  0],
                  [-3,  5],
                  [ 0,  0],
                  [ 0,  0],
                  [ 4, -3],
                  [ 5,  2],
                  [ 0,  0],
                  [ 0,  0]])
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        sz = math.floor(size * 1.4)
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = sz
        )
        X_data = np.column_stack((np.ones(sz), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = sz)
        ]

        cluster = rng.choice(2, size = sz, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(sz) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

        train_X = X_data[:size]
        test_X = X_data[size:]
        train_Y = Y_data[:size]
        test_Y = Y_data[size:]
        train_cluster = cluster[:size]
        test_cluster = cluster[size:]  
    else:
        p = true_Bk[0].shape[0] - 1
        sz = math.floor(size * 1.4)
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = sz
        )
        X_data = np.column_stack((np.ones(sz), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = sz)
        ]

        cluster = rng.choice(2, size = sz, replace = True, p = true_pi_)
        rb1 = cluster.copy()
        rb2 = np.ones(sz) - rb1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        Y_data = rb1_result + rb2_result

        train_X = X_data[:size]
        test_X = X_data[size:]
        train_Y = Y_data[:size]
        test_Y = Y_data[size:]
        train_cluster = cluster[:size]
        test_cluster = cluster[size:]  

    return {
        "train_X": train_X,
        "train_Y": train_Y,
        "train_cluster": train_cluster,
        "test_X": test_X,
        "test_Y": test_Y,
        "test_cluster": test_cluster,
        "true": true
    }


# In[196]:


### Section 5.3 : model 8 ###
### m = 2, p = 10, K = 4 ###
def data_generate_532(z, size, indep=True):

    ## The mixing proportions ##
    true_pi_ = [0.25, 0.25, 0.25, 0.25]

    ## The variance-covariance matrices of errors ##
    true_sigma = [
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]]),
        np.array([[1, 0.5],
                  [0.5, 1]])
    ]

    ## The coefficients matrices of each components ##
    true_Bk = [
        np.array([[-5,  2], 
                  [ 2, -4],
                  [ 4,  3], 
                  [ 0,  0], 
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0]]),
        np.array([[ 3,  4],
                  [ 0,  0],
                  [ 0,  0],
                  [-3,  5],
                  [ 0,  0],
                  [ 0,  0],
                  [ 4, -3],
                  [ 5,  2],
                  [ 0,  0],
                  [ 0,  0]]),
        np.array([[ 2, -4], 
                  [ 0,  0],
                  [ 0,  0], 
                  [ 0,  0], 
                  [-3,  3],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 5,  2],
                  [ 0,  0]]),
        np.array([[-4,  3],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [-2,  4],
                  [ 0,  0],
                  [ 0,  0],
                  [ 0,  0],
                  [ 3, -5]])
    ]

    true = {
        "pi": true_pi_,
        "Bk": true_Bk,
        "sigma": true_sigma
    }

    rng = np.random.default_rng(z)
    if indep:
        p = true_Bk[0].shape[0] - 1
        sz = math.floor(size * 1.4)
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = np.eye(p),
            size = sz
        )
        X_data = np.column_stack((np.ones(sz), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[2],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[3],
                size = sz)
        ]

        cluster = rng.choice(4, size = sz, replace = True, p = true_pi_)
        rb1 = np.zeros(sz)
        rb2 = np.zeros(sz)
        rb3 = np.zeros(sz)
        rb4 = np.zeros(sz)

        rb1[cluster == 0] = 1
        rb2[cluster == 1] = 1
        rb3[cluster == 2] = 1
        rb4[cluster == 3] = 1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb3_result = (X_data @ true_Bk[2] + error[2]).copy()
        rb4_result = (X_data @ true_Bk[3] + error[3]).copy()

        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        rb3_result[rb3 == 0] = 0
        rb4_result[rb4 == 0] = 0
        Y_data = rb1_result + rb2_result + rb3_result + rb4_result

        train_X = X_data[:size]
        test_X = X_data[size:]
        train_Y = Y_data[:size]
        test_Y = Y_data[size:]
        train_cluster = cluster[:size]
        test_cluster = cluster[size:]  
    else:
        p = true_Bk[0].shape[0] - 1
        sz = math.floor(size * 1.4)
        cov = np.eye(p)
        cov[cov == 0] = 0.5
        x_data = rng.multivariate_normal(
            mean = np.zeros(p),
            cov = cov,
            size = sz
        )
        X_data = np.column_stack((np.ones(sz), x_data))

        m = true_sigma[0].shape[0]
        error = [
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[0],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[1],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[2],
                size = sz),
            rng.multivariate_normal(
                mean = np.zeros(m),
                cov = true_sigma[3],
                size = sz)
        ]

        cluster = rng.choice(4, size = sz, replace = True, p = true_pi_)
        rb1 = np.zeros(sz)
        rb2 = np.zeros(sz)
        rb3 = np.zeros(sz)
        rb4 = np.zeros(sz)

        rb1[cluster == 0] = 1
        rb2[cluster == 1] = 1
        rb3[cluster == 2] = 1
        rb4[cluster == 3] = 1

        rb1_result = (X_data @ true_Bk[0] + error[0]).copy()
        rb2_result = (X_data @ true_Bk[1] + error[1]).copy()
        rb3_result = (X_data @ true_Bk[2] + error[2]).copy()
        rb4_result = (X_data @ true_Bk[3] + error[3]).copy()

        rb1_result[rb1 == 0] = 0
        rb2_result[rb2 == 0] = 0
        rb3_result[rb3 == 0] = 0
        rb4_result[rb4 == 0] = 0
        Y_data = rb1_result + rb2_result + rb3_result + rb4_result

        train_X = X_data[:size]
        test_X = X_data[size:]
        train_Y = Y_data[:size]
        test_Y = Y_data[size:]
        train_cluster = cluster[:size]
        test_cluster = cluster[size:]  

    return {
        "train_X": train_X,
        "train_Y": train_Y,
        "train_cluster": train_cluster,
        "test_X": test_X,
        "test_Y": test_Y,
        "test_cluster": test_cluster,
        "true": true
    }

