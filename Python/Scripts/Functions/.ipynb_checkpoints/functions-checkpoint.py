#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
from scipy import linalg


# In[3]:


def density_f(X, Y, B, sigma, inv_sigma):
    n = len(X)
    m = Y.shape[1]

    a = (2 * np.pi) ** (-m / 2)
    b = np.linalg.det(sigma) ** (-1/2)
    e = np.array([])
    for i in range(n):
        res = Y[i] - (B.T @ X[i])
        e = np.append(e, (-1/2) * (res.T @ inv_sigma @ res))

    return a * b * np.exp(e)


# In[3]:


def safe_regularize_sigma(S, alpha = 0.1, epsilon = 1e-3, prior = None):
    if prior is None:
        prior = np.eye(S.shape[1])

    if not np.all(np.isfinite(S)):
        return prior

    S = alpha * prior + (1 - alpha) * S

    eigvals, eigvecs = np.linalg.eigh(S)

    if not np.all(np.isfinite(eigvals)):
        return prior

    eigvals_safe = np.maximum(eigvals, epsilon)

    S_psd = eigvecs @ np.diag(eigvals_safe) @ eigvecs.T

    return S_psd

def safe_log(x):
    return np.log(np.maximum(x, np.finfo(float).tiny))
    

# In[4]:


def gradient_obj(X, Y, wk, Bk, inv_sigma_k):
    n = len(X)
    resid = Y - X @ Bk

    Z = resid @ inv_sigma_k
    Xw = wk * X.T
    grad = (-1/n) * (Xw @ Z)

    return grad


# In[5]:


def diff_norm(B, B_new):
    k = len(B)
    diff = []
    for i in range(k):
        diff.append(np.linalg.norm(B_new[i] - B[i], ord = 'fro'))

    return diff

def stopping_criteria(R, S, E):
    k = len(R)
    R_norm = [
        np.linalg.norm(R[j], ord = 'fro')
        for j in range(k)
    ]
    S_norm = [
        np.linalg.norm(S[j], ord = 'fro')
        for j in range(k)
    ]

    return any(R_norm > E[:, 0]) or any(S_norm > E[:, 1])


# In[6]:


def LASSO_thresholding(BB, p, lam, step_size):
    threshold = step_size * p * lam
    Bk_new = BB.copy()
    in_thres = (abs(BB) <= threshold)

    Bk_new = BB - (np.sign(BB) * threshold)
    Bk_new = np.where(in_thres, 0, Bk_new)

    return Bk_new

def SCAD_thresholding(BB, p, lam, step_size, a):
    threshold1 = step_size * p * lam
    threshold2 = (1 + step_size * p) * lam
    threshold3 = a * lam

    Bk_new = BB.copy()
    in_thres1 = (abs(BB) <= threshold1)
    in_thres2 = np.logical_and(abs(BB) > threshold1, abs(BB) <= threshold2)
    in_thres3 = np.logical_and(abs(BB) > threshold2, abs(BB) <= threshold3)

    Bk_new = np.where(in_thres1, 0, Bk_new)
    Bk_new = np.where(in_thres2, BB - np.sign(BB) * threshold1, Bk_new)
    Bk_new = np.where(in_thres3, (BB - np.sign(BB) * (a / (a - 1)) * threshold1) / (1 - step_size * p / (a - 1)), Bk_new)

    return Bk_new

def MCP_thresholding(BB, p, lam, step_size, a):
    threshold1 = step_size * p * lam
    threshold2 = a * lam

    Bk_new = BB.copy()
    in_thres1 = (abs(BB) <= threshold1)
    in_thres2 = np.logical_and(abs(BB) > threshold1, abs(BB) <= threshold2)

    Bk_new = np.where(in_thres1, 0, Bk_new)
    Bk_new = np.where(in_thres2, (BB - np.sign(BB) * threshold1) / (1 - step_size * p / a), Bk_new)

    return Bk_new


# In[1]:


# def e_step_w(pi_, density, l):
#     density_old = density[l-1]
#     pi_old = pi_[l-1]
#     n = len(density_old[0])
#     k = len(density_old)

#     density_tot = np.hstack(density_old).reshape(k, n).T

#     denom = np.repeat(np.sum(density_tot * pi_old, axis = 1), 2).reshape(n, k)

#     wk = (density_tot * pi_old) / denom

#     # return np.split(wk, k, axis = 1)
#     return [x.ravel() for x in np.split(wk, k, axis = 1)]

def e_step_w(pi_, density, l):
    density_old = density[l-1]
    pi_old = pi_[l-1]
    k = len(density_old)

    log_w = np.column_stack([np.log(pi_old[j]) + safe_log(density_old[j]) for j in range(k)])

    m = np.max(log_w, axis = 1, keepdims = True)

    w_mat = np.exp(log_w - m)
    w_mat /= np.sum(w_mat, axis = 1, keepdims = True)

    w_mat = np.clip(w_mat, 1e-3, 1 - 1e-3)

    w_mat /= np.sum(w_mat, axis = 1, keepdims = True)

    return [w_mat[:, j] for j in range(k)]


# In[8]:


def m_step_pi(w, l):
    w_new = w[l]
    n = len(w_new[0])

    pi_new = (1 / n) * np.sum(w_new, axis = 1)

    return pi_new


# In[78]:


def m_step_Bk_mvFMR(X, Y, Bk, w, l):
    Bk_old = Bk[l-1]
    w_new = w[l]

    n = len(X)
    p = X.shape[1] - 1
    m = Y.shape[1]
    k = len(Bk_old)

    WXX = [np.zeros((p+1, p+1)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            WXX[j] = WXX[j] + w_new[j][i] * np.outer(X[i, :], X[i, :])

    WXY = [np.zeros((p+1, m)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            WXY[j] = WXY[j] + w_new[j][i] * np.outer(X[i, :], Y[i, :])

    Bk_new = [np.linalg.solve(WXX[j], WXY[j]) for j in range(k)]

    return Bk_new


# In[9]:


def PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    inv_sigma_old = inv_sigma[l-1]
    w_new = w[l]
    k = len(pi_new)

    # Initialization #
    Bk_sol = [[]]
    Bk_sol[0] = Bk_old

    # Update Bk #
    Bk_sol.append(
        [
            LASSO_thresholding(Bk_sol[0][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[0][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val)
            for j in range(k)
        ]
    )

    t = 1
    while np.all(diff_norm(Bk_sol[t-1], Bk_sol[t]) >= np.repeat(1e-6, k)):
        Bk_sol.append(
            [
                LASSO_thresholding(Bk_sol[t][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[t][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val)
                for j in range(k)
            ]
        )

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[len(Bk_sol)-1]

    return Bk_new

def m_step_Bk_PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, l):
    Bk_new = PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, l)

    return Bk_new    


# In[10]:


def PGM_mvFMR_SCAD(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    inv_sigma_old = inv_sigma[l-1]
    w_new = w[l]
    k = len(pi_new)

    # Initialization #
    Bk_sol = [[]]
    Bk_sol[0] = Bk_old

    # Update Bk #
    Bk_sol.append(
        [
            SCAD_thresholding(Bk_sol[0][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[0][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val, a)
            for j in range(k)
        ]
    )

    t = 1
    while np.all(diff_norm(Bk_sol[t-1], Bk_sol[t]) >= np.repeat(1e-6, k)):
        Bk_sol.append(
            [
                SCAD_thresholding(Bk_sol[t][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[t][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val, a)
                for j in range(k)
            ]
        )

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[len(Bk_sol)-1]

    return Bk_new

def m_step_Bk_PGM_mvFMR_SCAD(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l):
    Bk_new = PGM_mvFMR_SCAD(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l)

    return Bk_new


# In[11]:


def PGM_mvFMR_MCP(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    inv_sigma_old = inv_sigma[l-1]
    w_new = w[l]
    k = len(pi_new)

    # Initialization #
    Bk_sol = [[]]
    Bk_sol[0] = Bk_old

    # Update Bk #
    Bk_sol.append(
        [
            MCP_thresholding(Bk_sol[0][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[0][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val, a)
            for j in range(k)
        ]
    )

    t = 1
    while np.all(diff_norm(Bk_sol[t-1], Bk_sol[t]) >= np.repeat(1e-6, k)):
        Bk_sol.append(
            [
                MCP_thresholding(Bk_sol[t][j] - eta_val * gradient_obj(X, Y, w_new[j], Bk_sol[t][j], inv_sigma_old[j]), pi_new[j], lamb[j], eta_val, a)
                for j in range(k)
            ]
        )

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[len(Bk_sol)-1]

    return Bk_new

def m_step_Bk_PGM_mvFMR_MCP(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l):
    Bk_new = PGM_mvFMR_MCP(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, a, l)

    return Bk_new


# In[12]:


def primal_residual(B, D):
    return B - D

def dual_residual(U, u, r):
    return (-r) * (U - u)

def epsilon_primal(B, C, e_abs, e_rel):
    e_abs_term = np.sqrt(len(B)) * e_abs
    e_rel_term = max(np.linalg.norm(B, ord = 'fro'), np.linalg.norm(C, ord = 'fro')) * e_rel

    return e_abs_term + e_rel_term

def epsilon_dual(U, e_abs, e_rel):
    e_abs_term = np.sqrt(len(U)) * e_abs
    e_rel_term = np.linalg.norm(U, ord = 'fro') * e_rel

    return e_abs_term + e_rel_term


# In[13]:


def bartels_stewart(A, B, C):
    return linalg.solve_sylvester(A, B, C)


# In[6]:


def ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    sigma_old = sigma[l-1]
    w_new = w[l]

    k = len(pi_new)
    n = X.shape[0]
    p = X.shape[1] - 1
    m = Y.shape[1]

    Bk_sol = [[]]
    Ck = [[]]
    Hk = [[]]
    U1 = [[]]
    epsilon = [np.zeros((k, 2))]

    # R is a list of primal residuals #
    # S is a list of dual residuals #
    R = [[]]
    S = [[]]

    # Initialization #
    Bk_sol[0] = Bk_old.copy()
    Ck[0] = Bk_old.copy()
    Hk[0] = [np.zeros((p+1, m)) for _ in range(k)]
    U1[0] = [np.zeros((p+1, m)) for _ in range(k)]
    R[0] = [np.zeros((p+1, p+1)) for _ in range(k)]
    S[0] = [np.zeros((p+1, p+1)) for _ in range(k)]

    A = [np.zeros((p+1, p+1)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            A[j] = A[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], X[i, :])

    E = [rho * sigma_old[j] for j in range(k)]

    W = [np.zeros((p+1, m)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            W[j] = W[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], Y[i, :])

    G = [W[j] + rho * (Ck[0][j] - U1[0][j]) @ sigma_old[j] for j in range(k)]

    # Update Bk by solving Sylvester's equation #
    Bk_sol.append([])
    Bk_sol[1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

    # Update Hk, Ck #
    Hk.append([])
    Ck.append([])
    Hk[1] = [alpha * Bk_sol[1][j] + (1 - alpha) * Ck[0][j] for j in range(k)]
    Ck[1] = [LASSO_thresholding(Hk[1][j] + U1[0][j], pi_new[j], lamb[j], rho) for j in range(k)]

    # Update U1 #
    U1.append([])
    U1[1] = [U1[0][j] + Hk[1][j] - Ck[1][j] for j in range(k)]

    # Calculating residuals #
    R.append([])
    S.append([])
    R[1] = [primal_residual(Bk_sol[1][j], Ck[1][j]) for j in range(k)]
    S[1] = [dual_residual(U1[1][j], U1[0][j], rho) for j in range(k)]

    epsilon.append([])
    epsilon[1] = np.array([[epsilon_primal(Bk_sol[1][j], Ck[1][j], e_abs, e_rel), epsilon_dual(U1[1][j], e_abs, e_rel)] for j in range(k)])

    t = 1
    while stopping_criteria(R[t], S[t], epsilon[t]):
        G = [W[j] + rho * (Ck[t][j] - U1[t][j]) @ sigma_old[j] for j in range(k)]

        # Update Bk by solving Sylvester's equation #
        Bk_sol.append([])
        Bk_sol[t+1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

        # Update Hk, Ck #
        Hk.append([])
        Ck.append([])
        Hk[t+1] = [alpha * Bk_sol[t+1][j] + (1 - alpha) * Ck[t][j] for j in range(k)]
        Ck[t+1] = [LASSO_thresholding(Hk[t+1][j] + U1[t][j], pi_new[j], lamb[j], rho) for j in range(k)]

        # Update U1 #
        U1.append([])
        U1[t+1] = [U1[t][j] + Hk[t+1][j] - Ck[t+1][j] for j in range(k)]

        # Calculating residuals #
        R.append([])
        S.append([])
        R[t+1] = [primal_residual(Bk_sol[t+1][j], Ck[t+1][j]) for j in range(k)]
        S[t+1] = [dual_residual(U1[t+1][j], U1[t][j], rho) for j in range(k)]

        epsilon.append([])
        epsilon[t+1] = np.array([[epsilon_primal(Bk_sol[t+1][j], Ck[t+1][j], e_abs, e_rel), epsilon_dual(U1[t+1][j], e_abs, e_rel)] for j in range(k)])

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[t]

    return Bk_new

def m_step_Bk_ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, l):
    Bk_new = ADMM_mvFMR_LASSO(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, l)

    return Bk_new


# In[15]:


def ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    sigma_old = sigma[l-1]
    w_new = w[l]

    k = len(pi_new)
    n = len(X)
    p = X.shape[1] - 1
    m = Y.shape[1]

    Bk_sol = [[]]
    Ck = [[]]
    Hk = [[]]
    U1 = [[]]
    epsilon = [np.zeros((k, 2))]

    # R is a list of primal residuals #
    # S is a list of dual residuals #
    R = [[]]
    S = [[]]

    # Initialization #
    Bk_sol[0] = Bk_old.copy()
    Ck[0] = Bk_old.copy()
    Hk[0] = [np.zeros((p+1, m)) for _ in range(k)]
    U1[0] = [np.zeros((p+1, m)) for _ in range(k)]
    R[0] = [np.zeros((p+1, p+1)) for _ in range(k)]
    S[0] = [np.zeros((p+1, p+1)) for _ in range(k)]

    A = [np.zeros((p+1, p+1)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            A[j] = A[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], X[i, :])

    E = [rho * sigma_old[j] for j in range(k)]

    W = [np.zeros((p+1, m)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            W[j] = W[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], Y[i, :])

    G = [W[j] + rho * (Ck[0][j] - U1[0][j]) @ sigma_old[j] for j in range(k)]

    # Update Bk by solving Sylvester's equation #
    Bk_sol.append([])
    Bk_sol[1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

    # Update Hk, Ck #
    Hk.append([])
    Ck.append([])
    Hk[1] = [alpha * Bk_sol[1][j] + (1 - alpha) * Ck[0][j] for j in range(k)]
    Ck[1] = [SCAD_thresholding(Hk[1][j] + U1[0][j], pi_new[j], lamb[j], rho, a) for j in range(k)]

    # Update U1 #
    U1.append([])
    U1[1] = [U1[0][j] + Hk[1][j] - Ck[1][j] for j in range(k)]

    # Calculating residuals #
    R.append([])
    S.append([])
    R[1] = [primal_residual(Bk_sol[1][j], Ck[1][j]) for j in range(k)]
    S[1] = [dual_residual(U1[1][j], U1[0][j], rho) for j in range(k)]

    epsilon.append([])
    epsilon[1] = np.array([[epsilon_primal(Bk_sol[1][j], Ck[1][j], e_abs, e_rel), epsilon_dual(U1[1][j], e_abs, e_rel)] for j in range(k)])

    t = 1
    while stopping_criteria(R[t], S[t], epsilon[t]):
        G = [W[j] + rho * (Ck[t][j] - U1[t][j]) @ sigma_old[j] for j in range(k)]

        # Update Bk by solving Sylvester's equation #
        Bk_sol.append([])
        Bk_sol[t+1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

        # Update Hk, Ck #
        Hk.append([])
        Ck.append([])
        Hk[t+1] = [alpha * Bk_sol[t+1][j] + (1 - alpha) * Ck[t][j] for j in range(k)]
        Ck[t+1] = [SCAD_thresholding(Hk[t+1][j] + U1[t][j], pi_new[j], lamb[j], rho, a) for j in range(k)]

        # Update U1 #
        U1.append([])
        U1[t+1] = [U1[t][j] + Hk[t+1][j] - Ck[t+1][j] for j in range(k)]

        # Calculating residuals #
        R.append([])
        S.append([])
        R[t+1] = [primal_residual(Bk_sol[t+1][j], Ck[t+1][j]) for j in range(k)]
        S[t+1] = [dual_residual(U1[t+1][j], U1[t][j], rho) for j in range(k)]

        epsilon.append([])
        epsilon[t+1] = np.array([[epsilon_primal(Bk_sol[t+1][j], Ck[t+1][j], e_abs, e_rel), epsilon_dual(U1[t+1][j], e_abs, e_rel)] for j in range(k)])

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[len(Bk_sol)-1]

    return Bk_new

def m_step_Bk_ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l):
    Bk_new = ADMM_mvFMR_SCAD(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l)

    return Bk_new


# In[16]:


def ADMM_mvFMR_MCP(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l):
    pi_new = pi_[l]
    Bk_old = Bk[l-1]
    sigma_old = sigma[l-1]
    w_new = w[l]

    k = len(pi_new)
    n = len(X)
    p = X.shape[1] - 1
    m = Y.shape[1]

    Bk_sol = [[]]
    Ck = [[]]
    Hk = [[]]
    U1 = [[]]
    epsilon = [np.zeros((k, 2))]

    # R is a list of primal residuals #
    # S is a list of dual residuals #
    R = [[]]
    S = [[]]

    # Initialization #
    Bk_sol[0] = Bk_old.copy()
    Ck[0] = Bk_old.copy()
    Hk[0] = [np.zeros((p+1, m)) for _ in range(k)]
    U1[0] = [np.zeros((p+1, m)) for _ in range(k)]
    R[0] = [np.zeros((p+1, p+1)) for _ in range(k)]
    S[0] = [np.zeros((p+1, p+1)) for _ in range(k)]

    A = [np.zeros((p+1, p+1)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            A[j] = A[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], X[i, :])

    E = [rho * sigma_old[j] for j in range(k)]

    W = [np.zeros((p+1, m)) for _ in range(k)]
    for j in range(k):
        for i in range(n):
            W[j] = W[j] + (1/n) * w_new[j][i] * np.outer(X[i, :], Y[i, :])

    G = [W[j] + rho * (Ck[0][j] - U1[0][j]) @ sigma_old[j] for j in range(k)]

    # Update Bk by solving Sylvester's equation #
    Bk_sol.append([])
    Bk_sol[1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

    # Update Hk, Ck #
    Hk.append([])
    Ck.append([])
    Hk[1] = [alpha * Bk_sol[1][j] + (1 - alpha) * Ck[0][j] for j in range(k)]
    Ck[1] = [MCP_thresholding(Hk[1][j] + U1[0][j], pi_new[j], lamb[j], rho, a) for j in range(k)]

    # Update U1 #
    U1.append([])
    U1[1] = [U1[0][j] + Hk[1][j] - Ck[1][j] for j in range(k)]

    # Calculating residuals #
    R.append([])
    S.append([])
    R[1] = [primal_residual(Bk_sol[1][j], Ck[1][j]) for j in range(k)]
    S[1] = [dual_residual(U1[1][j], U1[0][j], rho) for j in range(k)]

    epsilon.append([])
    epsilon[1] = np.array([[epsilon_primal(Bk_sol[1][j], Ck[1][j], e_abs, e_rel), epsilon_dual(U1[1][j], e_abs, e_rel)] for j in range(k)])

    t = 1
    while stopping_criteria(R[t], S[t], epsilon[t]):
        G = [W[j] + rho * (Ck[t][j] - U1[t][j]) @ sigma_old[j] for j in range(k)]

        # Update Bk by solving Sylvester's equation #
        Bk_sol.append([])
        Bk_sol[t+1] = [bartels_stewart(A[j], E[j], G[j]) for j in range(k)]

        # Update Hk, Ck #
        Hk.append([])
        Ck.append([])
        Hk[t+1] = [alpha * Bk_sol[t+1][j] + (1 - alpha) * Ck[t][j] for j in range(k)]
        Ck[t+1] = [MCP_thresholding(Hk[t+1][j] + U1[t][j], pi_new[j], lamb[j], rho, a) for j in range(k)]

        # Update U1 #
        U1.append([])
        U1[t+1] = [U1[t][j] + Hk[t+1][j] - Ck[t+1][j] for j in range(k)]

        # Calculating residuals #
        R.append([])
        S.append([])
        R[t+1] = [primal_residual(Bk_sol[t+1][j], Ck[t+1][j]) for j in range(k)]
        S[t+1] = [dual_residual(U1[t+1][j], U1[t][j], rho) for j in range(k)]

        epsilon.append([])
        epsilon[t+1] = np.array([[epsilon_primal(Bk_sol[t+1][j], Ck[t+1][j], e_abs, e_rel), epsilon_dual(U1[t+1][j], e_abs, e_rel)] for j in range(k)])

        if t == 500:
            break
        else:
            t = t + 1

    Bk_new = Bk_sol[len(Bk_sol)-1]

    return Bk_new

def m_step_Bk_ADMM_mvFMR_MCP(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l):
    Bk_new = ADMM_mvFMR_MCP(X, Y, pi_, Bk, sigma, w, e_abs, e_rel, lamb, rho, alpha, a, l)

    return Bk_new


def m_step_Bk_oracle(X, Y, Bk, inv_sigma, true_Bk, w, l):
    Bk_old = Bk[l-1]
    w_new = w[l]
    inv_sigma_old = inv_sigma[l-1]
    
    n = X.shape[0]
    p = X.shape[1] - 1
    q = Y.shape[1]
    k = len(Bk_old)

    Bk_new = []
    for j in range(k):

        w_k = w_new[j]

        XtWX = X.T @ (w_k[:, None] * X) + 1e-6 * np.eye(p + 1)
        XtWY = X.T @ (w_k[:, None] * Y)

        S_inv = inv_sigma_old[j]
        S_inv = 0.5 * (S_inv + S_inv.T)

        mask = (true_Bk[j] != 0)
        free_idx = np.where(mask.ravel(order='F'))[0]

        if len(free_idx) == 0:
            Bk_new.append(np.zeros((p + 1, q)))
            continue

        A = np.kron(S_inv, XtWX)
        rhs = (XtWY @ S_inv).ravel(order='F')

        A_ff = A[np.ix_(free_idx, free_idx)]
        rhs_f = rhs[free_idx]

        try:
            b_f = np.linalg.solve(A_ff, rhs_f)
        except np.linalg.LinAlgError:
            A_ff = A_ff + 1e-8 * np.eye(A_ff.shape[0])
            b_f = np.linalg.solve(A_ff, rhs_f)

        b_vec = np.zeros((p + 1) * q)
        b_vec[free_idx] = b_f

        B_hat = b_vec.reshape((p + 1, q), order='F')

        Bk_new.append(B_hat)

    return Bk_new

    
# In[66]:


def m_step_sigma(X, Y, Bk, w, l):
    n = len(X)
    m = Y.shape[1]

    Bk_new = Bk[l]
    w_new = w[l]
    k = len(Bk_new)

    numer = [np.zeros((m, m)) for _ in range(k)]
    denom = [sum(w_new[j]) for j in range(k)]

    res_mat = [Y - X @ Bk_new[j] for j in range(k)]
    for j in range(k):
        for i in range(n):
            numer[j] = numer[j] + w_new[j][i] * np.outer(res_mat[j][i], res_mat[j][i])

    sigma_new = [numer[j] / denom[j] for j in range(k)]

    return sigma_new


# In[108]:


def total_diff_norm(pi_, Bk, sigma, l):
    pi_new = pi_[l]
    Bk_new = Bk[l]
    sigma_new = sigma[l]

    pi_old = pi_[l-1]
    Bk_old = Bk[l-1]
    sigma_old = sigma[l-1]

    k = len(pi_new)

    pi_norm = np.sum((pi_new - pi_old) ** 2)

    Bk_norm = sum(np.sum((Bk_new[j] - Bk_old[j]) ** 2) for j in range(k))

    sigma_norm = sum(np.sum((sigma_new[j] - sigma_old[j]) ** 2) for j in range(k))

    return np.sqrt(pi_norm + Bk_norm + sigma_norm)


# In[7]:


def modified_BIC(pi_, Bk, density, w, l):
    pi_old = pi_[l - 1]
    w_old = w[l - 1]
    density_old = density[l - 1]
    Bk_old = Bk[l - 1]

    n = len(w_old[0])
    k = len(pi_old)

    log_pi_term = sum(
        np.sum(w_old[j]) * np.log(pi_old[j])
        for j in range(k)
    )

    log_density_term = sum(
        np.sum(w_old[j] * safe_log(density_old[j]))
        for j in range(k)
    )

    negative_ll = -2 * (log_pi_term + log_density_term)

    d_e = k + (k - 1) + sum(np.count_nonzero(b) for b in Bk_old)

    return negative_ll + np.log(n) * d_e

def modified_BIC_6(pi_, Bk, density, w, l):

    def floor_6(x):
        return np.trunc(x * 1e5) / 1e5

    pi_old = pi_[l - 1]
    w_old = w[l - 1]
    density_old = density[l - 1]
    Bk_old = Bk[l - 1]

    n = len(w_old[0])
    K = len(pi_old)

    # logpi_term <- sum(sapply(w_old, sum) * log(pi_old))
    logpi_term = np.sum(
        np.array([np.sum(w_old[j]) for j in range(K)]) * np.log(pi_old)
    )

    # logdensity_term
    logdensity_term = 0.0
    for j in range(K):
        logdensity_term += np.sum(
            w_old[j] * safe_log(density_old[j])
        )

    negative_ll = -2 * (logpi_term + logdensity_term)

    d_e = K + (K - 1)

    for j in range(K):
        B_floor = floor_6(Bk_old[j])
        d_e += np.sum(np.sum(B_floor != 0, axis=0))

    return negative_ll + np.log(n) * d_e


# In[ ]:


# def modified_BIC(pi_, Bk, density, w, l):
#     pi_old = pi_[l-1]
#     w_old = w[l-1]
#     density_old = density[l-1]
#     Bk_old = Bk[l-1]

#     n = len(w_old[0])
#     k = len(pi_old)

#     log_pi_term = [sum(w_old[j]) * np.log(pi_old[j]) for j in range(k)]
#     log_density_term = [sum(w_old[j] * np.log(density_old[j])) for j in range(k)]

#     negative_ll = (-2) * (sum(log_pi_term) + sum(log_density_term))
#     d_e = k + (k - 1) + np.count_nonzero(Bk_old)

#     return negative_ll + (np.log(n) * d_e)

# def modified_BIC_6(pi_, Bk, density, w, l):
#     def floor_6(x):
#         return math.trunc(x * 1e5) / 1e5
#     pi_old = pi_[l-1]
#     w_old = w[l-1]
#     density_old = density[l-1]
#     Bk_old = Bk[l-1]

#     n = len(w_old[0])
#     k = len(pi_old)

#     log_pi_term = [sum(w_old[j]) * np.log(pi_old[j]) for j in range(k)]
#     log_density_term = [np.sum(w_old[j] * np.log(density_old[j])) for j in range(k)]

#     negative_ll = (-2) * (sum(log_pi_term) + sum(log_density_term))

#     floor6_Bk = [b.copy() for b in Bk_old]
#     for j in range(k):
#         for i in range(len(Bk_old[0])):
#             for l in range(Bk_old[0].shape[1]):
#                 floor6_Bk[j][i, l] = floor_6(Bk_old[j][i, l])

#     d_e = k + (k - 1) + np.count_nonzero(floor6_Bk)

#     return negative_ll + (np.log(n) * d_e)

