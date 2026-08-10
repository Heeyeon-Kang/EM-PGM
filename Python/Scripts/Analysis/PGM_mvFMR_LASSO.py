#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


# In[3]:


from Data.simulation_seed_number import *
from Data.simulation_data import *
import Scripts.Functions.functions as fn

import importlib
importlib.reload(fn)


# In[16]:


lambda_s = [0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
            0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
            0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]


# In[9]:


def PGM_mvFMR_LASSO_nonfixedK(X, Y, eta_val, lambdas, init=False, maxiter=100):
    n = X.shape[0]
    p = X.shape[1] - 1
    m = Y.shape[1]

    BIC_LASSO_s = []
    total_LASSO = []
    optimal_OUTPUTS_LASSO = []
    optimal_w = []
    optimal_density = []
    total_w = []

    for k in range(1, 7):
        BIC_LASSO = []
        w_s = []
        density_s = []
        total_OUTPUTS_LASSO = []

        for z in range(len(lambdas)):
            try:
                density = [[]]
                w = [[]]
                Bk = [[]]
                sigma = [[]]
                inv_sigma = [[]]
                pi_ = [[]]

                ## Initialization ##
                Bk[0] = [np.zeros((p+1, m)) for _ in range(k)]

                if init:
                    pi_[0] = np.repeat(1/k, k)
                else:
                    if k == 1:
                        pi_[0] = [1]
                    elif k == 2:
                        pi_[0] = [0.45, 0.55]
                    elif k == 3:
                        pi_[0] = [0.27, 0.33, 0.4]
                    elif k == 4:
                        pi_[0] = [0.23, 0.24, 0.26, 0.27]
                    elif k == 5:
                        pi_[0] = [0.16, 0.18, 0.2, 0.22, 0.24]
                    elif k == 6:
                        pi_[0] = [0.14, 0.15, 0.16, 0.17, 0.18, 0.2]

                sigma[0] = [np.eye(m) for _ in range(k)]

                lamb = np.repeat(lambdas[z], k)

                theta_diff = [0]

                w[0] = [np.repeat(1/k, n) for _ in range(k)]

                ## The first E and M steps ##
                # E step #
                inv_sigma[0] = [np.linalg.inv(sigma[0][j]) for j in range(k)]
                density[0] = [fn.density_f(X, Y, Bk[0][j], sigma[0][j], inv_sigma[0][j]) for j in range(k)]
                for j in range(k):
                    density[0][j][density[0][j] == 0] = 1e-300

                w.append([])
                w[1] = fn.e_step_w(pi_, density, l = 1)

                # M step #
                pi_.append([])
                pi_[1] = fn.m_step_pi(w, l = 1)

                Bk.append([])
                new_Bk = fn.m_step_Bk_PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, l = 1)
                Bk[1] = new_Bk

                sigma.append([])
                new_sigma = fn.m_step_sigma(X, Y, Bk, w, l = 1)
                reg_new_sigma = [fn.safe_regularize_sigma(new_sigma[j]) for j in range(k)]
                sigma[1] = reg_new_sigma

                theta_diff.append([])
                theta_diff[1] = fn.total_diff_norm(pi_, Bk, sigma, l = 1)

                inv_sigma.append([])
                inv_sigma[1] = [np.linalg.inv(sigma[1][j]) for j in range(k)]
                density.append([])
                density[1] = [fn.density_f(X, Y, Bk[1][j], sigma[1][j], inv_sigma[1][j]) for j in range(k)]
                for j in range(k):
                    density[1][j][density[1][j] == 0] = 1e-300

                ## The iteration of E and M steps ##
                t = 1
                while theta_diff[t] >= 1e-6:
                    # E step #
                    w.append([])
                    w[t+1] = fn.e_step_w(pi_, density, l = t+1)

                    # M step #
                    pi_.append([])
                    pi_[t+1] = fn.m_step_pi(w, l = t+1)

                    Bk.append([])
                    new_Bk = fn.m_step_Bk_PGM_mvFMR_LASSO(X, Y, pi_, Bk, inv_sigma, w, lamb, eta_val, l = t+1)
                    Bk[t+1] = new_Bk

                    sigma.append([])
                    new_sigma = fn.m_step_sigma(X, Y, Bk, w, l = t+1)
                    reg_new_sigma = [fn.safe_regularize_sigma(new_sigma[j]) for j in range(k)]
                    sigma[t+1] = reg_new_sigma

                    theta_diff.append([])
                    theta_diff[t+1] = fn.total_diff_norm(pi_, Bk, sigma, l = t+1)

                    inv_sigma.append([])
                    inv_sigma[t+1] = [np.linalg.inv(sigma[t+1][j]) for j in range(k)]
                    density.append([])
                    density[t+1] = [fn.density_f(X, Y, Bk[t+1][j], sigma[t+1][j], inv_sigma[t+1][j]) for j in range(k)]
                    for j in range(k):
                        density[t+1][j][density[t+1][j] == 0] = 1e-300

                    output = {
                        "lambda": lambdas[z],
                        "pi": pi_[t+1],
                        "Bk": Bk[t+1],
                        "sigma": sigma[t+1]
                    }

                    # print(output)
                    # print(f"Iteration : {t + 1}")
                    # print(f"Lambda : {lambdas[z]}")
                    # print(f"k : {k}")

                    if theta_diff[t+1] == np.inf:
                        break
                    else:
                        if t == maxiter:
                            break
                        else:
                            t = t+1
            except Exception as e:
                    print(e)

            BIC_LASSO.append(fn.modified_BIC(pi_, Bk, density, w, t+1))

            w_s.append(w[t])
            density_s.append(density[t])
            total_OUTPUTS_LASSO.append(output)

        BIC_LASSO_s.append(np.where(np.isinf(BIC_LASSO), np.nan, BIC_LASSO))
        total_w.append(w_s)
        min_LASSO = np.nanargmin(BIC_LASSO_s[k-1])
        total_LASSO.append(total_OUTPUTS_LASSO)
        optimal_OUTPUTS_LASSO.append(total_OUTPUTS_LASSO[min_LASSO])
        optimal_w.append(w_s[min_LASSO])
        optimal_density.append(density_s[min_LASSO])

    optimal_BIC_LASSO_s = [np.min(BIC_LASSO_s[j]) for j in range(len(BIC_LASSO_s))]
    optimal_k = np.nanargmin(optimal_BIC_LASSO_s)
    optimal_OUTPUTS_LASSO_k = optimal_OUTPUTS_LASSO[optimal_k]
    optimal_w_k = optimal_w[optimal_k]
    optimal_density_k = optimal_density[optimal_k]

    return {
        "total": total_LASSO,
        "total_w": optimal_w,
        "total_density": optimal_density,
        "k": (optimal_k + 1),
        "optimal_w": optimal_w_k,
        "BIC": BIC_LASSO_s,
        "optimal": optimal_OUTPUTS_LASSO_k,
        "optimal_density": optimal_density_k,
        "optimal_BIC": optimal_BIC_LASSO_s
    }

