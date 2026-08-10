#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import numpy as np
import pandas as pd


# In[2]:


def matrix_norm(x):
    x = np.asarray(x)

    if x.ndim >= 2:
        return np.linalg.norm(x, ord = "fro")
    else:
        return np.linalg.norm(x)


def find_best_permutation(true_Bk, pred_Bk):
    k = len(true_Bk)
    permutations = list(itertools.permutations(range(k)))

    norm_values = []
    for perm in permutations:
        total_norm = sum(
            matrix_norm(
                np.asarray(true_Bk[j]) -
                np.asarray(pred_Bk[perm[j]])
            )
            for j in range(k)
        )

        norm_values.append(total_norm)

    best_perm = permutations[np.argmin(norm_values)]

    return best_perm


def align_prediction(true, pred):
    true_Bk = true["Bk"]
    pred_Bk = pred["Bk"]

    best_perm = find_best_permutation(true_Bk, pred_Bk)

    aligned_pred = pred.copy()

    aligned_pred["Bk"] = [pred_Bk[i] for i in best_perm]

    if "pi" in pred:
        aligned_pred["pi"] = np.asarray(pred["pi"])[list(best_perm)]

    if "sigma" in pred:
        aligned_pred["sigma"] = [pred["sigma"][i] for i in best_perm]

    return aligned_pred


def threshold_prediction(Bk, threshold = 5e-6):

    return [
        np.where(
            np.abs(np.asarray(B)) < threshold,
            0,
            np.asarray(B)
        )
        for B in Bk
    ]


# In[3]:


def TPR(true, pred):
    true_Bk = true["Bk"]

    aligned_pred = align_prediction(true, pred)
    new_pred_Bk = aligned_pred["Bk"]

    new_pred_Bk = threshold_prediction(new_pred_Bk)

    number_of_nonzeros = [np.count_nonzero(B) for B in true_Bk]

    tpr = []
    for k in range(len(true_Bk)):
        true_nonzero = np.asarray(true_Bk[k]) != 0
        selected = np.asarray(new_pred_Bk[k]) != 0

        numerator = np.count_nonzero(true_nonzero & selected)

        denominator = number_of_nonzeros[k]

        tpr.append(np.divide(numerator, denominator))

    return np.mean(tpr)


# In[4]:


def FPR(true, pred):
    true_Bk = true["Bk"]

    aligned_pred = align_prediction(true, pred)
    new_pred_Bk = aligned_pred["Bk"]

    new_pred_Bk = threshold_prediction(new_pred_Bk)

    number_of_zeros = [np.count_nonzero(B == 0) for B in true_Bk]

    fpr = []
    for k in range(len(true_Bk)):
        true_zero = np.asarray(true_Bk[k]) == 0
        selected = np.asarray(new_pred_Bk[k]) != 0

        numerator = np.count_nonzero(true_zero & selected)

        denominator = number_of_zeros[k]

        fpr.append(np.divide(numerator, denominator))

    return np.mean(fpr)


# In[5]:


def FDR(true, pred):
    true_Bk = true["Bk"]

    aligned_pred = align_prediction(true, pred)
    new_pred_Bk = threshold_prediction(aligned_pred["Bk"])

    fdr = []
    R = [np.count_nonzero(B != 0) for B in new_pred_Bk]

    for k in range(len(true_Bk)):
        if R[k] == 0:
            fdr.append(0)
        else:
            true_zero = np.asarray(true_Bk[k]) == 0
            selected = np.asarray(new_pred_Bk[k]) != 0

            numerator = np.count_nonzero(true_zero & selected)

            fdr.append(numerator / R[k])

    return np.mean(fdr)


# In[6]:


def MSE(true, pred):
    true_pi = np.asarray(true["pi"])
    true_Bk = true["Bk"]
    true_sigma = true["sigma"]

    aligned_pred = align_prediction(true, pred)

    new_pred_pi = np.asarray(aligned_pred["pi"])
    new_pred_Bk = aligned_pred["Bk"]
    new_pred_sigma = aligned_pred["sigma"]

    mse_pi = []
    mse_Bk = []
    mse_sigma = []
    for k in range(len(true_Bk)):
        mse_pi.append((true_pi[k] - new_pred_pi[k]) ** 2)

        Bk_diff = (np.asarray(true_Bk[k]) - np.asarray(new_pred_Bk[k]))

        mse_Bk.append(matrix_norm(Bk_diff) ** 2 / np.asarray(new_pred_Bk[k]).size)

        sigma_diff = (np.asarray(true_sigma[k]) - np.asarray(new_pred_sigma[k]))

        mse_sigma.append(matrix_norm(sigma_diff) ** 2 / np.asarray(new_pred_sigma[k]).size)

    mse = {
        "pi": np.mean(mse_pi),
        "Bk": np.mean(mse_Bk),
        "sigma": np.mean(mse_sigma)
    }

    return np.array([
        mse["pi"],
        mse["Bk"],
        mse["sigma"]
    ])


# In[7]:


def floor_6(x):
    return np.trunc(x * 1e5) / 1e5

def pred_ll_dens(dat, opt):
    opt_B = opt["Bk"]
    opt_sigma = opt["sigma"]

    opt_inv_sigma = [np.linalg.inv(sigma) for sigma in opt_sigma]

    test_x = dat["test_X"]
    test_y = dat["test_Y"]

    dens_list = [density_f(test_x, test_y, opt_B[k], opt_sigma[k], opt_inv_sigma[k]) for k in range(len(opt_B))]

    return dens_list

def predictive_ll(opt, dens, l, n):
    pi_ = np.asarray(opt[l]["pi"])
    dens_l = dens[l]

    D = np.column_stack(dens_l)

    in_log = D @ pi_

    pred_ll = (-2 * (1/n) * np.sum(np.log(in_log)))

    return pred_ll

