# -*- coding: utf-8 -*-
"""
@author: Kane
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv('Calcium_Data.csv') # this was the name of the file that contained all of the calcium fits
df = df[(df['Ca'] <= 1000) & (df['Width'] < 500) & (df['Length'] < 650)] # filtering out crysals that were larger than the .png.

# Create dataset index for each Vrf
vrf_values = np.sort(df['Vrf'].unique())
vrf_map = {vrf: i for i, vrf in enumerate(vrf_values)}
df["dataset"] = df["Vrf"].map(vrf_map)

# Extract arrays
N = df["Ca"].values
Vend = (df["Vend"] * df["eta"]).values
Vrf = df["Vrf"].values
z = df["Length"].values
dataset = df["dataset"].values

n_sets = len(vrf_values)

# -------------------------------
# MODEL
# L = [A Vrf^b (eta*Vend)^p + B_i] N^a
# -------------------------------
def global_model(X, *params):

    N, Vend, Vrf, dataset = X
    
    A = params[0]
    b = params[1]
    p = params[2]
    a = params[3]
    
    B = params[4:4+n_sets]
    
    result = np.zeros_like(N)
    
    for i in range(n_sets):
        mask = dataset == i
        result[mask] = (A * Vrf[mask]**b * Vend[mask]**p + B[i]) * N[mask]**a
        
    return result


# Initial guesses
A0 = 180
b0 = -0.11
p0 = -0.33
a0 = 0.33
B0 = [1.0]*n_sets

initial_guess = [A0, b0, p0, a0] + B0

# Fit
params, covariance = curve_fit(
    global_model,
    (N, Vend, Vrf, dataset),
    z,
    p0=initial_guess,
    maxfev=20000
)

param_errors = np.sqrt(np.diag(covariance))

# Extract parameters
A_fit, b_fit, p_fit, a_fit = params[:4]
B_fit = params[4:4+n_sets]

A_err, b_err, p_err, a_err = param_errors[:4]
B_fit_error = param_errors[4:4+n_sets]

# Print results
print("\nGLOBAL PARAMETERS")
print(f"A = {A_fit} ± {A_err}")
print(f"b = {b_fit} ± {b_err}")
print(f"p = {p_fit} ± {p_err}")
print(f"a = {a_fit} ± {a_err}")

print("\nDATASET PARAMETERS (B values)")
for i, vrf in enumerate(vrf_values):
    print(f"Vrf={vrf} | B={B_fit[i]} ± {B_fit_error[i]}")

# fitting parameter B as it varies with Vrf
x = vrf_values
x_err = np.zeros(len(x))
y = B_fit
y_err = B_fit_error

def pwr(x , c1 , m):
    return c1*(x)**m

# power fitting --------------------------------------------------------------
ip_pwr = [1 , 1]
popt , pcov = curve_fit(pwr , x , y , sigma = y_err , absolute_sigma=True, p0 = ip_pwr , maxfev = 1000000)

x_plot = np.linspace(vrf_values[0] , vrf_values[-1] , 100)
y_plot = pwr(x_plot , popt[0] , popt[1])
# power fitting --------------------------------------------------------------

# R^2 calcuation
y_true = B_fit
y_pred_power = pwr(vrf_values, *popt)
ss_res_power = np.sum((y_true - y_pred_power)**2)
ss_tot = np.sum((y_true - np.mean(y_true))**2)

r2_power = 1 - (ss_res_power / ss_tot)

print(f"R^2 for power law fit: {r2_power:.6f}")

print(popt)