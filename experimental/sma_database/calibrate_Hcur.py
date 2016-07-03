# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:23:31 2016

@author: Pedro Leal
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from xfoil_module import output_reader

def H_cur(sigma, H_min, delta_H, k):
    return H_min + delta_H*(1-np.exp(k*sigma))
    
stress= [50, 100, 150, 172, 200]

eps_t = []
plt.figure()
for stress_i in stress:
    stress_i = str(stress_i) + "MPa"
    raw_data = output_reader("filtered_data_" + stress_i + ".txt", separator=" ", 
                         header = ["Temperature", "Strain", "Stress"],)
    eps_t_i = max(raw_data["Strain"]) - min(raw_data["Strain"])
    eps_t.append(eps_t_i)
    eps = np.array(raw_data["Strain"]) - min(raw_data["Strain"])
    temperature = np.array(raw_data["Temperature"])
    plt.plot(temperature, eps, label = stress_i)

plt.legend()
print eps_t
popt = curve_fit(H_cur, stress, eps_t)

plt.figure()
plt.scatter(stress,eps_t)

