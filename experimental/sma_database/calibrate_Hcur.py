# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:23:31 2016

@author: Pedro Leal
"""
import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit

import matplotlib.pyplot as plt

from xfoil_module import output_reader

def H_cur(params, sigma, eps_t):
    H_min= params['H_min'].value
    delta_H = params['delta_H'].value
    k = params['k'].value
#    sigma_crit = params['sigma_crit'].value
    
#    model = H_min + delta_H*(1-np.exp(k*(sigma-sigma_crit)))
#    for i in range(len(model)):
#        if sigma[i] < sigma_crit:
#            model[i] = H_min
    model = H_min + delta_H*(1-np.exp(k*(sigma)))
    
    return model - eps_t
    
stress= np.array([50, 100, 150, 172, 200])
stress_fit = np.linspace(0,200,200)
# Martensite Young stiffness (MPa)
E_M = 38.03e3
# Austenite Young stiffness (MPa)
E_A = 18.46e3
# Martensite thermal expansion (1/C)
alpha_M = 0.#6.6e-6
# Austenite thermal expansion  (1/C)
alpha_A = 0.#11e-6
# Initial tem (T0)
Ti = 30.
# Final tem (Tf)
Tf = 140.
eps_t = []
delta_eps_list = []
plt.figure()
for stress_i in stress:
    stress_i_string = str(stress_i) + "MPa"
    raw_data = output_reader("filtered_data_" + stress_i_string + ".txt", separator=" ", 
                         header = ["Temperature", "Strain", "Stress"],)
    delta_eps= max(raw_data["Strain"]) - min(raw_data["Strain"])
    delta_eps_list.append(delta_eps)
    
    eps_t_i = delta_eps + alpha_A*(Tf - Ti) + stress_i*(E_M - E_A)/(E_M*E_A)
    
    print delta_eps,  alpha_A*(Tf - Ti),  stress_i*(E_M - E_A)/(E_M*E_A), eps_t_i
    eps_t.append(eps_t_i)
    eps = np.array(raw_data["Strain"]) - min(raw_data["Strain"])
    temperature = np.array(raw_data["Temperature"])
    plt.plot(temperature, eps, label = stress_i_string)

plt.legend()
plt.grid()
plt.xlabel("Temperature (C)")
plt.ylabel(r"$\Delta \epsilon$ (mm/mm)")
print eps_t

# create a set of Parameters
# 'value' is the initial condition
# 'min' and 'max' define your boundaries
params = Parameters()
params.add('H_min', value= min(eps_t), min=0., max=0.2)
params.add('delta_H', value= max(eps_t), min=0., max=.2)
params.add('k', value= 0, min=-0.01, max=0.)
#params.add('sigma_crit', value= 0, min=0., max=50.)

# do fit, here with leastsq model
result = minimize(H_cur, params, args=(stress, eps_t), method = "differential_evolution")

H_min= result.params['H_min'].value
delta_H = result.params['delta_H'].value
k = result.params['k'].value
#sigma_crit = result.params['sigma_crit'].value

print "H_min: ", H_min
print "H_max: ", H_min  + delta_H
print "k: ", k
#print "sigma_crit: ", sigma_crit

Hcur_fit = H_min + delta_H*(1-np.exp(k*stress_fit))
#Hcur_fit = H_min + delta_H*(1-np.exp(k*(stress_fit-sigma_crit)))
#for i in range(len(Hcur_fit)):
#    if stress_fit[i] < sigma_crit:
#        Hcur_fit[i] = H_min
            
plt.figure()
plt.scatter(stress, eps_t, label = "Experimental data")
plt.plot(stress_fit, Hcur_fit,"--", color="0.75", lw=2, label = "Exponential fit")
plt.xlabel("Stress (MPa)")
plt.ylabel(r"Max. Current Trans. Strain $H^{cur}$ (mm/mm)")
plt.grid()
plt.legend()

Hcur = H_min + delta_H*(1-np.exp(k*stress))
print "New values for eps_M - eps_A: ", Hcur - alpha_A*(Tf - Ti) - stress*(E_M - E_A)/(E_M*E_A)
print "Original values for eps_M - eps_A: ", delta_eps_list

plt.figure()
plt.scatter(stress, delta_eps_list, label = "Experimental data")
plt.plot(stress,  Hcur - alpha_A*(Tf - Ti) - stress*(E_M - E_A)/(E_M*E_A),"--", color="0.75", lw=2, label = "Exponential fit")
plt.xlabel("Stress (MPa)")
plt.ylabel(r"$\Delta \epsilon$ (mm/mm)")