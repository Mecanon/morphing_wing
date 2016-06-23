# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 18:42:21 2016

@author: Pedro Leal
"""

import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping
import numpy as np

from xfoil_module import output_reader

def cost(x):
    T_2 = x[0]
    T_3 = x[1]
    strain_2 = x[2]
    strain_3 = x[3]
    
    f = []
    for T_i in T:
        f_i = tangent_lines(T_i, T_1, T_2, T_3, T_4, strain_1, strain_2, strain_3, strain_4)
        f.append(f_i)
        
    f = np.array(f)
    strain_np = np.array(strain)
    
    rmse = np.sqrt(np.sum((f-strain_np)**2)/len(strain))
    
    return rmse
    
def tangent_lines(T, T_1, T_2, T_3, T_4, strain_1, strain_2, strain_3, strain_4):
    
    if T < T_2:
        return (strain_2 - strain_1)/(T_2 - T_1)*(T-T_1) + strain_1
    elif T < T_3:
        return (strain_3 - strain_2)/(T_3 - T_2)*(T-T_2) + strain_2
    else:
        return (strain_4 - strain_3)/(T_4 - T_3)*(T - T_4) + strain_4
        
raw_data = output_reader("filtered_data_100MPa.txt", header = ['Temperature',
                                               "Strain", "Stress"],)

temperature = raw_data['Temperature']
strain = raw_data['Strain']
stress = raw_data['Stress']


for i in range(len(temperature)):
    if stress[i] >=  100.:
        break

temperature = temperature[i:]
strain = strain[i:]

i = temperature.index(max(temperature))

temperature_A =  temperature[:i+1]
temperature_M =  temperature[i:]
strain_A =  strain[:i+1]
strain_M =  strain[i:]

#==============================================================================
# # Austenitic transformation
#==============================================================================
T = temperature_A
strain = strain_A

T_1  = temperature_A[0]
T_4  = temperature_A[-1]
strain_1 = strain_A[0]
strain_4 = strain_A[-1]

T_2 = 60
T_3 = 90
strain_2 = .05
strain_3 = .01


x0 = [T_2, T_3, strain_2, strain_3, ]
#cost(T_2, T_3, strain_1, strain_2, strain_3, strain_4)
#result = minimize(cost, x0, method = 'BFGS')
result = basinhopping(cost, x0, niter=200)
print "Austenite calibration"
print result.x

x = result.x
T_2 = x[0]
T_3 = x[1]
strain_2 = x[2]
strain_3 = x[3]

f = []
for T_i in T:
    f_i = tangent_lines(T_i, T_1, T_2, T_3, T_4, strain_1, strain_2, strain_3, strain_4)
    f.append(f_i)
    
f = np.array(f)
    
plt.plot(T,f)
plt.plot(temperature_A, strain_A)

#==============================================================================
#  Martenitic transformation
#==============================================================================
T = temperature_M
strain = strain_M

T_1  = temperature_M[0]
T_4  = temperature_M[-1]
strain_1 = strain_M[0]
strain_4 = strain_M[-1]

T_2 = 60
T_3 = 90
strain_2 = .05
strain_3 = .01

x0 = [T_2, T_3, strain_2, strain_3]
#cost(T_2, T_3, strain_1, strain_2, strain_3, strain_4)
#result = minimize(cost, x0, method = 'BFGS')
result = basinhopping(cost, x0, niter=200)
print "Martensite calibration"
print result.x

x = result.x
T_2 = x[0]
T_3 = x[1]
strain_2 = x[2]
strain_3 = x[3]

f = []
for T_i in T:
    f_i = tangent_lines(T_i, T_1, T_2, T_3, T_4, strain_1, strain_2, strain_3, strain_4)
    f.append(f_i)
    
f = np.array(f)
plt.plot(T,f)
plt.plot(temperature_M, strain_M)
