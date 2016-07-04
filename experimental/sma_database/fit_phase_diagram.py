# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 19:56:07 2016

@author: Pedro Leal
"""
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import numpy as np

def cost(x):
    Ts = x[0]
    Tf = x[1]
    C = x[2]
    
    sigma_list = np.array([ 50, 100, 150, 172, 200])
    Ts_array = sigma_list/C + Ts
    Tf_array = sigma_list/C + Tf
    
    rmse = np.sqrt(np.sum((Ts_exp-Ts_array)**2)/len(Ts_array)) + \
           np.sqrt(np.sum((Tf_exp-Tf_array)**2)/len(Tf_array)) 
    
    return rmse

#T_Ms_list = np.array([49.9, 73.55, 81.90, 83.6, 104.77])
#T_Mf_list = np.array([31.6, 47.71, 62.73, 78.10, 72.75])
#T_As_list = np.array([74.8, 79.04, 81.90, 100., 89.37])
#T_Af_list = np.array([79.8, 94.81,109.42, 105.6, 109.66])
  
#T_Ms_list = np.array([49.9, 80.48, 88.7, 95.05, 96.78, 104.00])
#T_Mf_list = np.array([31.6, 58.48, 64.3, 71.6, 71.73, 70.04])
#T_As_list = np.array([74.8,74.44, 79.42, 82.97, 84.39, 92.3])
#T_Af_list = np.array([79.8, 95.28, 101.74, 107.33, 109.82, 114.73])

# Results for DYN1 (plastified wire)
#T_Ms_list = np.array([ 80.48, 88.7, 95.05, 96.78, 104.00])
#T_Mf_list = np.array([58.48, 64.3, 71.6, 71.73, 70.04])
#T_As_list = np.array([74.44, 79.42, 82.97, 84.39, 92.3])
#T_Af_list = np.array([95.28, 101.74, 107.33, 109.82, 114.73])
#sigma_Ms_list = np.array([50, 100, 150, 172, 200])

# Results for DYN2 (plastified wire)
T_Ms_list = np.array([ 68.74, 75.71, 82.33, 84.77, 88.27])
T_Mf_list = np.array([ 57.74, 65.39, 71.29, 74.07, 77.88])
T_As_list = np.array([ 78.47, 83.82, 88.81, 91.38, 94.78])
T_Af_list = np.array([ 88.75, 95.02, 102.15, 105.12, 108.85])
sigma_Ms_list = np.array([50, 100, 150, 172, 200])
sigma_Mf_list = sigma_Ms_list
sigma_As_list = sigma_Ms_list
sigma_Af_list = sigma_Ms_list

plt.scatter(T_Ms_list, sigma_Ms_list, c='b', marker = 'o', label = '$T_{M_s}$')
plt.scatter(T_Mf_list, sigma_Mf_list, c='b', marker = '^', label = '$T_{M_f}$')
plt.scatter(T_As_list, sigma_As_list, c='b', marker = 's', label = '$T_{A_s}$')
plt.scatter(T_Af_list, sigma_Af_list, c='b', marker = 'p', label = '$T_{A_f}$')

Ts_exp = T_As_list
Tf_exp = T_Af_list

bounds = [(30.,90.),(50., 90.),(1.,15.)]
result_A = differential_evolution(cost, bounds)
print result_A.x

Ts_exp = T_Ms_list
Tf_exp = T_Mf_list

bounds = [(50.,90.),(0., 85.),(1.,15.)]
result_M = differential_evolution(cost, bounds, popsize = 400, maxiter =100)
print result_M.x


T_As = result_A.x[0]
T_Af = result_A.x[1]
T_Ms = result_M.x[0]
T_Mf = result_M.x[1]

C_M = result_M.x[2]
C_A = result_A.x[2]

print "C_M: ", C_M
print "C_A: ", C_A
print "T_Ms: ", T_Ms
print "T_Mf: ", T_Mf
print "T_As: ", T_As
print "T_Af: ", T_Af

T_Ms_fit = (sigma_Ms_list)/C_M + T_Ms
T_Mf_fit = (sigma_Mf_list)/C_M + T_Mf
T_As_fit = (sigma_As_list)/C_A + T_As
T_Af_fit = (sigma_Af_list)/C_A + T_Af

print "T_Ms at 172MPa: ", T_Ms_fit[3]
print "T_Mf at 172MPa: ", T_Mf_fit[3]
print "T_As at 172MPa: ", T_As_fit[3]
print "T_Af at 172MPa: ", T_Af_fit[3]


plt.plot(T_Ms_fit, sigma_Ms_list, c='g')
plt.plot(T_Mf_fit, sigma_Mf_list, c='g')
plt.plot(T_As_fit, sigma_As_list, c='r')
plt.plot(T_Af_fit, sigma_Af_list, c='r')
plt.legend(loc = 'upper left')
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")