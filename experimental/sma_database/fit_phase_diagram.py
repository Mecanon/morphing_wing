# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 19:56:07 2016

@author: Pedro Leal
"""
import matplotlib.pyplot as plt
import numpy as np

T_Ms_list = np.array([49.9, 73.55, 81.90, 83.6, 104.77])
T_Mf_list = np.array([31.6, 47.71, 62.73, 78.10, 72.75])
T_As_list = np.array([74.8, 79.04, 81.90, 100., 89.37])
T_Af_list = np.array([79.8, 94.81,109.42, 105.6, 109.66])

sigma_Ms_list = np.array([0, 50, 100, 172, 200])
sigma_Mf_list = sigma_Ms_list
sigma_As_list = sigma_Ms_list
sigma_Af_list = sigma_Ms_list

plt.scatter(T_Ms_list, sigma_Ms_list, c='b')
plt.scatter(T_Mf_list, sigma_Mf_list, c='g')
plt.scatter(T_As_list, sigma_As_list, c='k')
plt.scatter(T_Af_list, sigma_Af_list, c='r')

C_Ms = np.polyfit(T_Ms_list, sigma_Ms_list,1)
C_Mf = np.polyfit(T_Mf_list, sigma_Mf_list,1)
C_As = np.polyfit(T_As_list, sigma_As_list,1)
C_Af = np.polyfit(T_Af_list, sigma_Af_list,1)
print C_Ms, C_Mf
p_Ms = np.poly1d(C_Ms)
p_Mf = np.poly1d(C_Mf)
p_As = np.poly1d(C_As)
p_Af = np.poly1d(C_Af)

C_M = (C_Ms[0] + C_Mf[0])/2.
C_A = (C_As[0] + C_Af[0])/2.

print "C_M: ", C_M
print "C_A: ", C_A
print "T_Ms: ", -C_Ms[1]/C_M
print "T_Mf: ", -C_Mf[1]/C_M
print "T_As: ", -C_As[1]/C_A
print "T_Af: ", -C_Af[1]/C_A

T_Ms_fit = (sigma_Ms_list - C_Ms[1])/C_M
T_Mf_fit = (sigma_Mf_list - C_Mf[1])/C_M
T_As_fit = (sigma_As_list - C_As[1])/C_A
T_Af_fit = (sigma_Af_list - C_Af[1])/C_A

plt.plot(T_Ms_fit, sigma_Ms_list, c='g')
plt.plot(T_Mf_fit, sigma_Mf_list, c='g')
plt.plot(T_As_fit, sigma_As_list, c='r')
plt.plot(T_Af_fit, sigma_Af_list, c='r')