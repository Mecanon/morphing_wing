# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:31:45 2016

@author: Pedro Leal
"""
import numpy as np
#import matlab.engine
import math 

##Start Matlab engine
#eng = matlab.engine.start_matlab()
##Go to directory where matlab file is
#eng.cd('SMA_temperature_strain_driven')

H_max = 0.047
H_min = 0.
## Temperature
T_0 = 200.15
T_final = 200.15
 
#Initial martensitic volume fraction
MVF_init = 1.

#Initial transformation strain

sigma_0 = 300e6
sigma_crit = 140e6
E_M = 60E9
k = 0.021e-6

#Initial transformation strain
eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_0) - sigma_crit)))
#Define initial strain
eps_0 = eps_t_0 + sigma_0/E_M
    
# Linear strain for testing
n = 40
eps_inp = [eps_0, eps_0]

eps = np.linspace(eps_inp[0], eps_inp[1], n)
#Initial total strain is equal to transformation strain
#eps[0:1] = eps_t_0

def run(T_0, T_final, MVF_init, n, eps, eps_t_0, sigma_0, plot = 'True'):
    #Run SMA model
    for i in range(1,n):
        k = i+1
        if k == n:
            data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init, 
                                      eps_t_0, sigma_0, n, plot, nargout=4) 
        else:
            data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init,
                                      eps_t_0, sigma_0, n, 'False', nargout=4)
        if k != n:
            eps[k] = data[0][i][0]/E_M + data[3][i][0]#eps_t_0
    return data

data = run(T_0, T_final, MVF_init, n, eps, eps_t_0, sigma_0, plot = 'False')
#eps_inp = 