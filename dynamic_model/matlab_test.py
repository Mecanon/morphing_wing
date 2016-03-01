# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:31:45 2016

@author: Pedro Leal
"""
import numpy as np
import matlab.engine

#Start Matlab engine
eng = matlab.engine.start_matlab()
#Go to directory where matlab file is
eng.cd('SMA_temperature_strain_driven')

H_max = 0.0405
#sigma_crit = 

E_M = 31628.0e6
## Temperature
T_0 = 273.15
T_final = 273.15
 
#Initial martensitic volume fraction
MVF_init = 1.

#Initial transformation strain
eps_t_0 = 0. #H_max
# Linear strain for testing
n = 40
eps_inp = [0., eps_t_0]

eps = np.linspace(eps_inp[0], H_max, n)
#Initial total strain is equal to transformation strain
#eps[0:1] = eps_t_0

def run(T_0, T_final, MVF_init, n, eps, eps_t_0, plot = 'True'):
    #Run SMA model
    for i in range(1,n):
        k = i+1
        if k == n:
            data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init, 
                                      eps_t_0, plot, nargout=4) 
        else:
            data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init,
                                      eps_t_0, 'False', nargout=4)
    return data

#eps_inp = 