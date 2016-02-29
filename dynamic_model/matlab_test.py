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

## Temperature
T_0 =273.15 -55
T_final = 273.15 + 55

#Initial martensitic volume fraction
MVF_init = 1.

# Linear strain for testing
n = 40
eps_inp = [0, .0]

eps = np.linspace(eps_inp[0], eps_inp[1], n)

#Run SMA model
for i in range(1,n):
    k = i+1
    if k == n-1:
        data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init, 
                                  'True', nargout=3) 
    else:
        data = eng.OneD_SMA_Model(k, eps[i], T_0, T_final, MVF_init,
                                  'False', nargout=3)

#eps_inp = 