# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:32:12 2016

@author: Pedro Leal
"""

#!/usr/bin/env python
'''
Solves Langermann Multimodal Problem with Automatic Optimization Refinement.
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time
from scipy.optimize import minimize

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
import matlab.engine


#adding path to static model
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)

from static_model_B import run
# =============================================================================
# Objective function
# =============================================================================
def objfunc(x):
    x_J = .75
    length_steel = 0.05

    #SMA Pre-stress
    sigma_o = 100e6
    
    inputs = {'sma':{'x-': x_J - length_steel - x[0], 'y-':-x[2],
                     'x+': x_J - length_steel, 'y+':-x[2],
                     'pulley_position':'down'},
              'linear':{'x-':x_J - length_steel - x[1], 'y-':x[2],
                        'x+':x_J - length_steel, 'y+':x[2],
                       'actuator_type': 'wire',
                       'pulley_position':'up'},
              'sigma_o':sigma_o, 'R':x[2]}

    DataFile = open('opt_data.txt','a')
    for x_i in x:
        DataFile.write( '\t %.5f' % (x_i) )
    DataFile.close()
            
#    print inputs
    outputs = run(inputs = inputs, parameters = [eng])
    f = outputs['theta']
#    print f
    DataFile = open('opt_data.txt','a')
    DataFile.write( '\t %.5f' % (f) )    
    DataFile.write('\n')
    DataFile.close()
    
#    f = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]
    fail = 0
    g = []
    return f,g,fail
    
#==============================================================================
# Start Matlab engine
#==============================================================================
eng = matlab.engine.start_matlab()
#Go to directory where matlab file is
eng.cd('..')
eng.cd('SMA_temperature_strain_driven')

# =============================================================================
# 
# =============================================================================
chord = 1.
x_hinge = 0.75
safety = 0.005*chord

opt_prob = Optimization('Static model optimization', objfunc)
#l_s
opt_prob.addVar('x1', 'c', lower = 0.1 , upper = 0.5, value = 0.2)
#l_l
opt_prob.addVar('x2', 'c', lower = 0.1, upper = 0.5, value = 0.2)
#R
opt_prob.addVar('x5', 'c',  lower = 0.001, upper = 0.025, value = 0.02)
# #yl_n
# opt_prob.addVar('x6', 'c', lower = -.9, upper = 0.9, value = -4.593744e-001)
# #xl_p
# opt_prob.addVar('x7', 'c', lower = x_hinge + safety, upper = chord - safety, value =  8.187166e-001)
# #yl_p
# opt_prob.addVar('x8', 'c', lower = -.9, upper = 0., value = -5.719241e-001)

opt_prob.addObj('f')
#opt_prob.addCon('g', 'i')
print opt_prob

DataFile = open('opt_data.txt','w')
key_list = ['l_s', 'l_l', 'R']
output_list = ['theta']
for key in key_list + output_list:
    DataFile.write(key + '\t')
DataFile.write('\n')
DataFile.close()

# Global Optimization
nsga2 = NSGA2()
nsga2.setOption('PopSize', 40)
nsga2.setOption('maxGen', 50)
nsga2(opt_prob)
print opt_prob.solution(0)

# Local Optimization Refinement
#result = minimize(objfunc, [6.817445e-001, -5.216475e-001, 9.029895e-001, 
#                            8.726738e-001, 6.958111e-001, -4.593744e-001,
#                            8.187166e-001, -5.719241e-001 ], method='BFGS',
#                            options={'gtol': 1e-6, 'disp': True})
#slsqp = SLSQP()
#slsqp.setOption('MAXIT', 200)
#slsqp(opt_prob)
#print opt_prob.solution(0)