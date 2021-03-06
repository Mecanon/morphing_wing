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
import math
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

from static_model import run
# =============================================================================
# Objective function
# =============================================================================
def objfunc(x):
    inputs = {'sma':{'x-':x[0], 'y-':x[1], 'x+':x[2], 'y+':x[3]},
              'linear':{'x-':x[4], 'y-':x[5], 'x+':x[6], 'y+':x[7]}}

    DataFile = open('opt_data.txt','a')
    for x_i in x:
        DataFile.write( '\t %.5f' % (x_i) )
    DataFile.close()
            
#    print inputs
    outputs = run(inputs = inputs, parameters = [eng])

    k_nadir = 0
    k_utopia = 514.486
    
    normalized_theta = 1 + outputs['theta']/(math.pi/2.)
    normalized_k = (abs(outputs['k']) - k_nadir)/(k_utopia - k_nadir)
    f = 0.9*normalized_theta + 0.1*normalized_k
    
    DataFile = open('opt_data.txt','a')
    DataFile.write( '\t %.5f \t %.5f \t %.5f' % (outputs['theta'], outputs['k'], f) )      
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
safety = 0.05*chord

opt_prob = Optimization('Static model optimization',objfunc)
#xs_n
opt_prob.addVar('x1', 'c', lower = x_hinge/2. , upper = x_hinge - safety, value = 7.407724e-001)
#ys_n
opt_prob.addVar('x2', 'c', lower = -.9, upper = -0., value = -3.680615e-001)
#xs_p
opt_prob.addVar('x3', 'c', lower = x_hinge + safety, upper = chord - safety, value = 9.933211e-001)
#ys_p
opt_prob.addVar('x4', 'c', lower = 0., upper = .9, value = 6.004423e-001)
#xl_n
opt_prob.addVar('x5', 'c',  lower = x_hinge/2., upper = x_hinge - safety, value = 7.290939e-001)
#yl_n
opt_prob.addVar('x6', 'c', lower = -.9, upper = 0.9, value = -7.584186e-001)
#xl_p
opt_prob.addVar('x7', 'c', lower = x_hinge + safety, upper = chord - safety, value =  7.550874e-001)
#yl_p
opt_prob.addVar('x8', 'c', lower = -.9, upper = 0., value = -4.011175e-001)

opt_prob.addObj('f')
#opt_prob.addCon('g', 'i')
print opt_prob

DataFile = open('opt_data.txt','w')
key_list = ['xs-', 'ys-', 'xs+', 'ys+', 'xl-', 'yl-', 'xl+', 'yl+']
output_list = ['theta', 'k', 'f']
for key in  key_list + output_list:
    DataFile.write(key + '\t')
DataFile.write('\n')
DataFile.close()

# Global Optimization
nsga2 = NSGA2()
nsga2.setOption('PopSize', 40)
nsga2.setOption('maxGen', 50)
nsga2(opt_prob)
print opt_prob.solution(0)

## Local Optimization Refinement
#slsqp = SLSQP()
#slsqp.setOption('MAXIT', 200)
#slsqp(opt_prob.solution(0))
#print opt_prob.solution(0)