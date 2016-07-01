# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 19:18:23 2016

@author: Pedro Leal
"""

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

# =============================================================================
# Objective function
# =============================================================================
def objfunc_1(x):
    f = x[0] + x[1]
    fail = 0
    g = []
    return f,g,fail

def objfunc_2(x):
    f = x[0] - x[1]
    fail = 0
    g = []
    return f,g,fail
    
def objfunc_3(x):
    f1 = x[0] - x[1]
    f2 = x[0] + x[1]
    f = (f1, f2)
    fail = 0
    g = []
    return f,(g, g), (fail, fail)
# =============================================================================
# 
# =============================================================================
chord = 1.
x_hinge = 0.75
safety = 0.005*chord

opt_prob = Optimization('main', (objfunc_1, objfunc_2))
opt_prob.addObj("f1")
opt_prob.addObj("f2")
#xs_n
opt_prob.addVar('x1', 'c', lower = -1 , upper = 1, value = 6.817445e-001)
#ys_n
opt_prob.addVar('x2', 'c', lower = -1, upper = 1, value = -5.216475e-001)

#opt_prob.addObj('2', objfunc_2)
print opt_prob

# Global Optimization
nsga2 = NSGA2()
nsga2.setOption('PopSize', 10)
nsga2.setOption('maxGen', 10)
nsga2(opt_prob)
print opt_prob.solution(0)