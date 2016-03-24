# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:02:08 2016

@author: Pedro Leal
"""
import math
import matplotlib.pyplot as plt
from static_model import actuator

#linear =  {'x+': 0.280875, 'y+': -0.969125, 'y-': -0.969125, 'x-': 0.030875}

#sma =  {'x+': 2*math.sqrt(2), 'y+': 0, 'x-': -2, 'y-': 0}
J = {'x':0.5, 'y':0}
sma = {'x+': 0.030875 + J['x'], 'y+': 0.000654043478857, 'x-': -0.219125  + J['x'], 'y-': -0.0131357515566}


#y = {'u':1., 'l':-1}
y = {'u':0.0220568744237, 'l':-0.0220568744237}
#Linear actuator (s)
#l = actuator(linear, J, zero_stress_length = 0.06)
#Sma actuator (l)
s = actuator(sma, J, eps_0 = 0.1)
    
s.find_limits(y)
s.theta = s.max_theta
s.update()
s.plot_actuator()
plt.scatter([J['x']],[y['l']])

