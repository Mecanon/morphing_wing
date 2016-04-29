# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:31 2016

@author: endryws
"""


from actuator import actuator


if __name__ == '__main__':
    
#=======================================================================
# Design Variables
#=======================================================================
    radius = 0.2
    
    sma_geo_props = {'x-': -1, 'x+': -0.5, 'y-': 0, 'y+': radius}
    
    linear_geo_props = {'x-' : -1, 'x+': -0.5, 'y+': 0, 'y-': - radius}
    
    J = {'x': 0, 'y': 0}
    
    k = 0.5
    
    deformacao_inicial = 0.2
    
    design = "B"
    
#=======================================================================
# Creating instances of actuator
#=======================================================================
    
    # for the very first time the sma is modeled like a linear spring
    SMA = 'linear'
    sma = actuator(sma_geo_props, J, radius,k, SMA, design)
    
    linear = 'linear'
    linear = actuator(sma_geo_props, J, radius, k, linear, design)
