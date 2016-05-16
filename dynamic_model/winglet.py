# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:31 2016

@author: endryws
"""
import matplotlib.pyplot as plt

from actuator import actuator


if __name__ == '__main__':
    
#=======================================================================
# Design Variables
#=======================================================================
    radius = 0.2
    
    sma_geo_props = {'x-': -1, 'x+': -0.5, 'y-': radius, 'y+': radius, 
    'D': 0.1, 'N':10, 'actuator_type':'spring'}
    
    linear_geo_props = {'x-' : -1, 'x+': -0.5, 'y+': -radius, 
    'y-': -radius, 'D': 0.1, 'N':10, 'actuator_type': 'spring'}
    
    J = {'x': 0, 'y': 0}
    
    k = 0.5
    
    deformacao_inicial = 0.2
    
    design = "B"
    
#=======================================================================
# Creating instances of actuator
#=======================================================================
    
    # for the very first time the sma is modeled like a linear spring
    SMA = 'linear'
    sma = actuator(sma_geo_props, J, R = radius, material = SMA, 
                   design = design)
    
    linear = 'linear'
    linear = actuator(linear_geo_props, J, R = radius, material = linear, 
                      design = design)

#=======================================================================
# Frist plot 
#=======================================================================
# TODO Ver com o Pedro como colocar os dois plots em um único quadro
    sma.plot_actuator()
    #linear.plot_actuator()
#=======================================================================
