# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:31 2016

@author: endryws
"""
import math
import matplotlib.pyplot as plt

from actuator import actuator


if __name__ == '__main__':
    
#=======================================================================
# Design Variables
#=======================================================================
    radius = 0.2
    # To the actuator class work to give you a update, you have to give a
    # zero_stress_length or eps_0 != None
    
    sma_geo_props = {'x-': -1, 'x+': -0.6, 'y-': radius, 'y+': radius, 
    'D': 0.1, 'N':10, 'actuator_type':'spring'}
    
    linear_geo_props = {'x-' : -1, 'x+': -0.6, 'y+': -radius, 
    'y-': -radius, 'D': 0.1, 'N':10, 'actuator_type': 'spring'}
    
    J = {'x': 0, 'y': 0}
    
    k = 0.5
    
    deformacao_inicial = 0
    
    design = "B"
    
#=======================================================================
# Creating instances of actuator
#=======================================================================
    
    # for the very first time the sma is modeled like a linear spring
    SMA = 'linear'
    sma = actuator(sma_geo_props, J, R = radius, material = SMA, 
                   design = design, eps_0 = - deformacao_inicial)
    
    linear = 'linear'
    linear = actuator(linear_geo_props, J, R = radius, material = linear, 
                      design = design, eps_0 = deformacao_inicial,
                      actuator_type='spring')

#=======================================================================
# Frist plot 
#=======================================================================
    cir1 = plt.Circle((0,0), radius, alpha =.7, fc='k')
    cir2 = plt.Circle((0,0), (radius/2), alpha =.5, fc='w')
                                                            
    ax = plt.axes(aspect=1) # Create empty axes (aspect=1 it has to do 
                                    # with the scale
    
    ax.add_patch(cir1)                    
    ax.add_patch(cir2)
    
    cor = 'k--'
    sma.plot_actuator(cor = cor)
    linear.plot_actuator(cor = cor)
#=======================================================================
# Receiving the angle rotated the code delivery the new configuration 
#=======================================================================
    theta = math.radians(30)

    sma.update(theta = -theta)
    linear.update(theta = theta)
    sma.plot_actuator(cor = 'r')
    linear.plot_actuator(cor = 'b')
#========================================================================
# Receiving the deformation(épsilon) the code delivery the angle rotated
#========================================================================
    sma.calculate_theta()
    
    sma.plot_actuator(cor = 'k--')
    linear.plot_actuator(cor = 'y--')

