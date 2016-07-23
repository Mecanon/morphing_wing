# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 11:34:28 2016

@author: Pedro Leal
"""
import math

from optimization_tools import plot_generations

filename = "opt_data.txt"
n_generation = 50

outputs_plotted = ['theta', 'power'] 
output_constrained = ['theta']
#g1 = lambda x: True if x < 0. else False
#g = [g1]

#p1 = lambda x: math.degrees(x)
#p = [p1]

plot_generations(filename,n_generation = n_generation, 
                 outputs_plotted = outputs_plotted, last_best = False,
#                 output_constrained = output_constrained,
                 output_labels = [r'$\theta$'], plot_type = 'all',
                 units = [r'${}^{\circ}$'], label_size = [14,18])