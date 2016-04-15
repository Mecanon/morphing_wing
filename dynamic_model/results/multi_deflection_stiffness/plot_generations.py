# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 11:34:28 2016

@author: Pedro Leal
"""
import math

from optimization_tools import plot_generations

filename = "opt_data.txt"
n_generation = 50

outputs_plotted = ['f'] 
output_constrained = ['f']
g1 = lambda x: True if x != 0.5 else False
g = [g1]

p1 = lambda x: math.degrees(x)
p = None #[p1]

plot_generations(filename, g = g, p = p,n_generation = n_generation, 
                 outputs_plotted = outputs_plotted, last_best = False,
                 output_constrained = output_constrained,
                 output_labels = ['Cost'], plot_type = 'all and best',
                 units = None)