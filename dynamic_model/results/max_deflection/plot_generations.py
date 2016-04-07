# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 11:34:28 2016

@author: Pedro Leal
"""

from optimization_tools import plot_generations

filename = "opt_data.txt"
n_generation = 50

outputs_plotted = ['theta'] 
output_constrained = ['theta']
g1 = lambda x: True if x < 0. else False
g = [g1]

plot_generations(filename, g = g, n_generation = n_generation, 
                 outputs_plotted = outputs_plotted, last_best = False,
                 output_constrained = output_constrained,
                 output_labels = [r'$\theta$'], plot_type = 'all and best',
                 units = [r'${}^{\circ}$'])