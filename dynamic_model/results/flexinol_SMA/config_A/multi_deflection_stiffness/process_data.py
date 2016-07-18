# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:16:36 2016

@author: Pedro Leal
"""

import pickle
import matplotlib.pyplot as plt

pop = pickle.load( open( "pop.p", "rb" ) )
front = pickle.load( open( "front.p", "rb" ) )
#stats = pickle.load( open( "stats.p", "rb" ) )

# Current results are for radius = 0.00025
theta = []
power = []
for i in range(len(front)):
    theta.append(front[i][0])
    power.append((0.000381/0.00025)**2*front[i][1])
plt.scatter(theta, power)