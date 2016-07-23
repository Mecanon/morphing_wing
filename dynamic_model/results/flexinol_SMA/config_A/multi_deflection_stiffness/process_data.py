# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:16:36 2016

@author: Pedro Leal
"""
import math
import pickle
import matplotlib.pyplot as plt

pop = pickle.load( open( "pop.p", "rb" ) )
front = pickle.load( open( "front.p", "rb" ) )
#stats = pickle.load( open( "stats.p", "rb" ) )

# unziping data
deflection_frontier, power_frontier = zip(*front)
deflection_pop, power_pop = zip(*pop)

# converting to degrees
deflection_frontier = map(lambda x: math.degrees(x), deflection_frontier )
deflection_pop = map(lambda x: math.degrees(x), deflection_pop )

plt.scatter(deflection_pop, power_pop, color = 'r')
plt.scatter(deflection_frontier, power_frontier, color = 'b')