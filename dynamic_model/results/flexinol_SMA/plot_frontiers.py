# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 17:30:32 2016

@author: Pedro Leal
"""

import math
import pickle
import matplotlib.pyplot as plt

frontier_B = pickle.load( open( "B_front.p", "rb" ) )

# unziping data
deflection_B, power_B = zip(*frontier_B)

# converting to degrees
deflection_B = map(lambda x: math.degrees(x), deflection_B )

plt.plot(deflection_B, power_B, lw = 2, c = 'b', label = 'C')
plt.xlabel(r"Angle deflection ($^{\circ}$)", fontsize = 14)
plt.ylabel("Power (J)", fontsize = 14)
plt.legend()
plt.grid()