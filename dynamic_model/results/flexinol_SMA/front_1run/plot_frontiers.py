# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 17:30:32 2016

@author: Pedro Leal
"""

import math
import pickle
import matplotlib.pyplot as plt

frontier_A = pickle.load( open( "A_front.p", "rb" ) )
frontier_B = pickle.load( open( "B_front.p", "rb" ) )
frontier_C = pickle.load( open( "C_front.p", "rb" ) )

# unziping data
deflection_A, power_A = zip(*frontier_A)
deflection_B, power_B = zip(*frontier_B)
deflection_C, power_C = zip(*frontier_C)

# converting to degrees
deflection_A = map(lambda x: math.degrees(x), deflection_A )
deflection_B = map(lambda x: math.degrees(x), deflection_B )
deflection_C = map(lambda x: math.degrees(x), deflection_C )

plt.plot(deflection_A, power_A, lw = 2, c = 'r', label = 'Design A: free')
plt.plot(deflection_C, power_C, lw = 2, c = 'g', label = 'Design B: colinear')
plt.plot(deflection_B, power_B, lw = 2, c = 'b', label = 'Design C: pulley')

plt.xlabel(r"Angle deflection ($^{\circ}$)", fontsize = 14)
plt.ylabel("Power (J)", fontsize = 14)
plt.legend()
plt.grid()