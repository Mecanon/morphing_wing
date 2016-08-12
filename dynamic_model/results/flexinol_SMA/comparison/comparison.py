# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 19:37:45 2016

@author: Pedro Leal
"""

import math
import pickle
import matplotlib.pyplot as plt

data_A = pickle.load( open( "A_data.p", "rb" ) )
data_B = pickle.load( open( "B_data.p", "rb" ) )
data_C = pickle.load( open( "C_data.p", "rb" ) )

# converting to degrees
deflection_A = map(lambda x: math.degrees(x), data_A['theta'] )
deflection_B = map(lambda x: math.degrees(x), data_B['theta'] )
deflection_C = map(lambda x: math.degrees(x), data_C['theta'] )

# converting to force
F_A = map(lambda x: x*math.pi*(0.000381/2.)**2, data_A['sigma'] )
F_B = map(lambda x: x*math.pi*(0.000381/2.)**2, data_B['sigma'] )
F_C = map(lambda x: x*math.pi*(0.000381/2.)**2, data_C['sigma'] )

plt.figure()
plt.plot(data_A['T'], deflection_A, lw = 2, c = 'r', label = "Design A")
plt.plot(data_C['T'], deflection_C, lw = 2, c = 'g', label = "Design B")
plt.plot(data_B['T'], deflection_B, lw = 2, c = 'b', label = "Design C")
plt.grid()
plt.legend(loc = 3)
plt.xlabel("Temperature (K)")
plt.ylabel("Flap deflection (${}^{\circ}$)")

plt.figure()
plt.plot(data_A['T'], F_A, 'r', lw = 2, label = 'Design A SMA')
plt.plot(data_A['T'], data_A['F_l'], 'r--', lw = 2 , label = 'Design A linear')
plt.plot(data_C['T'], F_C, 'g', lw = 2, label = 'Design B SMA' )
plt.plot(data_C['T'], data_C['F_l'], 'g--', lw = 2, label = 'Design B linear' )
plt.plot(data_B['T'], F_B, 'b', lw = 2, label = 'Design C SMA' )
plt.plot(data_B['T'], data_B['F_l'], 'b--', lw = 2, label = 'Design C linear' )
plt.ylabel("Force (N)")
plt.xlabel("Temperature (K)")
plt.grid()
plt.legend(loc = 2)