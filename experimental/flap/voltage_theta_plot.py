# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:33:31 2016

@author: Pedro Leal
"""
import math
import numpy as np
from scipy import polyfit, polyval
import matplotlib.pyplot as plt

from xfoil_module import output_reader

# Wire properties
rho = 3.55041728247e-06
A = math.pi*(0.000381/2.)**2
L = math.sqrt((0.168 + 0.018)**2 +  (0.06 + 0.01)**2)

# Number of point to ignore at each voltage
N = 20
# Import data from IR camera
filename = "voltage_angle_run4.txt"
Data = output_reader(filename, separator='\t', output=None, rows_to_skip=1,
                     header=['Date', 'Time', 'Voltage', 'X',  'Z'], column_types = [str, str, float, 
                     float, float])

deflection = Data['Z']
voltage =  np.array(Data['Voltage'])
delta_t = 0.05 # seconds
time = delta_t * np.array(range(len(Data['Z'])))

#Filter drift
#(drift, drift_0) = polyfit(time[0:7000], deflection[0:7000],1)
#
#deflection = deflection - drift*time

plt.figure()
plt.scatter(voltage, deflection)
plt.xlabel("Voltage (mV)")
plt.ylabel("Flap deflection (${}^{\circ}$)")
plt.grid()

fig, ax1 = plt.subplots()
ax1.plot(time, deflection, 'b-')
ax1.set_xlabel('Time (s)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel("Flap deflection (${}^{\circ}$)", color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

current = A*voltage/(rho*L)

ax2 = ax1.twinx()
ax2.plot(time, voltage, 'r')
ax2.set_ylabel("Voltage (mV)", color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()

plt.figure()
plt.scatter(voltage, deflection)
plt.xlabel("Voltage (mV)")
plt.ylabel("Flap deflection (${}^{\circ}$)")
plt.grid()

# Plot with current
plt.figure()
plt.scatter(current, deflection)
plt.xlabel("Current (mA)")
plt.ylabel("Flap deflection (${}^{\circ}$)")
plt.grid()