# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:33:31 2016

@author: Pedro Leal
"""
import numpy as np
from scipy import polyfit, polyval
import matplotlib.pyplot as plt

from xfoil_module import output_reader

# Number of point to ignore at each voltage
N = 20
# Import data from IR camera
filename = "voltage_angle.txt"
Data = output_reader(filename, separator='\t', output=None, rows_to_skip=1,
                     header=['Voltage', 'X', 'Y', 'Z'])
VoltageInput = 2800

deflection = Data['Z']
voltage = VoltageInput * np.array(Data['Voltage'])

delta_t = 0.05 # seconds

time = delta_t * np.array(range(len(Data['Z'])))

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


ax2 = ax1.twinx()
ax2.plot(time, voltage, 'r')
ax2.set_ylabel("Voltage (mV)", color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()