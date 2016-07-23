# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 16:10:33 2016

@author: Pedro Leal
"""

import matplotlib.pyplot as plt

sigma_0 = [100, 200, 300, 400]
H = [15.2, 14.28, 13.78, 13.62]
W = [1.36, 1.22, 1.01, 0.85]
efficiency = [8.9474, 8.5188, 7.3449, 6.02]
deflection = [-38.91, -38.91, -38.30, -38.51]

plt.figure()
fig, ax1 = plt.subplots()
ax1.plot(sigma_0, H, 'b--', label = 'Heating', lw = 2)
ax1.plot(sigma_0, W, 'b', label = "Work", lw = 2)
ax1.set_xlabel(r'$\sigma_o$ (MPa)', fontsize = 18)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel("Energy (J)", color = 'b', fontsize = 14)
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(sigma_0, efficiency, 'r', lw = 2)
ax2.set_ylabel("Efficiency (%)", color='r', fontsize = 14)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()
ax1.legend(loc = 'best')