# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 15:30:04 2016

@author: Pedro Leal
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

# Generate data: for N=1e6, the triangulation hogs 1 GB of memory
N = 1000
x, y = 10 * np.random.random((2, N))

x = list(x)
y = list(y)
for i in range(N):
    if y[i] >= x[i]:
        y[i] = None
        x[i] = None
x = filter(None, x)
y = filter(None, y)
x = np.array(x)
y = np.array(y)
rho = np.sin(3*x) + np.cos(7*y)**3

x_lim = [0,1,1]
y_lim = [0,0,1]
rho_lim = [0,0,0]

x = np.append(x, x_lim)
y = np.append(y, y_lim)
rho = np.append(rho, rho_lim)

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 300), np.linspace(y.min(), y.max(), 300)
xi, yi = np.meshgrid(xi, yi)

# Interpolate; there's also method='cubic' for 2-D data such as here
zi = scipy.interpolate.griddata((x, y), rho, (xi, yi), method='cubic')

plt.imshow(zi, vmin=rho.min(), vmax=rho.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.show()
#

x, y = 10 * np.random.random((2, N))

