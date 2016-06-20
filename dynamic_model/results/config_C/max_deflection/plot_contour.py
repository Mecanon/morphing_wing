# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 15:30:04 2016

@author: Pedro Leal
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

from xfoil_module import output_reader
from airfoil_module import Naca00XX, create_x

filename = 'opt_data.txt'
safety = 0.05
x_J = 0.75
Data = output_reader(filename)

x = Data['xl-'] + Data['xl+'] + Data['xs-'] + Data['xs+']
y = Data['yl-'] + Data['yl+'] + Data['ys-'] + Data['ys+']
theta = Data['theta'] + Data['theta'] + Data['theta'] + Data['theta']

#Need to denormalize y:
for i in range(len(x)):
    y_airfoil = Naca00XX(1., 0.12, [x[i]], return_dict = 'y')
    thickness = y_airfoil['u'] - y_airfoil['l']
    y[i] = y[i]*thickness/2.

x_data = x[::]
y_data = y[::]
#==============================================================================
# Adding point with nothing
#==============================================================================
def add_nothing_points(x_start, x_end, x, y, z, n = 10, m = 10):
    x_safety = []
    y_safety = []
    x_list = list(np.linspace(x_start, x_end, n))
    y_list = Naca00XX(1., 0.12, x_list, return_dict = 'y', for_xfoil = False)
    thicknesses = []
    for i in range(len(x_list)):
        thicknesses.append(y_list['u'][i] - y_list['l'][i])
    
    for i in range(m):
        x_safety = x_safety + x_list
        for j in range(n):
            y_safety.append(((i+1)/float(m))*thicknesses[j] - thicknesses[j]/2.)
    x = x + x_safety
    y = y + y_safety
    z = z + [0]*len(x_safety)
    return x, y, z

#Near joint
x, y, theta = add_nothing_points(x_J - safety, x_J + safety, x, y, theta, m = 20)
#Near trailing edge
x, y, theta = add_nothing_points(1 - safety, 1, x, y, theta)
#Near leading edge
x, y, theta = add_nothing_points(0, x_J/2., x, y, theta, n = 40)
# Airfoil outer mold
x_outer = []
y_outer = []
x_list = list(create_x(1., n = 100, distribution = 'polar'))
y_list = Naca00XX(1., 0.12, x_list, return_dict = 'y', for_xfoil = False)
x = x + x_list + x_list
y = y + y_list['u'] + y_list['l']
theta = theta + 2*[0]*len(y_list['u'])
#space around the outer_mold at static
y_outer = []
x_list = list(np.linspace(x_J/2., x_J-safety, 20))
y_list = Naca00XX(1., 0.12, x_list, return_dict = 'y', for_xfoil = False)
thicknesses = []
for i in range(len(x_list)):
    thicknesses.append(y_list['u'][i] - y_list['l'][i])

for i in range(3):
    for j in range(20):
        y_outer.append(((i+1)/float(3))*0.1*thicknesses[j]/2. + 0.9*thicknesses[j]/2.)
    for j in range(20):
        y_outer.append(-((i+1)/float(3))*0.1*thicknesses[j]/2. - 0.9*thicknesses[j]/2.)
x = x + 6*x_list
y = y + y_outer
theta = theta + [0]*len(y_outer)
print 'x_list', x_list
print 'y_outer', y_outer
print 'thicknesses', thicknesses

#space around the outer_mold at flap
y_outer = []
x_list = list(np.linspace(x_J + safety, 1.-safety, 20))
y_list = Naca00XX(1., 0.12, x_list, return_dict = 'y', for_xfoil = False)
thicknesses = []
for i in range(len(x_list)):
    thicknesses.append(y_list['u'][i] - y_list['l'][i])

for i in range(3):
    for j in range(20):
        y_outer.append(((i+1)/float(3))*0.1*thicknesses[j]/2. + 0.9*thicknesses[j]/2.)
    for j in range(20):
        y_outer.append(-((i+1)/float(3))*0.1*thicknesses[j]/2. - 0.9*thicknesses[j]/2.)
x = x + 6*x_list
y = y + y_outer
theta = theta + [0]*len(y_outer)

print 'x_list', x_list
print 'y_outer', y_outer
print 'thicknesses', thicknesses

for i in range(len(x_data)):
    plt.scatter(x_data[i],y_data[i], c= ((1.-float(i)/len(x_data), 
                                                float(i)/len(x_data),
                                                0, 1.)))
## Generate data: for N=1e6, the triangulation hogs 1 GB of memory
#N = 1000
#x, y = 10 * np.random.random((2, N))
#
#x = list(x)
#y = list(y)
#for i in range(N):
#    if y[i] >= x[i]:
#        y[i] = None
#        x[i] = None
#x = filter(None, x)
#y = filter(None, y)
#x = np.array(x)
#y = np.array(y)
#rho = np.sin(3*x) + np.cos(7*y)**3

#x_lim = [0,1,0.75,1]
#y_lim = [0,0,1]
#rho_lim = [0,0,0]

#x = np.append(x, x_lim)
#y = np.append(y, y_lim)
#rho = np.append(rho, rho_lim)

# Set up a regular grid of interpolation points
#xi, yi = np.linspace(min(x), max(x), 500), np.linspace(min(y), max(y), 500)
#xi, yi = np.meshgrid(xi, yi)
#
## Interpolate; there's also method='cubic' for 2-D data such as here
#zi = scipy.interpolate.griddata((x, y), theta, (xi, yi), method='cubic')
#
#plt.imshow(zi, vmin=min(theta), vmax=max(theta), origin='lower',
#           extent=[min(x), max(x), min(y), max(y)])
#plt.colorbar()
#plt.show()

#==============================================================================
# Density part
#==============================================================================
#H, xedges, yedges = np.histogram2d(y_data, x_data, bins = (50,50), range = [[0,1],[-.5,.5]])
#plt.contour(xedges[:-1], yedges[:-1], H)