# -*- coding: utf-8 -*-
"""
Created on Wed May 18 14:28:55 2016

@author: Pedro Leal
"""

import matplotlib.pyplot as plt
import math

weights = [0.6, 0.7, 0.8, 0.9, 1.]
theta_list = [-0.14832, -0.30202, -0.34233, -0.39834, -0.76897]
k_list = [346.3982, 386.24, 389.09436, 464.9625, 514.486]
markers = ['o', 's', '1', '2', '*']
#convert to degrees
for i in range(len(theta_list)):
    theta_list[i] = math.degrees(theta_list[i])

plt.figure()
plt.plot(k_list, theta_list, linewidth = 2.)
for i in range(len(theta_list)):
    plt.scatter(k_list[i], theta_list[i],
                label = r'$w_{\theta}$= %.1f, $w_k$= %1.f' % (weights[i],1-weights[i]),
                color = (( 0, float(i)/(len(theta_list)-1), 1 - float(i)/(len(theta_list)-1))),
                linewidth = 4., marker = 's')
plt.xlabel('Stiffness coefficient (N/m)')
plt.ylabel(r'$\theta ({}^{\circ})$')
plt.grid()
plt.legend(loc='lower left')