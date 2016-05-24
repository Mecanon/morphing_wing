# -*- coding: utf-8 -*-
"""
Created on Tue May 24 10:50:05 2016

@author: Pedro Leal
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from actuator import actuator

def plotter(theta, m):
    print theta
    m.theta = theta
    m.update()
    if m.material == 'linear':
        colour = 'b'
    elif m.material == 'SMA':
        colour = 'r'
    image = plt.plot([m.x_n + m.x_J, m.x_n + m.x_J + m.r_1], 
             [m.y_n + m.y_J, m.y_n + m.y_J + m.r_2], colour)
    image.append(plt.scatter([m.x_n + m.x_J, m.x_n + m.x_J + m.r_1], 
                [m.y_n + m.y_J, m.y_n + m.y_J + m.r_2], c=colour))
    image.append(plt.scatter([m.x_J],[m.y_J], c = colour))
    return image
geo_props = {'x-':0, 'y-':0, 'x+':1, 'y+':0}
J = {'x':0., 'y':0.}
mechanism = actuator(geo_props, J, zero_stress_length=1.)

fig = plt.figure()

n = 10
theta_list = np.linspace(0, math.pi/2., n)

ims = []
for theta in theta_list: 
    ims.append(plotter(theta, mechanism))

im_ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=3000,
                                   blit=True)
#im_ani.save('im.mp4', metadata={'artist':'Guido'})

plt.show()
