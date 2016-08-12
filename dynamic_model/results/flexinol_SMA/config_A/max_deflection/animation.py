# -*- coding: utf-8 -*-
"""
Created on Tue May 24 10:50:05 2016

@author: Pedro Leal
"""
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import xfoil_module as xf
import airfoil_module as af
from actuator import actuator

#plt.rcParams['animation.ffmpeg_path'] = 'C:\ffmpeg\bin'

Data = pickle.load(open( "data.p", "rb" ))

def plotter(i, database, s, l, x, y, J_x):
    theta = database['theta'][i]
    T = database['T'][i]
    
    #plot SMA actuator
    s.theta = theta
    s.update()
    image = plt.plot([s.x_n + s.x_J, s.x_n + s.x_J + s.r_1], 
             [s.y_n + s.y_J, s.y_n + s.y_J + s.r_2], 'r')
    image.append(plt.scatter([s.x_n + s.x_J, s.x_n + s.x_J + s.r_1], 
                [s.y_n + s.y_J, s.y_n + s.y_J + s.r_2], c = 'r'))
    image.append(plt.scatter([s.x_J],[s.y_J], c = 'r'))

    #Plot linear actuator
    l.theta = theta
    l.update()
    image2 = plt.plot([l.x_n + l.x_J, l.x_n + l.x_J + l.r_1], 
             [l.y_n + l.y_J, l.y_n + l.y_J + l.r_2], 'b')
    image += image2
    image.append(plt.scatter([l.x_n + l.x_J, l.x_n + l.x_J + l.r_1], 
                [l.y_n + l.y_J, l.y_n + l.y_J + l.r_2], c = 'b'))
    image.append(plt.scatter([l.x_J],[l.y_J], c = 'g'))
    
    image.append(plt.text(0, -0.16, 'T = %.0f K' % T))
    image.append(plt.text(0, -0.19, r'$\theta$ = %.1f$^{\circ}$' % math.degrees(theta)))
    
    #plot flap
    image = plot_flap(x, y, J_x, theta = - theta, image = image)
    
    return image

def plot_flap(x, y, x_J, theta, image):
    """
    Plot flap with actuators. theta is clockwise positive.
    @Author: Endryws (modified by Pedro Leal)
    
    Created on Fri Mar 18 14:26:54 2016
    """
    x_dict, y_dict = af.separate_upper_lower(x, y)
    
    # Below I create the dictionarys used to pass to the function find_hinge
    upper = {'x': x_dict['upper'], 'y': y_dict['upper']} # x and y upper points
    lower = {'x': x_dict['lower'], 'y': y_dict['lower']} # x and y lower points
    hinge = af.find_hinge(x_J, upper, lower) 
    
    #=======================================================================
    # With the Joint (hinge) point, i can use the find flap function to
    # found the points of the flap in the airfoil.
    #=======================================================================
    
    data = {'x': x, 'y': y}
    static_data, flap_data = af.find_flap(data, hinge)
    R = hinge['y_upper']
    theta_list = np.linspace(3*math.pi/2, math.pi/2, 50)
    x_circle_list = hinge['x'] + R*np.cos(theta_list)
    y_circle_list = hinge['y'] + R*np.sin(theta_list)

    n_upper = len(flap_data['x'])/2
    
    # Ploting the flap in the original position
#    plt.plot(flap_data['x'][:n_upper],flap_data['y'][:n_upper],'k--')
#    plt.plot(flap_data['x'][n_upper:],flap_data['y'][n_upper:],'k--')         
    # Rotate and plot
    upper = {'x': np.concatenate((flap_data['x'][:n_upper], x_circle_list)),
             'y': np.concatenate((flap_data['y'][:n_upper], y_circle_list))}
    lower = {'x':(flap_data['x'][n_upper:]),
             'y':(flap_data['y'][n_upper:])}
            
    rotated_upper, rotated_lower = af.rotate(upper, lower, hinge, theta, 
                                             unit_theta = 'rad')
    image2 = plt.plot(static_data['x'], static_data['y'],'k')
    
    image3 = plt.plot(rotated_upper['x'], rotated_upper['y'],'k')
    image4 = plt.plot(rotated_lower['x'], rotated_lower['y'],'k')

    image += image2 + image3 + image4
    return image

def thickness(x, t, chord):
    y = af.Naca00XX(chord, t, [x], return_dict = 'y')
    thickness_at_x = y['u'] - y['l']
    return thickness_at_x 

#Positions
J = {'x':0.75, 'y':0.}

sma = {'x-': 4.389066e-001, 'y-': -8.311361e-001, 
       'x+': 7.990382e-001, 'y+': 6.039162e-002}
linear = {'x-': 7.323110e-001, 'y-': 7.573718e-001, 
       'x+': 8.543053e-001, 'y+': -2.499118e-001}
  
airfoil = "naca0012"
chord = 1.#0.6175
t = 0.12*chord

J = {'x':0.75, 'y':0.}

# need to transform normalized coordiantes in to global coordinates
sma['y+'] = sma['y+']*thickness(sma['x+'], t, chord)/2.
sma['y-'] = sma['y-']*thickness(sma['x-'], t, chord)/2.

linear['y+'] = linear['y+']*thickness(linear['x+'], t, chord)/2.
linear['y-'] =  linear['y-']*thickness(linear['x-'], t, chord)/2.

#Generate x,y coordinates for flap plots
filename = xf.file_name(airfoil, output='Coordinates')
Data_xy = xf.output_reader(filename, output='Coordinates',
                        header = ['x','y'])
x = Data_xy['x']
y = Data_xy['y']

#Create actuators
eps_0 = Data['eps_s'][0]
s = actuator(sma, J, eps_0 = eps_0, material = 'SMA')
l = actuator(linear, J, zero_stress_length = 0.2023, material = 'linear')

fig = plt.figure()
fig.gca(xlim=[0,1], ylim=[-0.2,0.1])
fig.gca().set_aspect('equal', adjustable='box')
plt.grid()
plt.xlabel('${}_{I}x + x_J$')
plt.ylabel('${}_{I}y$')
        
ims = []
for i in range(len(Data['T'])): 
    ims.append(plotter(i, Data, s, l, x, y, J['x']))

im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
                                   blit=True)
print "Compiling video"
im_ani.save('im.mp4',dpi=400)
plt.show()
