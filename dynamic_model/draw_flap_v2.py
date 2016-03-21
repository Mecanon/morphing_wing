# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 14:26:54 2016

@author: endryws
"""
#=======================================================================
# Importing librarys
#=======================================================================

import math
import numpy as np          
import pylab as plt

#=======================================================================
# This function calculates the points we need to plot the inicial 
# NACA00XX series.
# This code can be founded in Aeropy repository in Github.  
# To our aplication we will use return_dict = 'xy', to receive two 
# dictionarys x points and y points
# in the dictionarys we have two keys u and l, u to upper points, and
# l to lower points.
#=======================================================================

def Naca00XX(c, t, x_list, return_dict = 'y'):
    """
    Generates a simetric NACA airfoil.
    Inputs:
    :param c: chord
    :param t: max thickness
    :param x_list: points between 0 and c (chord lenght) 
    :param return_dict: return dictionary y('y') or x and y ('xy')
    Returns dictionary with keys 'u'pper and 'l'ower
    
    The Naca function can be found in: 
    https://en.wikipedia.org/wiki/NACA_airfoil  

    Created on Wed Feb 03 12:50:52 2016
    
    @author: Endryws and Pedro Leal
    """
    y_upper = []
    y_lower = []
    for x in x_list:
        xc= x/c # Is just for increase speed and 
                # facilitate future changes.
        a1 = 5*t*c
        t1 = 0.2969*(math.sqrt(xc))
        t2 = -0.1260*(xc)
        t3 = -0.3516*((xc)**2)
        t4 = 0.2843*((xc)**3)
        t5 = -0.1015*((xc)**4)
        y = (a1*(t1+t2+t3+t4+t5))
        y_upper.append(y)
        y_lower.append(y*(-1)) # is just for pick the y axis
                                       # negative numbers
    y = {'u': y_upper[:-1],
         'l':y_lower[-2::-1]}
    if return_dict == 'y':
        return y
        
    elif return_dict == 'xy':
        x = {'u':x_list[:-1],
             'l':x_list[-2::-1]}
        return x, y
#=======================================================================
# This function can be founded in Aeropy.
# This function find the joint point, and returns a dictionary with 
# the points x and y of the joint.
# Can be used only for Naca00XX type of airfoil.
#=======================================================================

def find_hinge(x_hinge, upper, lower):
    """From the points of upper and lower surface find the y coordinate
    of the hinge at x_hinge
    
    :param x_hinge: float x-coordinate of the hinge
    :param upper: dictionary with keys x and y, coordiantes of 
                  upper surface
    :param lower: dictionary with keys x and y, coordiantes of 
                  lower surface    
    """
    def find_closest(data, x_hinge, fwd, aft):
        """ Find the points afterwards(aft) and forwards(fwd) that are 
        closest to x_hinge. Retuns dictionaries aft and fwd, each with 
        keys x and y."""
        for i in range(len(data['x'])):
            xi = data['x'][i]
            yi = data['y'][i]
            error = abs(xi - x_hinge)
    
            #If equal just save and break it
            if x_hinge == xi:
                aft['x'] = xi
                aft['y'] = yi
                aft['error'] = error
                fwd['x'] = xi
                fwd['y'] = yi
                fwd['error'] = error
                break
            
            if xi > x_hinge and error < aft['error']:
                aft['x'] = xi
                aft['y'] = yi
                aft['error'] = error
    
            if xi < x_hinge and error < fwd['error']:
                fwd['x'] = xi
                fwd['y'] = yi
                fwd['error'] = error
        return fwd, aft
        
    #Find y hinge
    upper_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    upper_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    hinge = {'x': x_hinge, 'y': 9e9}
    
    upper_fwd, upper_aft = find_closest(upper, x_hinge, upper_fwd, upper_aft)
    lower_fwd, lower_aft = find_closest(lower, x_hinge, lower_fwd, lower_aft) 
    
    #Interpolate
    hinge['y_upper'] =  upper_fwd['y'] + (
                        hinge['x'] - upper_fwd['x'])*(
                        upper_aft['y'] - upper_fwd['y'])/(
                        upper_aft['x'] - upper_fwd['x'])
    
    hinge['y_lower'] =  lower_fwd['y'] + (
                        hinge['x'] - lower_fwd['x'])*(
                        lower_aft['y'] - lower_fwd['y'])/(
                        lower_aft['x'] - lower_fwd['x'])
    hinge['y']  = (hinge['y_upper'] + hinge['y_lower'])/2.
    
    return hinge
        
#=======================================================================
#This function can be founded in Aeropy code.
# Used to find the flap points.
#=======================================================================

def find_flap(data, hinge, extra_points = None):
    """Create the static airfoil and flap dictionaries containing the 
    outer mold coordinates of both. Because it is necessary to have an
    intersection between the static lower surface and the flap lower 
    surface, it is sometimes interesting to have an extra poin in the 
    static surface to garantee the intersection.
    
    :param data: dictionary with x and y cooridnates of the whole outer
    mold
    
    :param hinge: dictionary with x and y coordinates of hinge. Can be 
    found via find_hinge.
                  
    :param extra_points: include extra points to surface, extending it 
    more than it will actually be (usefull for intersecting lines)"""
    
    #Generate empty dictionaries which will store final data.
    flap_data = {'x':[], 'y':[]}
    static_data = {'x':[], 'y':[]}
    
    type = None
    
    for i in range(len(data['x'])):
        xi = data['x'][i]
        yi = data['y'][i]
    #Because of the way that xfoil works, the upper list will always
    #begin from the trailing edge    
        if xi > hinge['x']:
            if type == None:
                type = 'upper'
            elif type == 'lower':
                flap_data['x'].append(hinge['x'])
                flap_data['y'].append(hinge['y_' + type])
                static_data['x'].append(hinge['x'])
                static_data['y'].append(hinge['y_' + type])
                # this will make the at the hinge to be included 
                # only once
                type = 'upper'
                #If the extra point is unecessary in the lower
                if extra_points == 'lower':
                    static_data['x'].append(data['x'][i+1])
                    static_data['y'].append(data['y'][i+1])

            flap_data['x'].append(xi)
            flap_data['y'].append(yi)
            
    #Because of the way that xfoil works, the lower list will always
    #begin from the  leading edge  
        else:
            if type == None:
                type = 'lower'
            elif type == 'upper':
                #If the extra point is unecessary in the upper
                if extra_points == 'upper' and len(static_data['x']) == 0:
                    static_data['x'].append(data['x'][i-2])
                    static_data['y'].append(data['y'][i-2])
                    static_data['x'].append(data['x'][i-1])
                    static_data['y'].append(data['y'][i-1])
                flap_data['x'].append(hinge['x'])
                flap_data['y'].append(hinge['y_' + type])
                static_data['x'].append(hinge['x'])
                static_data['y'].append(hinge['y_' + type])
                # this will make the at the hinge to be included only 
                # once
                type = 'lower'
                
            static_data['x'].append(xi)
            static_data['y'].append(yi)   
    return static_data, flap_data
    
#=======================================================================
# This function can be founded in Aeropy.
# I understood that this function rotates the points given in theta.
#=======================================================================
    
def rotate(upper, lower, origin, theta, unit_theta = 'deg'):
    """
    :param upper: dictionary with keys x and y, each a list
    
    :param lower: dictionary with heys x and y, each a list
    
    :param origin: dictionary with keys x and y, each a float
    
    :param theta: float representing angle in degrees clock-wise
    """
    output = []
    
    #For trigonometric relations in numpy, theta must be in radians
    if unit_theta == 'deg':
        theta = theta * np.pi/180.
    # Rotation transformation Matrix
    T = [[np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]]
        
    for coordinates in [upper, lower]:
        rotated_coordinates = {'x':[], 'y':[]}
        for i in range(len(coordinates['x'])):
            # The rotation must take place with the center of rotation
            # as the origin
            cx = coordinates['x'][i] - origin['x']
            cy = coordinates['y'][i] - origin['y']
            # Rotate
            rot_x = T[0][0]*cx + T[0][1]*cy
            rot_y = T[1][0]*cx + T[1][1]*cy
            # Store and add back the values of the origin
            rotated_coordinates['x'].append(rot_x + origin['x'])
            rotated_coordinates['y'].append(rot_y + origin['y'])
        output.append(rotated_coordinates)
    return output[0], output[1]
    
#=======================================================================
# Here comes the plot of the initial Naca 00XX airfoil.
#=======================================================================

# Parameters to the Naca00XXfunction
chord = 1
t = 0.12 # max thickness
numberOfPoints = 5000 # number of points between the inicial and final
                      # points of the chord
x_list = np.linspace(0,chord,numberOfPoints)

xPoints, yPoints = Naca00XX(chord,t,x_list, 'xy')


    
plt.plot(xPoints['u'],yPoints['u'],'k')
plt.plot(xPoints['l'],yPoints['l'],'k')

#=======================================================================
# Here I use the find hinge function to found the joint point.
#=======================================================================

Jx = float(chord*(0.65)) # The x point of the joint (hinge).
                 # I use 2/3 of the chord, but thisparameter can be 
                 # changed

# Below I create the dictionarys used to pass to the function find_hinge

upper = {'x': xPoints['u'], 'y': yPoints['u']} # x and y upper points

lower = {'x': xPoints['l'], 'y': yPoints['l']} # x and y lower points

hinge = find_hinge(Jx, upper, lower) 

#=======================================================================
# With the Joint (hinge) point, i can use the find flap function to
# found the points of the flap in the airfoil.
#=======================================================================

data = {'x': np.concatenate((xPoints['u'], xPoints['l'])),
        'y': np.concatenate((yPoints['u'], yPoints['l']))}


static_data, flap_data = find_flap(data, hinge)

R = hinge['y_upper']
 
theta_list = np.linspace(math.pi/2, 3*math.pi/2,1000)

x_circle_list = hinge['x'] + R*np.cos(theta_list)
y_circle_list = hinge['y'] + R*np.sin(theta_list)

# Ploting the flap in the original position

plt.plot( np.concatenate((x_circle_list, x_circle_list,flap_data['x'])),
         np.concatenate((y_circle_list,y_circle_list*(-1),flap_data['y']
         )),'k')
         

#=======================================================================
# Using the rotate function to generate the points rotated.
#=======================================================================

theta = 30 # angle to rotate in degrees
nPositive = len(flap_data['x'])/2

upper = {'x': np.concatenate((x_circle_list, flap_data['x'][:nPositive])),
         'y': np.concatenate((y_circle_list,flap_data['y'][:nPositive]))}
    
lower = {'x':(flap_data['x'][nPositive:]),
         'y':(flap_data['y'][nPositive:])}
        
rotatedPositive,rotatedNegative = rotate(upper, lower, hinge, theta)


plt.plot(rotatedPositive['x'], rotatedPositive['y'],'k')
plt.plot(rotatedNegative['x'], rotatedNegative['y'],'k')

plt.axes().set_aspect('equal')