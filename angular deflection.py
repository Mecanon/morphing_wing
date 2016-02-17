# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:01:34 2015

Rigid body movement, can be couple via l_FA
To be fixed:
-This code is only valid if the lines are like in the paper, on the crontrary it
 will not work.
@author: leal & etancredo
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Used Functions
def plotter(color='r'):
    plt.scatter([F['x'],A['x'],J['x']],[F['y'],A['y'],J['y']], color=color)
    plt.plot([F['x'],A['x'],J['x'],F['x']],[F['y'],A['y'],J['y'],F['y']],color)

def distance(C1,C2):
    dx=abs(C1['x']-C2['x'])
    dy=abs(C1['y']-C2['y'])
    l = np.sqrt(dx**2 + dy**2)

    return l

c=1.0
# Defining the positions
# Front attachment point
F={'x':c*-0.1,'y':1.}
# Aft attachment point
A={'x':c*0.3,'y':0}
# Joint attachment point
J={'x':c*0.0,'y':0}

# Calculating the lengths
l_FJ=distance(F,J)
l_FA=distance(F,A)
l_JA=distance(J,A)

plotter()
theta0=math.atan2(A['y']-J['y'],A['x']-J['x'])


print "inicial theta: ", np.degrees(theta0)


print 'incial l_FA:', l_FA


print '  '

for i in range(30):
    alfa = -0.01*(1+i)
    A={'x':A['x']*np.cos(alfa) - A['y']*np.sin(alfa),
       'y':A['x']*np.sin(alfa) + A['y']*np.cos(alfa)}
    plotter(color='b')
    
    
    theta1=theta0
    theta0=math.atan2(A['y']-J['y'],A['x']-J['x'])
    print "theta0:" , np.degrees(theta0)

   

    # Calculating the lengths
    l_FJ=distance(F,J)
    l_FA=distance(F,A)
    l_JA=distance(J,A)
    print 'l_FA:', l_FA
    print "  "

print 'theta0:' , np.degrees(theta0)    
print 'l_FA:', l_FA




