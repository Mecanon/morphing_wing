# -*- coding: utf-8 -*-
"""
this code calculates de angular deflection of a rigid winglet with 2 articulations points 
@author: etancredo
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Used Functions
def plotter(color='r'):
    plt.scatter([F['x'],M['x'],T['x']],[F['y'],M['y'],T['y']], color=color)
    plt.plot([F['x'],M['x'],T['x'],],[F['y'],M['y'],T['y'],],color)

def distance(C1,C2):
    dx=abs(C1['x']-C2['x'])
    dy=abs(C1['y']-C2['y'])
    l = np.sqrt(dx**2 + dy**2)

    return l

c=1.0
# Defining the positions
# Front attachment point
F={'x':c*0,'y':0.}
# Midle attachment point
M={'x':c*0.3,'y':0}
# Tip point
T={'x':c*1,'y':0}

# Calculating the lengths
l_FT=distance(F,T)
l_FM=distance(F,M)
l_TM=distance(T,M)


plotter()

#bearing radius
R = 0.1



for i in range(1):
    #in all interations the contration of the spring is the same
    #calculating alfa, beta and gama based on SMA spring contraction
    contractionFM = 0.01 * c
    contractionMT = 0.02 * c
    alfa = (contractionFM)/R
    gama = (contractionMT)/R
 
    alfa = alfa
    M={'x':M['x']*np.cos(alfa) - M['y']*np.sin(alfa),
       'y':M['x']*np.sin(alfa) + M['y']*np.cos(alfa)}

    beta = alfa + gama
    T={'x':T['x']*np.cos(beta) - T['y']*np.sin(beta),
       'y':T['x']*np.sin(beta) + T['y']*np.cos(beta)}

    plotter(color='b')
    


print np.degrees(beta)




'''   
    # Calculating the lengths
    l_FT=distance(F,T)
    l_FM=distance(F,M)
    l_TM=distance(T,M)
    print 'l_FM:', l_FM
    print "  "

print 'l_FM:', l_FM

'''


