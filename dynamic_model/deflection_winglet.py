# -*- coding: utf-8 -*-
"""
this code calculates de angular deflection of a rigid winglet with 2 articulations points 
@author: etancredo
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Used Functions
def plotter(M, T, fraction_M = 0., fraction_T = 0.):
    """Plot multiactuated winglet/flap
    :param M: dict with coordinates x and y of the middle joint
    :param T: dict with coordinates x and y of the trailing joint
    :param fraction_M: float [0,1] for color of M
    :param fraction_T: float [0,1] for color of T"""
    
    plt.scatter([F['x']], [F['y']], color = 'b')
    
    c = [fraction_M, 0.0, 1 - fraction_M]
    plt.scatter([M['x']], [M['y']], color = c)
    plt.plot([F['x'],M['x']],[F['y'],M['y']], color = c)
    
    c = [fraction_T, 0.0, 1 - fraction_T]
    plt.scatter([T['x']], [T['y']], color = c)
    plt.plot([M['x'], T['x'],],[M['y'], T['y'],], color = c)

def distance(C1,C2):
    dx=abs(C1['x']-C2['x'])
    dy=abs(C1['y']-C2['y'])
    l = np.sqrt(dx**2 + dy**2)

    return l

c = 1.0
# Defining the positions
# Front attachment point
F = {'x':c*0, 'y':0.}
# Midle attachment point
M = {'x':c*0.3, 'y':0}
# Tip point
T = {'x':c*1, 'y':0}

# Calculating the lengths
l_FT = distance(F,T)
l_FM = distance(F,M)
l_TM = distance(T,M)

plotter(M, T)

#bearing radius
R = 0.1

#Max contractions
max_contractionFM = 0.01*c
max_contractionMT = 0.02*c

list_contractionsFM = np.linspace(0, max_contractionFM, 10)
list_contractionsMT = np.linspace(0, max_contractionMT, 10)
for contractionFM in list_contractionsFM:
    for contractionMT in list_contractionsMT:
        #in all interations the contration of the spring is the same
        #calculating alfa, beta and gama based on SMA spring contraction
#        contractionFM = 0.01 * c
#        contractionMT = 0.02 * c
        alfa = (contractionFM)/R
        gamma = (contractionMT)/R
     
        alfa = alfa
        M_alpha = {'x':M['x']*np.cos(alfa) - M['y']*np.sin(alfa),
           'y':M['x']*np.sin(alfa) + M['y']*np.cos(alfa)}
    
        beta = alfa + gamma
        T_beta = {'x':T['x']*np.cos(beta) - T['y']*np.sin(beta),
           'y':T['x']*np.sin(beta) + T['y']*np.cos(beta)}
    
        plotter(M_alpha, T_beta, fraction_M = contractionFM/max_contractionFM,
                fraction_T = contractionMT/max_contractionMT)
#TODO: change contraction to strain. Strain is negative when contracting and
# positive when expanding. Usually strain = - contraction/length_o
#TODO: plots of wingtip displacement versus deflection (alpha,gamma)
#TODO: plots of wingtip displacement versus strain (contractionFM, contractionMT)
#TODO: plots of deflections (relative, alpha, gamma) versus strains
'''   
    # Calculating the lengths
    l_FT=distance(F,T)
    l_FM=distance(F,M)
    l_TM=distance(T,M)
    print 'l_FM:', l_FM
    print "  "

print 'l_FM:', l_FM

'''


