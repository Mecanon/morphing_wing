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
#dimensions obtained from the experimental model, in centimeters
c = 51.5
# Defining the positions
# Front attachment point
F = {'x':c*0, 'y':0.}
# Midle attachment point
M = {'x':c*0.214, 'y':0}
# Tip point
T = {'x':c*1, 'y':0}

# Calculating the lengths
l_FT = distance(F,T)
l_FM = distance(F,M)
l_MT = distance(T,M)

plotter(M, T)

#angular variables lists

list_alpha=[]
list_gamma = []
list_beta = []

#bearing radius
R = 1.8
#-----------------------------------------------------------------------------#
#Max strains
max_strainFM = -(0.01*c)/l_FM
max_strainMT = -(0.02*c)/l_MT

list_strainsFM = np.linspace(0, max_strainFM, 2)
list_strainsMT = np.linspace(0, max_strainMT, 2)
for strainFM in list_strainsFM:
    alpha = -(strainFM)*l_FM/R
    list_alpha.append(np.degrees(alpha))
    
    for strainMT in list_strainsMT:
        #in all interations the contration of the spring is the same
        #calculating alpha, beta and gama based on SMA spring strain
        
        gamma = -(strainMT)*l_MT/R
        list_gamma.append(np.degrees(gamma))
        
        M_alpha = {'x':M['x']*np.cos(alpha) - M['y']*np.sin(alpha),
           'y':M['x']*np.sin(alpha) + M['y']*np.cos(alpha)}
    
        beta = alpha + gamma
        list_beta.append(np.degrees(beta))
        
        T_alpha = {'x':T['x']*np.cos(alpha) - T['y']*np.sin(alpha),
           'y':T['x']*np.sin(alpha) + T['y']*np.cos(alpha)}
        
        T_beta = {'x':M_alpha['x'] + 
                      (T_alpha['x']-M_alpha['x'])*np.cos(beta) - 
                      (T_alpha['y']-M_alpha['y'])*np.sin(beta),
                  'y':M_alpha['y'] + 
                      (T_alpha['x']-M_alpha['x'])*np.sin(beta) + 
                      (T_alpha['y']-M_alpha['y'])*np.cos(beta)}
           
        plotter(M_alpha, T_beta, fraction_M = strainFM/max_strainFM,
                fraction_T = strainMT/max_strainMT)
        print distance(F,M_alpha), distance(T_beta,M_alpha)
plt.grid()


#DONE?#TODO: change strain to strain. Strain is negative when contracting and
# positive when expanding. Usually strain = - strain/length_o
#TODO: plots of wingtip displacement versus deflection (alpha,gamma)
#TODO: plots of wingtip displacement versus strain (strainFM, strainMT)
#TODO: plots of deflections (relative, alpha, gamma) versus strains



print  list_alpha
print  list_beta
print  list_gamma    
'''   
    # Calculating the lengths
    l_FT=distance(F,T)
    l_FM=distance(F,M)
    l_MT=distance(T,M)
    print 'l_FM:', l_FM
    print "  "

print 'l_FM:', l_FM

'''

