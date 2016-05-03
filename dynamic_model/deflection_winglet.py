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

#angular variables lists

list_gamma=[]
list_beta = []
list_delta = []
aux_gamma = []
aux_beta = []

#wingtip height list

list_wingtipheightgamma = []
list_wingtipheightbeta = []


#bearing radius
R = 1.8
#-----------------------------------------------------------------------------#
#Max strains
max_strainFM = -(0.01*c)/l_FM
max_strainMT = -(0.02*c)/l_MT

list_strainsFM = np.linspace(0, max_strainFM, 10)
list_strainsMT = np.linspace(0, max_strainMT, 10)
plt.figure()

for strainFM in list_strainsFM:
    
    gamma = -(strainFM)*l_FM/R
    aux_gamma.append(np.degrees(gamma))

    
    for strainMT in list_strainsMT:
        #in all interations the contration of the spring is the same
        #calculating gamma, beta and delta based on SMA spring strain
        
        gamma = -(strainFM)*l_FM/R
        list_gamma.append(np.degrees(gamma))        
        
                  
        beta = -(strainMT)*l_MT/R                
        list_beta.append(np.degrees(beta))
        
        if strainFM == 0 :
            aux_beta.append(np.degrees(beta))
         
        
        M_gamma = {'x':M['x']*np.cos(gamma) - M['y']*np.sin(gamma),
           'y':M['x']*np.sin(gamma) + M['y']*np.cos(gamma)}
    
        
        T_gamma = {'x':T['x']*np.cos(gamma) - T['y']*np.sin(gamma),
           'y':T['x']*np.sin(gamma) + T['y']*np.cos(gamma)}
        
        T_beta = {'x':M_gamma['x'] + 
                      (T_gamma['x']-M_gamma['x'])*np.cos(beta) - 
                      (T_gamma['y']-M_gamma['y'])*np.sin(beta),
                  'y':M_gamma['y'] + 
                      (T_gamma['x']-M_gamma['x'])*np.sin(beta) + 
                      (T_gamma['y']-M_gamma['y'])*np.cos(beta)}
        
        if strainMT == 0: 
            list_wingtipheightgamma.append(T_beta['y'])
            
        if strainFM == 0:
            list_wingtipheightbeta.append(T_beta['y'])
        
        delta = np.arctan((T_beta['y']/T_beta['x']))
        list_delta.append(np.degrees(delta))
        
        
        
        plotter(M_gamma, T_beta, fraction_M = strainFM/max_strainFM,
                fraction_T = strainMT/max_strainMT)

plt.grid()


#DONE: plots of wingtip displacement versus deflection (gamma,beta)
def deflectionbygamma():
    #calculates for beta = 0
    plt.figure()
    plt.grid()
    plt.scatter(aux_gamma,list_wingtipheightgamma,color = 'c')
    plt.plot(aux_gamma, list_wingtipheightgamma,color = 'c' )
    plt.xlabel('gamma (considering beta = 0) ')
    plt.ylabel('tip height')
    
def deflectionbybeta():
    #calculates for gamma = 0
    plt.figure()
    plt.grid()
    plt.scatter(aux_beta,list_wingtipheightbeta ,color = 'g')
    plt.plot(aux_beta, list_wingtipheightbeta ,color = 'g' )
    plt.xlabel('beta (considering gamma = 0)')
    plt.ylabel('tip height')


#DONE: plots of wingtip displacement versus strain (strainFM, strainMT)


def deflectionbystrainFM():
    #calculates for beta = 0
    plt.figure()
    plt.grid()
    plt.scatter(list_strainsFM,list_wingtipheightgamma)
    plt.plot(list_strainsFM,list_wingtipheightgamma)
    plt.xlabel('FM strains')
    plt.ylabel('wingtipheight (for beta = 0)')
    
def deflectionbystrainMT():
    #calculates for gamma = 0
    plt.figure()
    plt.grid()
    plt.scatter(list_strainsMT,list_wingtipheightbeta, color = [0.8,0.2,0.7])
    plt.plot(list_strainsMT,list_wingtipheightbeta, color = [0.8,0.2,0.7])
    plt.xlabel('MT strains')
    plt.ylabel('wingtipheight (for gamma = 0)')    



#TODO: plots of deflections (relative, gamma, beta) versus strains



#plotting functions

deflectionbygamma()
deflectionbybeta()
deflectionbystrainFM()
deflectionbystrainMT()

'''
print  'gamma: ', list_gamma
print  'beta: ', list_beta
print 'delta: ' , list_delta
'''

