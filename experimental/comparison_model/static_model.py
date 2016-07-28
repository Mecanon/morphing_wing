# -*- coding: utf-8 -*-
"""
- dynamics of a flap with two actuators in different positions
- can take in to account initial strain
- calculates necessary spring stiffness for a given linear actuator length
  and for defined positions for the actuators ends.

Will have:
- coupling with Edwin matlab code
- coupling with Aeropy
Created on Wed Feb 17 13:10:30 2016

@author: Pedro Leal
"""
import math
import numpy as np
import pickle

import airfoil_module as af
from flap import flap

def run(inputs, parameters = None):
    """Function to be callled by DOE and optimization. Design Variables are 
        the only inputs.
        
        :param inputs: {'sma', 'linear', 'sigma_o'}"""
    def thickness(x, t, chord):
        y = af.Naca00XX(chord, t, [x], return_dict = 'y')
        thickness_at_x = y['u'] - y['l']
        return thickness_at_x 

    if parameters != None:
        eng = parameters[0]
        import_matlab = False
    else:
        eng = None
        import_matlab = True
        
    sma = inputs['sma']
    linear = inputs['linear']
    
    spring = inputs['spring']
           
    airfoil = "naca0012"
    chord = 1.#0.6175
    t = 0.12*chord

    J = {'x':0.75, 'y':0.}
    
#    # need to transform normalized coordiantes in to global coordinates
#    sma['y+'] = sma['y+']*thickness(sma['x+'], t, chord)/2.
#    sma['y-'] = sma['y-']*thickness(sma['x-'], t, chord)/2.
#    
#    linear['y+'] = linear['y+']*thickness(linear['x+'], t, chord)/2.
#    linear['y-'] =  linear['y-']*thickness(linear['x-'], t, chord)/2.
    
    #Adding the area key to the dictionaries
    sma['area'] = math.pi*(0.000381/2)**2
    linear['area'] = 0.001
    
    # Design constants   
    #arm length to center of gravity
    r_w = 0.1
    
    #Aicraft weight (mass times gravity)
    W = 0.0523*9.8 #0.06*9.8
    alpha = 0.
    V = 10 #m/s
    altitude = 10000. #feet
    
    # Temperature
    T_0 = 273.15 + 30
    T_final = 273.15 + 140
     
    #Initial martensitic volume fraction
    MVF_init = 1.
    
    # Number of steps and cycles
    n = 200
    n_cycles = 0
    #~~~~~~~~~~~~~~~~~~~~~bb~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Parameters to select how to output stuff
    all_outputs = True
    save_data = True
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if all_outputs:
#        print "SMA: real y dimensions: ", sma['y-'], sma['y+'],sma['y+']*thickness(sma['x+'], t, chord)/2., sma['y-']*thickness(sma['x-'], t, chord)/2.
#        print "linear: real y dimensions: ", linear['y-'], linear['y+'], linear['y+']*thickness(linear['x+'], t, chord)/2., linear['y-']*thickness(linear['x-'], t, chord)/2.
        eps_s, eps_l, theta, sigma, MVF, T, eps_t, theta, F_l, k, L_s = flap(airfoil, 
                               chord, J, sma, linear, spring, 
                               W, r_w, V, altitude, alpha, T_0, 
                               T_final, MVF_init, n, all_outputs = True,
                               import_matlab = import_matlab, eng=eng,
                               n_cycles = n_cycles, aero_loads = False)

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(np.rad2deg(theta), eps_s, lw=2., label = "$\epsilon_s$")
        plt.plot(np.rad2deg(theta), eps_l, 'b--',lw=2, label = "$\epsilon_l$")
#        plt.scatter(theta, eps_s, c = 'b')
#        plt.scatter(theta, eps_l, c = 'b')
        plt.ylabel('$\epsilon$', fontsize=24)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        print len(T), len(eps_s), len(eps_l), len(theta), len(eps_t)
        plt.figure()
        plt.plot(np.rad2deg(theta), eps_t, lw=2.)
#        plt.scatter(theta, eps_t, c = 'b')
        plt.ylabel('$\epsilon_t$', fontsize=24)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        plt.figure()
        plt.plot(np.rad2deg(theta), MVF, lw=2.)
#        plt.scatter(theta, MVF, c = 'b')
        plt.ylabel('$MVF$', fontsize=24)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()

        plt.figure()
        plt.plot(T, MVF, lw=2.)
#        plt.scatter(T, MVF, c = 'b')
        plt.ylabel('$MVF$', fontsize=24)
        plt.xlabel('$T (K)$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()

        plt.figure()
        plt.plot(T, sigma, lw=2.)
#        plt.scatter(T, sigma, c = 'b')
        plt.ylabel('$\sigma$', fontsize=24)
        plt.xlabel('$T (K)$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        plt.figure()
        plt.plot(T, eps_s, 'b', lw=2., label = "$\epsilon_s$")
        plt.plot(T, eps_l, 'b--',lw=2, label = "$\epsilon_l$")
#        plt.scatter(T, eps_s, c = 'b')
#        plt.scatter(T, eps_l, c = 'b')
        plt.xlabel('$T (K)$', fontsize=20)
        plt.ylabel('$\epsilon$', fontsize=24)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        plt.figure()
        plt.plot(T, np.rad2deg(theta), lw=2.)
#        plt.scatter(T, theta, c = 'b')
        plt.xlabel('$T (K)$', fontsize=20)
        plt.ylabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.grid()
        
        F_s = []
        for i in range(len(sigma)):
            F_s.append(sigma[i]*sma['area'])
#        sigma_MPa = []
#        for sigma_i in sigma:
#            sigma_MPa.append(sigma_i/1e6)
        plt.figure()
        plt.plot(np.rad2deg(theta), F_s, 'b', lw=2., label = "$F_s$")
        plt.plot(np.rad2deg(theta), F_l, 'b--', lw=2., label = "$F_l$")
#        plt.scatter(theta, F_s, c = 'b')
#        plt.scatter(theta, F_l, c = 'b')
        plt.ylabel('$F (N)$', fontsize=20)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()        
    else:
        theta, k= flap(airfoil, chord, J, sma, linear, spring, 
                       W, r_w, V, altitude, alpha, T_0, 
                       T_final, MVF_init, n, all_outputs = False,
                       import_matlab = import_matlab, eng=eng,
                       n_cycles = n_cycles, aero_loads = False)
    
    if save_data == True:
        Data = {'theta': theta, 'eps_s': eps_s, 'eps_l': eps_l, 
                'sigma': sigma, 'xi': MVF, 'T': T, 'eps_t': eps_t,
                'F_l': F_l, 'k': k, 'L_s':L_s}
        pickle.dump(Data, open( "data.p", "wb" ) )
    
    return {'theta': theta, 'k': k}

if __name__ == '__main__':
    J = {'x':0.75, 'y':0.}

    # Real positioning on flap model (y coordinates are not normalized)          
    sma = {'x-': J['x'] - 0.018, 'y-': -0.02, 
           'x+': J['x'] + 0.152, 'y+': 0.006}
    linear = {'x-': J['x'] - 0.018, 'y-': 0.03, 
           'x+': J['x'] + 0.152, 'y+': -0.006}

     
    # Spring stiffness (N/m)
    k = 127.33 
    
    spring = {'k': k, 'L_0': 0.0588, 'L_solid': 0.05030}
    data = run({'sma':sma, 'linear':linear, 'spring':spring})
    print  'k: ', data['k'], 'theta:', data['theta'][-1]
    DataFile = open('data.txt','a')
							
##==============================================================================
## Run withou run function
##==============================================================================
#    #Hole positioning
#    J = {'x':0.25, 'y':0.}
#    #y coordinates are percentual
#    sma = {'x-': J['x'], 'y-': -0.02*2, 'x+': 0.1225 + J['x'],
#           'y+': 0.0135*2, 'area':math.pi*0.00025**2}
#    linear = {'x-': J['x'], 'y-': 0.032*2, 'x+': 0.146 + J['x'], 
#              'y+': -0.0135*2, 'area':0.001}
#    
#    #original bias spring length
#    length_l = 0.06 #
#    
#    #arm length to center of gravity
#    r_w = 0.15
#    
#    #Aicraft weight (mass times gravity)
#    W = 0.06*9.8
#    alpha = 0.
#    V = 10 #m/s
#    altitude = 10000. #feet
#    
#    airfoil = "naca0012"
#    chord = 0.6175
#    
#    ## Temperature
#    T_0 = 220.15
#    T_final = 400.15
#     
#    #Initial martensitic volume fraction
#    MVF_init = 1.
#    
#    # Number of steps
#    n = 200
#    
#    data = flap(airfoil, chord, J, sma, linear, sigma_o, length_l, W, r_w,
#                V, altitude, alpha, T_0, T_final, MVF_init, n,
#                all_outputs = True)
    