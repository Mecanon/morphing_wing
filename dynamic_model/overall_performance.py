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

def run(inputs, parameters = None, W = 0.2*9.8, alpha = 0., V = 10):
    """Function to be callled by DOE and optimization. Design Variables are 
        the only inputs.
        
        :param inputs: {'sma', 'linear', 'sigma_o'}
        :param W: weight (N)
        :param alpha: angle of attack(radians)
        :param V: velocity (m/s)
        """
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
    sigma_o = 400e6
           
    airfoil = "naca0012"
    chord = 1.#0.6175
    t = 0.12*chord

    J = {'x':0.75, 'y':0.}
    
    # need to transform normalized coordiantes in to global coordinates
    sma['y+'] = sma['y+']*thickness(sma['x+'], t, chord)/2.
    sma['y-'] = sma['y-']*thickness(sma['x-'], t, chord)/2.
    
    linear['y+'] = linear['y+']*thickness(linear['x+'], t, chord)/2.
    linear['y-'] =  linear['y-']*thickness(linear['x-'], t, chord)/2.
    
    #Adding the area key to the dictionaries
    sma['area'] = math.pi*0.00025**2
    linear['area'] = 0.001
    
    # Design constants   
    #arm length to center of gravity
    r_w = 0.15
    
    #Aicraft weight (mass times gravity)
    altitude = 10000. #feet
    
    # Temperature
    T_0 = 220.
    T_final = 420.
     
    #Initial martensitic volume fraction
    MVF_init = 1.
    
    # Number of steps and cycles
    n = 200
    n_cycles = 0
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Parameters to select how to output stuff
    all_outputs = False
    save_data = False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if all_outputs:
        eps_s, eps_l, theta, sigma, MVF, T, eps_t, theta, F_l, k, L_s = flap(airfoil, 
                               chord, J, sma, linear, sigma_o, 
                               W, r_w, V, altitude, alpha, T_0, 
                               T_final, MVF_init, n, all_outputs = True,
                               import_matlab = import_matlab, eng=eng,
                               n_cycles = n_cycles)

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(theta, eps_s, lw=2., label = "$\epsilon_s$")
        plt.plot(theta, eps_l, 'b--',lw=2, label = "$\epsilon_l$")
#        plt.scatter(theta, eps_s, c = 'b')
#        plt.scatter(theta, eps_l, c = 'b')
        plt.ylabel('$\epsilon$', fontsize=24)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        print len(T), len(eps_s), len(eps_l), len(theta), len(eps_t)
        plt.figure()
        plt.plot(theta, eps_t, lw=2.)
#        plt.scatter(theta, eps_t, c = 'b')
        plt.ylabel('$\epsilon_t$', fontsize=24)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()
        
        plt.figure()
        plt.plot(theta, MVF, lw=2.)
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
        plt.plot(T, theta, lw=2.)
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
        plt.plot(theta, F_s, 'b', lw=2., label = "$F_s$")
        plt.plot(theta, F_l, 'b--', lw=2., label = "$F_l$")
#        plt.scatter(theta, F_s, c = 'b')
#        plt.scatter(theta, F_l, c = 'b')
        plt.ylabel('$F (N)$', fontsize=20)
        plt.xlabel(r'$\theta ({}^{\circ})$', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()        
    else:
        theta, k= flap(airfoil, chord, J, sma, linear, sigma_o, 
                       W, r_w, V, altitude, alpha, T_0, 
                       T_final, MVF_init, n, all_outputs = False,
                       import_matlab = import_matlab, eng=eng,
                       n_cycles = n_cycles)
    
    if save_data == True:
        Data = {'theta': theta, 'eps_s': eps_s, 'eps_l': eps_l, 
                'sigma': sigma, 'xi': MVF, 'T': T, 'eps_t': eps_t,
                'F_l': F_l, 'k': k, 'L_s':L_s}
        pickle.dump(Data, open( "data.p", "wb" ) )
    
    return {'theta': theta, 'k': k}
    
if __name__ == '__main__':
    J = {'x':0.75, 'y':0.}
    # Position coordinates from holes. y coordinates are a fraction of thickness/2.

    #Optimal for max deflection
#    sma = {'x-': 6.817445e-001, 'y-': -5.216475e-001, 
#           'x+': 9.029895e-001, 'y+': 8.726738e-001}
#    linear =  {'x-': 6.958111e-001, 'y-': -4.593744e-001, 
#               'x+': 8.187166e-001, 'y+': -5.719241e-001}
               
#    linear = {'x+': 0.8746717751201235, 'y+': -0.37669650699461776, 'y-': 0.7445609404288144, 'x-': 0.7394757699162359}
#    sma = {'x+': 0.8793685619240611, 'y+': 0.713701519873574, 'y-': -0.7509778773146336, 'x-': 0.6785498036887645}
#    #Optimal multiobjective               
#    sma = {'x-': 6.161543e-001, 'y-': -6.631015e-001, 
#           'x+': 8.697452e-001, 'y+': 3.962915e-001}
#    linear =  {'x-': 4.593649e-001, 'y-': -7.127816e-001, 
#               'x+': 8.269874e-001, 'y+': -1.587640e-001}

#    #Extension test           
#    sma = {'x-': 0.6, 'y-': -0.6, 
#           'x+': 0.8, 'y+': -0.6}
#    linear =  {'x-': 0.6, 'y-': 0.6, 
#               'x+': 0.8, 'y+': 0.6}

    #Compression test with n_real = n              
    sma = {'x-': 7.407724e-001, 'y-': -3.680615e-001, 
           'x+': 9.933211e-001, 'y+': 6.004423e-001}
    linear = {'x-': 7.290939e-001, 'y-': -7.584186e-001, 
           'x+': 7.550874e-001, 'y+': -4.011175e-001}

    #SMA Pre-stress
    sigma_o = 400e6

    V_list = np.linspace(0,35,10)
    alpha_list = np.linspace(0,math.radians(10.),10)
    W_list = np.linspace(0,5,10)
    
    data_V = {'theta':[],'k':[]}
    for V_i in V_list:
        data = run({'sma':sma, 'linear':linear, 'sigma_o':sigma_o},
                   V = V_i)
        data_V['theta'].append(data['theta'])
        data_V['k'].append(data['k'])  
    
    data_alpha = {'theta':[],'k':[]}
    for alpha_i in alpha_list:
        data = run({'sma':sma, 'linear':linear, 'sigma_o':sigma_o},
                   alpha = alpha_i)
        data_alpha['theta'].append(data['theta'])
        data_alpha['k'].append(data['k'])   
    
    data_W = {'theta':[],'k':[]}
    for W_i in W_list:
        data = run({'sma':sma, 'linear':linear, 'sigma_o':sigma_o},
                   W = W_i)
        data_W['theta'].append(data['theta'])
        data_W['k'].append(data['k'])
        
    pickle.dump([data_V, data_alpha, data_W], open( "data_overall.p", "wb" ) )    
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
    