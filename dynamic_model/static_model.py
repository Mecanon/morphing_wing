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
from scipy.interpolate import interp1d

import airfoil_module as af
from flap import flap
from flap_multiobjective import flap_multiobjective

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
    sigma_o = 100e6
           
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
    sma['area'] = math.pi*(0.000381/2.)**2
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
    n = 1000
    n_cycles = 1
    #~~~~~~~~~~~~~~~~~~~~~bb~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Parameters to select how to output stuff
    all_outputs = True
    save_data = False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if all_outputs:
#        print "SMA: real y dimensions: ", sma['y-'], sma['y+'],sma['y+']*thickness(sma['x+'], t, chord)/2., sma['y-']*thickness(sma['x-'], t, chord)/2.
#        print "linear: real y dimensions: ", linear['y-'], linear['y+'], linear['y+']*thickness(linear['x+'], t, chord)/2., linear['y-']*thickness(linear['x-'], t, chord)/2.
        eps_s, eps_l, theta, sigma, MVF, T, eps_t, theta, F_l, k, L_s, H_cur = flap(airfoil, 
                               chord, J, sma, linear, sigma_o, 
                               W, r_w, V, altitude, alpha, T_0, 
                               T_final, MVF_init, n, all_outputs = True,
                               import_matlab = import_matlab, eng=eng,
                               n_cycles = n_cycles)

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
        plt.plot(T, H_cur, lw=2.)
        plt.ylabel('$H_{cur}$', fontsize=24)
        plt.xlabel('T', fontsize=20)
        plt.legend(loc = 'best', fontsize = 'x-large')
        plt.grid()

        sigma_crit = 0
        H_max = 0.0550
        H_min = 0.0387
        k = 4.6849e-09
    
        plt.figure()
        plt.plot(sigma, H_cur, lw=2.)
        plt.plot(sigma, H_min + (H_max - H_min)*(1. - np.exp(-k*(abs(np.array(sigma)) - \
                                                     sigma_crit))), 'k', lw=2.)
        plt.ylabel('$H_{cur}$', fontsize=24)
        plt.xlabel('$\sigma$', fontsize=20)
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

def run_multiobjective(inputs, parameters = None):
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
    sigma_o = inputs['sigma_o']
           
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
    sma['area'] = math.pi*(0.000381/2.)**2
    linear['area'] = 0.001
    
    # Design constants   
    #arm length to center of gravity
    r_w = 0.10
    
    #Aicraft weight (mass times gravity)
    W = 0.0523*9.8 #0.06*9.8
    alpha = 0.
    V = 10 #m/s
    altitude = 10000. #feet
    
    # Temperature
    T_0 = 273.15 + 30
    T_final = inputs['T_f']
     
    #Initial martensitic volume fraction
    MVF_init = 1.
    
    # Number of steps and cycles
    n = 200
    n_cycles = 0
    #~~~~~~~~~~~~~~~~~~~~~bb~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Parameters to select how to output stuff
    all_outputs = True
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if all_outputs:
        eps_s, eps_l, theta, sigma, MVF, T, eps_t, theta, F_l, k, L_s = flap_multiobjective(airfoil, 
                               chord, J, sma, linear, sigma_o, 
                               W, r_w, V, altitude, alpha, T_0, 
                               T_final, MVF_init, n, all_outputs = True,
                               import_matlab = import_matlab, eng=eng,
                               n_cycles = n_cycles)

        return theta, sigma, T, MVF, eps_s, L_s

def power(delta_t, sigma, T, xi, eps_s, L_s, output = "all"):
    """
    Calculate work, power and current.
    
    - output: defines what is the function output (Power or all)
    """
    sigma_o = 100e6
    r = 0.000381/2.
    d = 2*r
    T_o = 273.15 + 30
    
    alpha = 0.           #set to zero on purpose
    c = 320.              #invented
    rho = 6450.
    
    #Transformation strain properties
    H_max = 0.0550
    H_min = 0.0387
    sigma_crit = 0
    k = 4.6849e-09
    
    rho_E_M = 0.8e-6         #Dynalloy
    rho_E_A = 1.0e-6         #Dynalloy
    E_A = 3.7427e+10
    E_M = 8.8888e+10
    C_A = 7.9498e+06
    C_M = 7.1986e+06
    M_s = 363.5013
    M_f = 297.9735
    A_s = 324.6427
    A_f = 385.0014
    n1 = 0.1752
    n2 = 0.1789
    n3 = 0.1497
    n4 = 0.2935
    sigma_cal = 200E6
    
    #==============================================================================
    # # Heat Transfer parameters
    #==============================================================================
    # Gravity:
    g = 9.8 #ms-2
    # Atmospheric pressure
    P_air = 101325. # Pa
    # Molar
    M = 0.0289644  #kg/mol
    # Ideal gas constant
    R = 8.31447  #J/(mol K)
    # Air density:
    rho_air = P_air*M / (R*T_o)
    # Sutherland's law coefficients
    C1 = 1.458e-6 #kg/m.s.sqrt(K)
    C2 = 110.4 #K
    # Air dynamic viscosity:
    mu_air = (C1 * T_o**(3./2)) / (T_o+C2)
    # Air kinematic viscosity:
    nu_air = mu_air/rho_air
    # Air specific heat at constant pressure
    CP_list = [1.0038, 1.0049, 1.0063, 1.0082, 1.0106, 1.0135, 1.0206]
    T_list = [275., 300., 325., 350., 375., 400., 450.]
    Cp_f = interp1d(T_list, CP_list)
    # Air conductivity
    k_list = [2.428e-5, 2.624e-5, 2.816e-5, 3.003e-5, 3.186e-5, 3.365e-5, 3.710e-5]
    k_f = interp1d(T_list, k_list)
    
    # Nusselt number coefficients
    alpha_1 = 1.
    alpha_2 = 0.287
    
    #==============================================================================
    # Calculate Power and current
    #==============================================================================
    I_list = []
    P_list = []
    W_list = []
    n = len(eps_s)
    for i in range(1, n):
        delta_sigma = sigma[i] - sigma[i-1]
        delta_T = T[i] - T[i-1]
        delta_eps = eps_s[i] - eps_s[i-1]
        delta_xi = xi[i] - xi[i-1]
        
        T_avg = (T[i] + T[i-1])/2.
        Cp_air = Cp_f(T_avg)
        k_air = k_f(T_avg)
        # Grashof number for external flow around a cylinder
        Gr = 2*abs(T[i] - T_o)/(T[i] + T_o)*(g*d**3)/(nu_air**2)
        # Prandtl number definition
        Pr = mu_air*Cp_air/k_air
        # Nusselt number and parameter
        Nu = (alpha_1 + alpha_2*(Gr*Pr/(1 + (0.56/Pr)**(9./16))**(16./9))**(1./6))**2
        # Calculate convection coefficient h from definition of Nusselt number
        h = k_air*Nu/d
        
        rho_E = rho_E_M*xi[i] + (1-xi[i])*rho_E_A
        
        if abs(sigma[i]) <= sigma_crit:
            dH_cur = 0
        else:
            dH_cur = k*(H_max-H_min)*math.exp(-k*(abs(sigma[i])-sigma_crit))*np.sign(sigma[i])
        H_cur = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - sigma_crit)))
        H_cur_cal = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_cal) - sigma_crit)))
        
        rho_delta_s0 = (-2*(C_M*C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/(C_M + C_A)
        a1 = rho_delta_s0*(M_f - M_s)
        a2 = rho_delta_s0*(A_s - A_f)
        a3 = -a1/4 * (1 + 1/(n1+1) - 1/(n2+1)) + a2/4 * (1+1/(n3+1) - 1/(n4+1))
        Y_0_t = rho_delta_s0/2*(M_s - A_f) - a3
        D = ((C_M - C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/((C_M + C_A)*(H_cur_cal+ sigma_cal*dH_cur))
    
        pi_t = Y_0_t + D*abs(sigma[i])*H_cur
    
        #constant h
        P = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma + \
            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        
        P_list.append(P)
        
        if output == 'all':
            I = r*math.pi*math.sqrt((r/rho_E)*((r/delta_t)*((T[i]*alpha*delta_sigma + \
                rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) ) + \
                2.*h*(T[i] - T_o))))
        
            
            dW = math.pi*r**2*L_s[0]*0.5*(sigma[i]+sigma[i-1])*delta_eps
            
            I_list.append(I)
            W_list.append(dW)
        
    Total_power = 0
    for i in range(len(P_list)-1):
        Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
    if output == 'all':
        return I_list, P_list, W_list, Total_power
    elif output == "power":
        return Total_power
        
if __name__ == '__main__':
    J = {'x':0.75, 'y':0.}
    # Position coordinates from holes. y coordinates are a fraction of thickness/2.

    #Optimal A from max deflection             
#    sma = {'x-': 4.389066e-001, 'y-': -8.311361e-001, 
#           'x+': 7.990382e-001, 'y+': 6.039162e-002}
#    linear = {'x-': 7.323110e-001, 'y-': 7.573718e-001, 
#           'x+': 8.543053e-001, 'y+': -2.499118e-001}
																 
#    #Optimal C from max deflection             
    sma = {'x-': 3.941320e-001, 'y-': -8.647118e-001, 
           'x+': 8.116175e-001, 'y+': 3.137898e-002 }
    linear = {'x-': 3.941320e-001, 'y-': -8.647118e-001, 
           'x+': 8.116175e-001, 'y+': 3.137898e-002 }

#    sma = {'x-': 0.72316, 'y-': -0.75730, 
#           'x+': 0.75844, 'y+': 0.06584}
#    linear = {'x-': 0.43045, 'y-': 0.32455, 
#           'x+': 0.81779, 'y+': -0.09255 }
     
    T_f = 358.66849
    data = run({'sma':sma, 'linear':linear})
#    theta, sigma, T, MVF, eps_s, L_s= run_multiobjective({'sma':sma, 'linear':linear, 'T_f':T_f})
#    print 'theta: ', theta[-1], 'T:', T[-1]
#    delta_t = 0.05   
#    
#    P = power(delta_t, sigma, T, MVF, eps_s, L_s, output = "power")
#    print  'P: ', P
    print 'theta:', data['theta'][-1]							
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
    