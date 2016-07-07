# -*- coding: utf-8 -*-
"""
Analyze the heating, current and power usage of teh actuation
Created on Thu Apr 28 09:56:23 2016

@author: Pedro Leal
"""
import math
import numpy as np
import matplotlib.pyplot as plt

def power(delta_t, sigma, T, xi, eps_s, L_s, output = "all"):
    """
    Calculate work, power and current.
    
    - output: defines what is the function output (Power or all)
    """
    sigma_o = 100e6
    r = 0.00025
    d = 2*r
    T_o = 200.
    
    alpha = 0.           #set to zero on purpose
    c = 837.36              #invented
    rho = 6450.
    
    #Transformation strain properties
    H_max = 0.1209
    H_min = 0.0924
    sigma_crit = 0
    k = 5.9713e-09 
    
    rho_E_M = 0.8e-6         #Dynalloy
    rho_E_A = 1.0e-6         #Dynalloy
    E_A = 2.1496e+10
    E_M = 3.3453e+10
    C_A = 8.0370e+06
    C_M = 7.1233e+06
    M_s = 362.5851
    M_f = 297.4771
    A_s = 318.3625
    A_f = 386.8458
    n1 = 0.1919
    n2 = 0.1823
    n3 = 0.1623
    n4 = 0.2188
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
    Cp_air = 1.005
    # Air conductivity
    k_air = 0.0264
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
        I = r*math.pi*math.sqrt((r/rho_E)*((r/delta_t)*((T[i]*alpha*delta_sigma + \
            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) ) + \
            2.*h*(T[i] - T_o))))
    
        P = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma + \
            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        
        dW = math.pi*r**2*L_s[0]*0.5*(sigma[i]+sigma[i-1])*delta_eps
        
        I_list.append(I)
        P_list.append(P)
        W_list.append(dW)
        
    Total_power = 0
    for i in range(len(P_list)-1):
        Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
    if output == 'all':
        return I_list, P_list, W_list, Total_power
    elif output == "power":
        return Total_power
    
if __name__ == '__main__':
    import pickle
    #Load data
    Data = pickle.load(open( "data.p", "rb" ))
    
    sigma = Data['sigma']
    T = Data['T']
    xi = Data['xi']
    eps_s = Data['eps_s']
    L_s = Data['L_s']
    
    #Time step
    delta_t = 0.05    
    
    I, P, W, Total_power = power(delta_t, sigma, T, xi, eps_s, L_s, output = "all")
    n = len(eps_s)
    t = np.linspace(0,(n-2)*delta_t, n-1)
    
    plt.figure()
    plt.plot(t, I, 'b')
    plt.scatter(t, I, c = 'b')
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
              min(I) - 0.02*(max(I)-min(I)),
              max(I) + 0.02*(max(I)-min(I))])
    plt.grid()
    
    plt.figure()
    plt.plot(t, P, 'b')
    plt.scatter(t, P, c = 'b')
    plt.xlabel('Time (s)')
    plt.ylabel('Power (W)')
    plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
              min(P) - 0.02*(max(P)-min(P)),
              max(P) + 0.02*(max(P)-min(P))])
    plt.grid()
    
    print 'Total power is %f Joules' % Total_power