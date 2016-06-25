# -*- coding: utf-8 -*-
"""
Analyze the heating, current and power usage of teh actuation
Created on Thu Apr 28 09:56:23 2016

@author: Pedro Leal
"""
import math
import numpy as np
import pickle
import matplotlib.pyplot as plt

#Time step
delta_t = 0.05

sigma_o = 400e6
r = 0.00025
d = 2*r
T_o = 200.

alpha = 0.           #set to zero on purpose
c = 0.               #invented
rho = 1.             #invented

#Transformation strain properties
H_max = 0.04
H_min = 0.
sigma_crit = 140e6
k = 0.021e-6  

rho_E_M = 1.         #invented
rho_E_A = 1.         #invented
E_A = 60E9
E_M = 60E9
C_A = 7.8E6
C_M = 7.3e6
M_s = 333.
M_f = 220
A_s = 274.
A_f = 370.
n1 = 0.06
n2 = 0.06
n3 = 0.06
n4 =0.06
sigma_cal = 300E6

#Load data
Data = pickle.load(open( "data.p", "rb" ))

sigma = Data['sigma']
T = Data['T']
xi = Data['xi']
eps_s = Data['eps_s']
L_s = Data['L_s']

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
    
t = np.linspace(0,(n-2)*delta_t, n-1)

plt.figure()
plt.plot(t, I_list, 'b')
plt.scatter(t, I_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
          min(I_list) - 0.02*(max(I_list)-min(I_list)),
          max(I_list) + 0.02*(max(I_list)-min(I_list))])
plt.grid()

plt.figure()
plt.plot(t, P_list, 'b')
plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
          min(P_list) - 0.02*(max(P_list)-min(P_list)),
          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()

Total_power = 0
for i in range(len(P_list)-1):
    Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
print 'Total power is %f Joules' % Total_power