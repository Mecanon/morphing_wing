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

#Load data
Data = pickle.load(open( "data.p", "rb" ))

sigma = Data['sigma']
T = Data['T']
xi = Data['xi']
eps_s = Data['eps_s']
L_s = Data['L_s']

n = len(eps_s)

h = 10.               #invented (adiabatic)
h_list = np.linspace(0,100., 5)
P_h_list = []
total_power_list = []
for j in range(len(h_list)):
    h = h_list[j]
    P_list = []
    I_list = []
    for i in range(1, n):
        delta_sigma = sigma[i] - sigma[i-1]
        delta_T = T[i] - T[i-1]
        
        delta_xi = xi[i] - xi[i-1]
        
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
        I_list.append(I)
        P_list.append(P)
        
    P_h_list.append(P_list)
    Total_power = 0
    for i in range(len(P_list)-1):
        Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
    total_power_list.append(Total_power)

t = np.linspace(0,(n-2)*delta_t, n-1)

#plt.figure()
#plt.plot(t, I_list, 'b')
#plt.scatter(t, I_list, c = 'b')
#plt.xlabel('Time (s)')
#plt.ylabel('Current (A)')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(I_list) - 0.02*(max(I_list)-min(I_list)),
#          max(I_list) + 0.02*(max(I_list)-min(I_list))])
#plt.grid()

plt.figure()
for i in range(len(h_list)):
    color=((1.-float(i)/(len(h_list)-1), float(i)/(len(h_list)-1),0, 1.))
    plt.plot(t, P_h_list[i], 'b', label = 'h = ' + str(h_list[i]), color = color)
#plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.legend(loc= 'upper left')

plt.figure()
plt.plot(h_list, total_power_list)
plt.xlabel('Convection coefficient')
plt.ylabel('Total power consumption (J)')
plt.grid()

plt.figure()
plt.plot(h_list, 100.*total_power_list[0]/np.array(total_power_list))
plt.xlabel('Convection coefficient $h$ ')
plt.ylabel('Percentage of power spent on transformation (%)')
plt.grid()

print 'Total power is %f Joules' % Total_power


h = 10.               #invented (adiabatic)
T_list = np.linspace(200,300., 5)
P_T_list = []
total_power_list = []
for j in range(len(T_list)):
    T_o = T_list[j]
    P_list = []
    I_list = []
    for i in range(1, n):
        delta_sigma = sigma[i] - sigma[i-1]
        delta_T = T[i] - T[i-1]
        
        delta_xi = xi[i] - xi[i-1]
        
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
#        I = r*math.pi*math.sqrt((r/rho_E)*((r/delta_t)*((T[i]*alpha*delta_sigma + \
#            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) ) + \
#            2.*h*(T[i] - T_o))))
#    
        P = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma + \
            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        I_list.append(I)
        P_list.append(P)
        
    P_T_list.append(P_list)
    Total_power = 0
    for i in range(len(P_list)-1):
        Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
    total_power_list.append(Total_power)
    
plt.figure()
for i in range(len(T_list)):
    color = ((1.-float(i)/(len(T_list)-1), float(i)/(len(T_list)-1),0, 1.))
    plt.plot(t, P_T_list[i], 'b', label = '$T_o$ = ' + str(T_list[i]), color = color)
#plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.legend(loc= 'upper left')

plt.figure()
plt.plot(T_list, total_power_list)
plt.xlabel('Temperature (K)')
plt.ylabel('Total power consumption (J)')
plt.grid()