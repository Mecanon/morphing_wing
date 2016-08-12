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

r = 0.000381/2.
d = 2*r

alpha = 0.           #set to zero on purpose
c = 320.           
rho = 6450.

#Transformation strain properties
H_max = 0.0550
H_min = 0.0387
sigma_crit = 0
k = 4.6849e-09

rho_E_M = 0.8352e-6
rho_E_A = 0.9343e-6
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
#Load data
Data = pickle.load(open( "data.p", "rb" ))

sigma = Data['sigma']
T = Data['T']
xi = Data['xi']
eps_s = Data['eps_s']
L_s = Data['L_s']

T_o = T[0]

n = len(eps_s)

#==============================================================================
# Calculate output work
#==============================================================================
W_list = []
deltaW_list = []
Total_work = 0
total_work_list = []
for i in range(1, n-1):
    delta_eps_f = abs(eps_s[i+1] - eps_s[i])
    delta_eps_b = abs(eps_s[i] - eps_s[i-1])
    delta_sigma_f = abs(sigma[i+1] - sigma[i-1])
    delta_sigma_b = abs(sigma[i] - sigma[i-1])
#    avg_eps = abs(eps_s[i] + eps_s[i-1])/2.
#    avg_sigma = abs(sigma[i] + sigma[i-1])/2.
    
    dW = math.pi*r**2*L_s[0]*0.5*(eps_s[i+1]*delta_sigma_f  + eps_s[i]*delta_sigma_b + \
                                  sigma[i+1]*delta_eps_f + sigma[i]*delta_eps_b)/delta_t

#    deltaW = math.pi*r**2*L_s[0]*(eps_s[i]*delta_sigma/delta_t + sigma[i]*delta_eps/delta_t)
    W_list.append(dW*delta_t)
#    deltaW_list.append(deltaW)
#    Total_work += deltaW
    total_work_list.append(Total_work)
Total_work = sum(W_list)
Total_delta = Total_work
deltaW_list = W_list
#print Total_delta

#==============================================================================
# Calculate input heat for different h
#==============================================================================
h_list = np.linspace(0,100., 6)
P_h_list = []
total_power_list = []

thermoelasctic_comp = []
specific_heat_comp = []
transformation_force_comp = []
transformation_entropy_comp = []
xi_rate = []
for j in range(len(h_list)):
    h = h_list[j]
    P_list = []
    if j == 0:
        I_list = []
    a = 0
    b = 0
    for i in range(1, n-1):
        delta_sigma_f = sigma[i+1] - sigma[i]
        delta_sigma_b = sigma[i] - sigma[i-1]
        
        delta_T_f = T[i+1] - T[i]
        delta_T_b = T[i] - T[i-1]
        
        delta_xi_f = xi[i+1] - xi[i]
        delta_xi_b = xi[i] - xi[i-1]
        rho_E = rho_E_M*xi[i] + (1-xi[i])*rho_E_A
        
        if abs(sigma[i]) <= sigma_crit:
            dH_cur = 0
        else:
            dH_cur = k*(H_max-H_min)*math.exp(-k*(abs(sigma[i])-sigma_crit))*np.sign(sigma[i])
        H_cur = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma[i]) - sigma_crit)))
        H_cur_cal = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_cal) - sigma_crit)))
        
        rho_delta_s0 = (-2*(C_M*C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/(C_M + C_A)
        a1 = rho_delta_s0*(M_f - M_s)
        a2 = rho_delta_s0*(A_s - A_f)
        a3 = -a1/4 * (1 + 1/(n1+1) - 1/(n2+1)) + a2/4 * (1+1/(n3+1) - 1/(n4+1))
        Y_0_t = rho_delta_s0/2*(M_s - A_f) - a3
        D = ((C_M - C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/((C_M + C_A)*(H_cur_cal+ sigma_cal*dH_cur))
    
        pi_t = Y_0_t + D*abs(sigma[i])*H_cur
    
        #constant h

    
        H_f = math.pi*r**2*L_s[i+1]*((T[i+1]*alpha*delta_sigma_f + \
            rho*c*delta_T_f + delta_xi_f*(-pi_t + rho_delta_s0*T[i+1]) )/delta_t + \
            2.*(h/r)*(T[i+1] - T_o))

        H_b = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma_b + \
            rho*c*delta_T_b + delta_xi_b*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        
        P = 0.5*(H_f + H_b)*delta_t
        
        I = P*rho_E*L_s[i]/((r**2)*math.pi)
        print P

#        print a,b

        if j == 0:
            I_list.append(I)
            front = math.pi*r**2*L_s[i+1]*(T[i+1]*alpha*delta_sigma_f)/delta_t
            back = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma_b))/delta_t
            thermoelasctic_comp.append(0.5*(front + back)*delta_t)
            
            front = math.pi*r**2*L_s[i+1]*(rho*c*delta_T_f)/delta_t
            back = math.pi*r**2*L_s[i]*(rho*c*delta_T_b)/delta_t
            specific_heat_comp.append(0.5*(front + back)*delta_t)
            
            front = math.pi*r**2*L_s[i+1]*(delta_xi_f*(-pi_t)/delta_t)
            back = math.pi*r**2*L_s[i]*(delta_xi_b*(-pi_t) )/delta_t
            transformation_force_comp.append(0.5*(front + back)*delta_t)
            
            front = math.pi*r**2*L_s[i+1]*(delta_xi_f*(rho_delta_s0*T[i+1]))/delta_t
            back = math.pi*r**2*L_s[i]*(delta_xi_b*(rho_delta_s0*T[i]) )/delta_t
            transformation_entropy_comp.append(0.5*(front + back)*delta_t)
            
            xi_rate.append((xi[i+1]-xi[i-1])/delta_t)
        P_list.append(P)
        
    P_h_list.append(P_list)
    Total_power = sum(P_list)
    total_power_list.append(Total_power)

t = np.linspace(0,(n-2)*delta_t, n-1)

plt.figure()
plt.plot(t[1:], I_list, 'b')
plt.scatter(t[1:], I_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
          min(I_list) - 0.02*(max(I_list)-min(I_list)),
          max(I_list) + 0.02*(max(I_list)-min(I_list))])
plt.grid()

plt.figure()
for i in range(len(h_list)):
    color=((1.-float(i)/(len(h_list)-1), float(i)/(len(h_list)-1),0, 1.))
    plt.plot(t[1:], P_h_list[i],  label = 'h = ' + str(h_list[i]), color = color)
plt.plot(t[1:], deltaW_list, 'b', label = '$\dot{W}$')
#plt.plot(t, W_list, 'b', label = '$\dot{W}$')
#plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.legend(loc= 'best')

plt.figure()
plt.subplot(211)
for i in range(len(h_list)):
    color=((1.-float(i)/(len(h_list)-1), float(i)/(len(h_list)-1),0, 1.))
    plt.plot(t[1:], P_h_list[i],  label = 'h = ' + str(h_list[i]), color = color)
plt.plot(t[1:], deltaW_list, 'b', label = '$\dot{W}$')
#plt.plot(t, W_list, 'b', label = '$\dot{W}$')
#plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
plt.grid()
plt.subplot(212)
plt.plot(t,xi[:-1])
plt.xlabel('Time (s)')
plt.ylabel(r'$\xi$')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.legend(loc= 'best')

plt.figure()
plt.subplot(711)
plt.plot(t[1:], P_h_list[0])
plt.plot(t[1:], deltaW_list, 'b', label = '$\dot{W}$')
#plt.plot(t, W_list, 'b', label = '$\dot{W}$')
#plt.scatter(t, P_list, c = 'b')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')
plt.grid()
plt.subplot(712)
plt.plot(t,xi[:-1])
plt.xlabel('Time (s)')
plt.ylabel(r'$\xi$')
plt.grid()
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.subplot(713)
plt.plot(t[1:],xi_rate)
plt.xlabel('Time (s)')
plt.ylabel(r'$\dot{\xi}$')
plt.grid()

plt.subplot(714)
plt.plot(t[1:],thermoelasctic_comp)
plt.xlabel('Time (s)')
plt.ylabel(r'$A L_s T \alpha_T \dot{\sigma}$')
plt.grid()
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.subplot(715)
plt.plot(t[1:],specific_heat_comp)
plt.xlabel('Time (s)')
plt.ylabel(r'$A L_s \rho c \dot{T}$')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.subplot(716)
plt.plot(t[1:],transformation_force_comp)
plt.xlabel('Time (s)')
plt.ylabel(r'$ - A L_s \pi^t \dot{\xi}$')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.subplot(717)
plt.plot(t[1:],transformation_entropy_comp)
plt.xlabel('Time (s)')
plt.ylabel(r'$A L_s \rho \Delta s_o T \dot{\xi}$')
#plt.axis([min(t) - 0.02*(max(t)-min(t)), max(t)+ 0.02*(max(t)-min(t)),
#          min(P_list) - 0.02*(max(P_list)-min(P_list)),
#          max(P_list) + 0.02*(max(P_list)-min(P_list))])
plt.grid()
plt.legend(loc= 'best')

plt.figure()
plt.plot(h_list, total_power_list)
plt.xlabel('Convection coefficient')
plt.ylabel('Total power consumption (J)')
plt.grid()

plt.figure()
plt.plot(h_list, 100.*Total_delta/np.array(total_power_list))
plt.xlabel('Convection coefficient $h$ ')
plt.ylabel('Efficiency (%)')
plt.grid()

print 'Total adiabatic power is %f Joules' % total_power_list[0]
print 'Total work is %f Joules' % Total_delta
print 'Adiabatic efficiency is %f ' % (Total_delta/total_power_list[0])

#==============================================================================
# Calculate input heat for different delta_t
#==============================================================================
delta_t_list = np.linspace(0.001,0.05, 50)
#h = 10.
h_dt_power_list = []

for i in range(len(h_list)):
    h = h_list[i]
    total_power_list = []
    for j in range(len(delta_t_list)):
        delta_t = delta_t_list[j]
        P_list = []
        I_list = []
        a = 0
        b = 0
        for i in range(1, n-1):
            delta_sigma_f = sigma[i+1] - sigma[i]
            delta_sigma_b = sigma[i] - sigma[i-1]
            
            delta_T_f = T[i+1] - T[i]
            delta_T_b = T[i] - T[i-1]
            
            delta_xi_f = xi[i+1] - xi[i]
            delta_xi_b = xi[i] - xi[i-1]
            rho_E = rho_E_M*xi[i] + (1-xi[i])*rho_E_A
            
            if abs(sigma[i]) <= sigma_crit:
                dH_cur = 0
            else:
                dH_cur = k*(H_max-H_min)*math.exp(-k*(abs(sigma[i])-sigma_crit))*np.sign(sigma[i])
            H_cur = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma[i]) - sigma_crit)))
            H_cur_cal = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_cal) - sigma_crit)))
            
            rho_delta_s0 = (-2*(C_M*C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/(C_M + C_A)
            a1 = rho_delta_s0*(M_f - M_s)
            a2 = rho_delta_s0*(A_s - A_f)
            a3 = -a1/4 * (1 + 1/(n1+1) - 1/(n2+1)) + a2/4 * (1+1/(n3+1) - 1/(n4+1))
            Y_0_t = rho_delta_s0/2*(M_s - A_f) - a3
            D = ((C_M - C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/((C_M + C_A)*(H_cur_cal+ sigma_cal*dH_cur))
        
            pi_t = Y_0_t + D*abs(sigma[i])*H_cur
        
            #constant h
    
        
            H_f = math.pi*r**2*L_s[i+1]*((T[i+1]*alpha*delta_sigma_f + \
                rho*c*delta_T_f + delta_xi_f*(-pi_t + rho_delta_s0*T[i+1]) )/delta_t + \
                2.*(h/r)*(T[i+1] - T_o))
    
            H_b = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma_b + \
                rho*c*delta_T_b + delta_xi_b*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
                2.*(h/r)*(T[i] - T_o))
            
            P = 0.5*(H_f + H_b)*delta_t
            
            I = (r**2)*math.pi*P/rho_E/L_s[i]
#            print a,b
            I_list.append(I)
            P_list.append(P)
            
        Total_power = sum(P_list)
        total_power_list.append(Total_power)
    h_dt_power_list.append(total_power_list)
t = np.linspace(0,(n-2)*delta_t, n-1)

plt.figure()
for i in range(len(h_list)):
#    print len(h_dt_power_list)
    color=((1.-float(i)/(len(h_list)-1), float(i)/(len(h_list)-1),0, 1.))
    plt.plot((n-1)*np.array(delta_t_list), 100.*Total_delta/np.array(h_dt_power_list[i]),
             color = color, label= 'h = %.f' % h_list[i])
plt.xlabel('Time (s) ')
plt.ylabel('Efficiency (%)')
plt.grid()
plt.legend(loc='best')
#==============================================================================
# Calculate heat input for different T_o
#==============================================================================
delta_t = 0.05
h = 10.               #invented (adiabatic)
T_list = np.linspace(200,300., 5)
P_T_list = []
total_power_list = []
for j in range(len(T_list)):
    T_o = T_list[j]
    P_list = []
    I_list = []
    for i in range(1, n-1):
        delta_sigma_f = sigma[i+1] - sigma[i]
        delta_sigma_b = sigma[i] - sigma[i-1]
        
        delta_T_f = T[i+1] - T[i]
        delta_T_b = T[i] - T[i-1]
        
        delta_xi_f = xi[i+1] - xi[i]
        delta_xi_b = xi[i] - xi[i-1]
        rho_E = rho_E_M*xi[i] + (1-xi[i])*rho_E_A
        
        if abs(sigma[i]) <= sigma_crit:
            dH_cur = 0
        else:
            dH_cur = k*(H_max-H_min)*math.exp(-k*(abs(sigma[i])-sigma_crit))*np.sign(sigma[i])
        H_cur = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma[i]) - sigma_crit)))
        H_cur_cal = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_cal) - sigma_crit)))
        
        rho_delta_s0 = (-2*(C_M*C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/(C_M + C_A)
        a1 = rho_delta_s0*(M_f - M_s)
        a2 = rho_delta_s0*(A_s - A_f)
        a3 = -a1/4 * (1 + 1/(n1+1) - 1/(n2+1)) + a2/4 * (1+1/(n3+1) - 1/(n4+1))
        Y_0_t = rho_delta_s0/2*(M_s - A_f) - a3
        D = ((C_M - C_A)*(H_cur_cal + sigma_cal*dH_cur + sigma_cal*(1/E_M - 1/E_A)))/((C_M + C_A)*(H_cur_cal+ sigma_cal*dH_cur))
    
        pi_t = Y_0_t + D*abs(sigma[i])*H_cur
    
        #constant h

    
        H_f = math.pi*r**2*L_s[i+1]*((T[i+1]*alpha*delta_sigma_f + \
            rho*c*delta_T_f + delta_xi_f*(-pi_t + rho_delta_s0*T[i+1]) )/delta_t + \
            2.*(h/r)*(T[i+1] - T_o))

        H_b = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma_b + \
            rho*c*delta_T_b + delta_xi_b*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        
        P = 0.5*(H_f + H_b)*delta_t
        
        I = (r**2)*math.pi*P/rho_E/L_s[i]
        
        I_list.append(I)
        P_list.append(P)
        
    P_T_list.append(P_list)
    Total_power = sum(P_list)
    total_power_list.append(Total_power)
    
plt.figure()
for i in range(len(T_list)):
    color = ((1.-float(i)/(len(T_list)-1), float(i)/(len(T_list)-1),0, 1.))
    plt.plot(t[1:], P_T_list[i], 'b', label = '$T_o$ = ' + str(T_list[i]), color = color)
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