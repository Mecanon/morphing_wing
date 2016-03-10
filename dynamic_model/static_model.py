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
from scipy.optimize import newton

class actuator():
    """
    Actuator object where inputs:
    - -(n): is a dictionary with keys 'x' and 'y', the coordinates of 
      the forwards end in the global coordinate system (origin at the
      leading edge)
    - +(p): is a dictionary with keys 'x' and 'y', the coordinates of 
      the afterwards end in the global coordinate system (origin at the
      leading edge)
    - J: is a dictionary with keys 'x' and 'y', the coordinates of 
      the joint in the global coordinate system (origin at the
      leading edge)        
    - material: linear or SMA
    """
    #
    def __init__(self, positions, J, area = None, zero_stress_length = None,
                 eps_0 = None, k = None, material = 'linear'):
        """
        Initiate class and it's basic atributes:
        - x_n and y_n: are the x & y coordinates of the forwards end in the
          local coordinate system (origin at the joint). n stands for
          negative
        - x_p and y_p: are the x & y coordinates of the afterwards end in the
          local coordinate system (origin at the joint). p stand for 
          positive
        - eps_0: inital strain (epsilon), if not defined, it is zero
        - k: linear spring elastic coefficient
        """
        #Storing inputs in local coordinate system
        self.x_n= positions['x-'] - J['x']
        self.y_n = positions['y-'] - J['y']
        self.x_p = positions['x+'] - J['x']
        self.y_p = positions['y+'] - J['y']
        self.x_J = J['x']
        self.y_J = J['y']
        self.material = material
        
        # Projections of original wire length r
        self.r_1_0 = self.x_p - self.x_n
        self.r_2_0 = self.y_p - self.y_n
        self.length_r_0 = math.sqrt(self.r_1_0**2 + self.r_2_0**2)    
        
        # Initializing current wire length r
        self.r_1 = self.r_1_0
        self.r_2 = self.r_2_0
        self.length_r = self.length_r_0
        
        #define initial values for r, theta and force
        self.theta = 0
        self.F = 0.
        #Cross section area
        self.area = area
        self.sigma = None
        
        #In case zero_stress_length is not defined and initial strain is
        #known, calculate the zero stress length. If defined, calculate
        #the initial strain
        
        if zero_stress_length == None and eps_0 != None:
            self.eps_0 = eps_0
            self.eps = self.eps_0
            self.zero_stress_length = self.length_r_0/(1 + self.eps)
        elif zero_stress_length != None and eps_0 == None:
            self.zero_stress_length = zero_stress_length
            self.eps_0 = self.length_r/self.zero_stress_length - 1.
            self.eps = self.eps_0            
        
        if material == 'linear':
            self.k = k
    def calculate_theta(self, theta_0 = 0.):
        """
        Calculate angle for given deformation epsilon.
        """
        def diff_eq(theta):
            sin = math.sin(theta)
            cos = math.cos(theta)
            
            diff = (2.*self.x_p*cos - 2.*self.y_p*sin)*(self.x_p*sin + \
                    self.y_p*cos - self.y_n) - (2.*self.x_p*sin + \
                    2.*self.y_p*cos)*(self.x_p*cos - self.x_n - self.y_p*sin) 
            
            return diff
        
        def eq(theta):
            eta = (self.eps - self.eps_0) + 1.
            sin = math.sin(theta)
            cos = math.cos(theta)
            
            r_1 = self.x_p*cos - self.y_p*sin - self.x_n
            r_2 = self.y_p*cos + self.x_p*sin - self.y_n

            return r_1**2 + r_2**2 - (eta*self.length_r_0)**2

        self.theta = newton(eq, theta_0, diff_eq, maxiter = 200)
        return self.theta
        
    def update(self):
        """If vector length r or theta is changed, new coordinates are 
        calculated"""
        self.r_1 = self.x_p*math.cos(self.theta) - \
                   self.y_p*math.sin(self.theta) - self.x_n
        self.r_2 = self.y_p*math.cos(self.theta) + \
                   self.x_p*math.sin(self.theta) - self.y_n
        self.length_r = math.sqrt(self.r_1**2 + self.r_2**2)
        self.eps = self.length_r/self.zero_stress_length - 1. 
     
    def calculate_force(self, source = 'strain'):
        
        if source == 'strain':
            if self.material == 'linear':
                self.F = self.k*self.eps*self.length_r
            elif self.material == 'SMA':
                print "Put Edwin code here"
        #Calculate force from stress and cross section
        elif source == 'sigma':
            self.F = self.area * self.sigma
                
    def calculate_torque(self):
        """Calculate torque given the actuator force: r \times F (where a is 
        global coordinates)"""

        #calculate components of force vector        
        F_1 = - self.F*self.r_1/self.length_r
        F_2 = - self.F*self.r_2/self.length_r
#        print x_n, y_n
        
        #calculate torque
        self.torque = (self.x_p*math.cos(self.theta) - \
                       self.y_p*math.sin(self.theta))*F_2 - \
                      (self.y_p*math.cos(self.theta) + \
                       self.x_p*math.sin(self.theta))*F_1    
        return self.torque
    def plot_actuator(self):
        if self.material == 'linear':
            colour = 'b'
        elif self.material == 'SMA':
            colour = 'r'
        plt.figure(1)
        plt.axes().set_aspect('equal')
        plt.scatter([self.x_n + self.x_J, self.x_p + self.x_J], 
                    [self.y_n + self.y_J, self.y_p + self.y_J], 
                    c=colour)
        plt.scatter([self.x_J],[self.y_J], c = 'g')
        plt.plot([self.x_n + self.x_J, self.x_p + self.x_J], 
                 [self.y_n + self.y_J, self.y_p + self.y_J],
                 colour)
#    def calculate_initial(self, tension = "sigmaf"):
#        """Calculate initial displacement for SMA and linear"""
#        if self.material == "linear":
#            if tension == "sigmaf":
#                self.eps = self.area * self.sigmaf *(self.)
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    chord = 0.6175

#==============================================================================
# Design variables
#==============================================================================
    #Hole positioning
    J = {'x':0.25, 'y':0.}
    sma = {'x-': J['x'], 'y-': -0.02, 'x+': 0.1225 + J['x'], 'y+': 0.0135 }
    linear = {'x-': J['x'], 'y-': 0.032, 'x+': 0.146 + J['x'], 'y+': -0.0135}

    #SMA Pre-stress
    sigma_o = 300e6
#==============================================================================
# Design constants
#==============================================================================
    #Areas
    area_l = 0.001
    area_s = math.pi*0.00025**2
    
    #original bias spring length
    length_l = 0.06 #
    
    #SMA properties
    E_M = 60E9
    sigma_crit = 140e6
    H_max = 0.04
    H_min = 0.
    k = 0.021e-6
    
    #If I change the value of the material young modulus, need to update k
    k = k*E_M/60E9
    #arm length to center of gravity
    r_w = 0.15
    
    #Aicraft weight (mass times gravity)
    W = 0.06*9.8
    
#==============================================================================
# Initial conditions   
#==============================================================================
    #Initial transformation strain
    eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - sigma_crit)))
    #Define initial strain
    eps_0 = eps_t_0 + sigma_o/E_M
    #Linear actuator (s)
    l = actuator(linear, J, area_l, zero_stress_length = length_l, material = 'linear')
    #Sma actuator (l)
    s = actuator(sma, J, area_s, eps_0 = eps_0, material = 'SMA')

#    plt.figure(1)
#    s.plot_actuator()
#    l.plot_actuator()

#    print 'linear', l.length_r, l.eps 
#    print 'sma', s.length_r, s.eps
#    
#    #Input initial stress   
#    s.sigma = sigma_o
#    s.calculate_force(source = 'sigma')
#    
#    # TODO: For now K is calculated in function of the rest. Afterwards it will
#    # be an input
#
#    l.k = ((s.F/s.length_r)*(s.x_p*s.r_2 - s.y_p*s.r_1) + r_w*W)/(l.eps*(l.y_p*l.r_1 - l.x_p*l.r_2))
#    l.calculate_force(source = 'strain')
#    
#    s.theta = l.calculate_theta()
#
#    #Calculate initial torques
#    s.calculate_torque()    
#    l.calculate_torque()
#    
##==============================================================================
## Test of deformation per theta
##==============================================================================
#    import numpy as np
#    import matplotlib.pyplot as plt
#    
#    theta_list = np.linspace(0, -math.pi/4.)
#    eps_s_list = []
#    eps_l_list = []
#    
#    for theta in theta_list:
#        s.theta = theta
#        l.theta = theta
#        
#        s.update()
#        l.update()
#        
##        s.calculate_force()
#        l.calculate_force()
#        
#        eps_s_list.append(s.eps)
#        eps_l_list.append(l.eps)
##    plt.figure()    
##    plt.plot(np.degrees(theta_list), eps_s_list, 'r', np.degrees(theta_list), eps_l_list, 'b')  
##    plt.xlabel('$\\theta (degrees)$')
##    plt.ylabel('$\epsilon$')
##    BREAK
#    eps_s_list_test = eps_s_list
#    theta_list_test = theta_list
###==============================================================================
## Matlab simulation
###==============================================================================
##    import matlab.engine
##    from scipy.optimize import newton
##    #If the derivate for the newton function is not defined, it uses the
##    #secant method
##    
##    #Start Matlab engine
##    eng = matlab.engine.start_matlab()
##    #Go to directory where matlab file is
##    eng.cd('SMA_temperature_strain_driven')
#
#    def run(T_0, T_final, MVF_init, i, n, eps, eps_t_0, sigma_0 = 0,
#            eps_0 = 0, plot = 'True'):
#        """Run SMA model
#        
#        - all inputs are scalars"""
#        k = i+1
#        if k == n:
#            data = eng.OneD_SMA_Model(k, eps, T_0, T_final, MVF_init, 
#                                      eps_t_0, sigma_0, eps_0, n, plot, nargout=5) 
#        else:
#            data = eng.OneD_SMA_Model(k, eps, T_0, T_final, MVF_init,
#                                      eps_t_0, sigma_0, eps_0, n, 'False', nargout=5)
#        return data
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
#    def equlibrium(eps_s, s, l, T_0, T_final, MVF_init, sigma_0,
#                   i, n, r_w, W):
#        """Calculates the moment equilibrium. Function used for the 
#        secant method.
#        """
#        
#        #calculate new theta for eps_s and update all the parameter
#        #of the actuator class
#        s.eps = eps_s
#        s.calculate_theta()
#        s.update()
#        
#        l.theta = s.theta
#        l.update()
#        
#        #SMA (Constitutive equation: coupling via sigma)
#        data = run(T_0, T_final, MVF_init, i, n, eps_s, eps_t_0, sigma_0,
#                   s.eps_0, plot = 'False')
#        
#        s.sigma = data[0][i][0]
#        s.calculate_force(source = 'sigma')
#        tau_s = s.calculate_torque()
#        
#        #Linear (Geometric equation: coupling via theta)
#        l.calculate_force()
#        tau_l = l.calculate_torque()
#        
#        #weight (Geometric equation: coupling via theta)
#        tau_w = - r_w*W*math.cos(l.theta)
#        
#        return tau_s + tau_l + tau_w
#    eps_s = eps_0
#    eps_s_list = [eps_s]
#    theta_list = [s.theta]
#    for i in range(1, n):
#        eps_s = newton(equlibrium, x0 = eps_s, args = ((s, l, T_0, T_final, MVF_init, sigma_o,
#                                   i, n, r_w, W,)), maxiter = 500, tol = 1.0e-6)
#        eps_s_list.append(eps_s)
##        s.eps = eps_s
##        theta_list.append(math.degrees(s.calculate_theta()))
##        print i, eps_s, s.theta
#    #Extra run with prescribed deformation (which has already been calculated)
#    # to get all the properties
#    for i in range(1, n):
#        data = run(T_0, T_final, MVF_init, i, n, eps_s_list[i], eps_t_0, sigma_0 = sigma_o,
#                eps_0 = eps_0, plot = 'False')
#        
#    for eps_s in eps_s_list:
#        s.eps = eps_s
#        s.calculate_theta()
#        s.update()
#        theta_list.append(s.theta)
#        print eps_s, math.degrees(s.theta)
###        result = equlibrium(eps_s, s, l, T_0, T_final, MVF_init, sigma_o,
###                   i, n, r_w, W)
##
###        #Calculate new linear strain and theta for new force
###        conv_tol = 1e-1
###        conv_error = 1.
###        it_counter = 0
###        dampner = 0.
###        while conv_error > conv_tol:
###            prev_sigma = s.sigma  
###            
###            #calculate new SMA stress
###            data = run(T_0, T_final, MVF_init, i, n, s.eps, eps_t_0, sigma_o, eps_0, plot = 'False')
###            s.sigma = round(data[0][i][0], 1)
###            print 'SMA stress: ', s.sigma
###            s.calculate_force(source = 'sigma')
###            s.calculate_torque()  
###            
###            prev_theta = l.theta
###            prev_eps = l.eps
###            prev_eps_s = s.eps
###           
###            cos = math.cos(l.theta)
###            sin = math.sin(l.theta)
###            l.eps = dampner*prev_eps + (1-dampner)*((s.F/s.length_r)*((s.x_p*cos - s.y_p*sin)*s.r_2 - \
###                    (s.x_p*sin + s.y_p*cos)*s.r_1) + r_w*W)/(l.k* \
###                    ((l.x_p*sin + l.y_p*cos)*l.r_1 - (l.x_p*cos - \
###                    l.y_p*sin)*l.r_2)) 
###
###            #Calculate new theta
###            l.calculate_theta()
###            s.theta = l.theta
###            #update all values related to theta
###
###            l.update()
###            s.update()
###            s.calculate_torque()
###            
###            s.eps = s.eps
###            conv_error = abs(prev_sigma- s.sigma)
###            print it_counter, l.eps, s.eps
###            it_counter += 1
###            if it_counter == 10:
###                break
###
###        print i, 's_eps: ', s.eps, 'l_eps: ', l.eps, 'theta: ', l.theta, 'sigma: ', s.sigma, 'error: ', conv_error, it_counter
##        
##        
##        #update angle, r and strain at SMA actuator
###        theta = l.calculate_theta()
###        s.theta = theta
###        s.update()
###        l.update()
###        #strain of SMA actuator
###        s.calculate_force()
###        s.calculate_torque()
###        counter += 1
###        
###        #calculating error with final value 
###        error_eps = abs(previous_l_eps - l.eps)
###        #updating previous strain
###        previous_l_eps = l.eps
####        damping_counter +=1
###        print 'iteration: ', counter
###        print 'strain: ', s.eps, l.eps
###        print 'angle: ', s.theta, l.theta
###        print 'torque: ', s.torque
##        print "error: ", error_eps