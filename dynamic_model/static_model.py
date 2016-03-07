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

class actuator():
    """
    Actuator object where inputs:
    - F: is a dictionary with keys 'x' and 'y', the coordinates of 
      the forwards end in the global coordinate system (origin at the
      leading edge)
    - A: is a dictionary with keys 'x' and 'y', the coordinates of 
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
    def calculate_theta(self):
        """
        Calculate angle for given deformation epsilon.
        """
        #deformation = strain + 1
        eta = (self.eps - self.eps_0) + 1.

        A1 = eta*self.r_1_0 + self.x_n
        
        A2 = eta*self.r_2_0 + self.y_n
        cos = (A2 + A1*(self.x_p/self.y_p))/(self.y_p +(self.x_p/self.y_p)*self.x_p)
        sin = (self.x_p*cos - A1) / self.y_p
#        print sin, cos
        self.theta = math.degrees(math.atan2(sin,cos))
#        self.update()
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
                  
    def calculate_torque(self):
        """Calculate torque given the actuator force: a \times F (where a is 
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
     
    def calculate_force(self, source = 'strain'):
        
        if source == 'strain':
            if self.material == 'linear':
                self.F = self.k*self.eps*self.length_r
            elif self.material == 'SMA':
                print "Put Edwin code here"
        #Calculate force from stress and cross section
        elif source == 'sigma':
            self.F = self.area * self.sigma
#    def calculate_initial(self, tension = "sigmaf"):
#        """Calculate initial displacement for SMA and linear"""
#        if self.material == "linear":
#            if tension == "sigmaf":
#                self.eps = self.area * self.sigmaf *(self.)
if __name__ == '__main__':
    chord = 1.

#==============================================================================
# Design variables
#==============================================================================
    #Hole positioning
    sma = {'x-': -1*chord, 'y-': -2*chord, 'x+': 1*chord, 'y+': 1*chord }
    linear = {'x-': -1.*chord, 'y-':-2*chord, 'x+': 1*chord, 'y+': -1*chord}
    J = {'x':0., 'y':0.}
    #SMA Pre-stress
    sigma_o = 300e6
#==============================================================================
# Design constants
#==============================================================================
    #Areas
    area_l = 0.001
    area_s = 0.001
    
    #original bias spring length
    length_l = 2.
    
    #Spring coefficient
    k = 100.
    
    #SMA properties
    E_M = 60E9
    sigma_crit = 140e6
    H_max = 0.04
    H_min = 0.
    k = 0.021e-6
    
    #arm length to center of gravity
    a_w = 1.
    
    #Aicraft weight
    W = 1.
#==============================================================================
# Initial conditions   
#==============================================================================
    #Define initial strain
    eps_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - sigma_crit)))
    #Linear actuator (s)
    l = actuator(linear, J, area_l, zero_stress_length = length_l, material = 'linear')
    #Sma actuator (l)
    s = actuator(sma, J, area_s, eps_0 = eps_0, material = 'SMA')

    print 'linear', l.length_r, l.eps 
    print 'sma', s.length_r, s.eps
    
    #Input initial stress   
    s.sigma = sigma_o
    s.calculate_force(source = 'sigma')
    
    # TODO: For now K is calculated in function of the rest. Afterwards it will
    # be an input
    s.calculate_torque()
    print s.torque
    l.k = ((s.F/s.length_r)*(s.x_p*s.r_2 - s.y_p*s.r_1) + a_w*W)/(l.eps*(l.y_p*l.r_1 - l.x_p*l.r_2))
    l.calculate_force(source = 'strain')
    
    s.theta = l.calculate_theta()

    print s.calculate_torque()    
    print l.calculate_torque()
    
    
#        #strain (epsilon) of linear actuator
#        l.eps = dampner*previous_l_eps + (1.-dampner)*(1./k)/(l.x_p*math.sin(theta) + l.y_p*math.cos(theta)*l.r_1 - \
#              (l.x_p*math.cos(theta) - l.y_p*math.sin(theta))*l.r_2)*(s.torque + \
#              a_w*W*math.cos(theta))
        #update angle, r and strain at SMA actuator
#        theta = l.calculate_theta()
#        s.theta = theta
#        s.update()
#        l.update()
#        #strain of SMA actuator
#        s.calculate_force()
#        s.calculate_torque()
#        counter += 1
#        
#        #calculating error with final value 
#        error_eps = abs(previous_l_eps - l.eps)
#        #updating previous strain
#        previous_l_eps = l.eps
##        damping_counter +=1
#        print 'iteration: ', counter
#        print 'strain: ', s.eps, l.eps
#        print 'angle: ', s.theta, l.theta
#        print 'torque: ', s.torque
#        print "error: ", error_eps