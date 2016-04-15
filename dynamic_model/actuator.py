# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 17:27:40 2016

@author: Pedro Leal
"""
import math
from scipy.optimize import newton
import numpy as np

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
    def __init__(self, geo_props, J, area = None, zero_stress_length = None,
                 eps_0 = None, k = None, material = 'linear', design = 'A'):
        """
        Initiate class and it's basic atributes:
        - geoprops: dictionary with the following keys:
            - x_n and y_n: are the x & y coordinates of the forwards end in the
              local coordinate system (origin at the joint). n stands for
              negative
            - x_p and y_p: are the x & y coordinates of the afterwards end in the
              local coordinate system (origin at the joint). p stand for 
              positive
            - area
        - design:
            - A: em X
            - B: pulley
            - C: co-linear
        - eps_0: inital strain (epsilon), if not defined, it is zero
        - k: linear spring elastic coefficient
        """
        #Storing inputs in local coordinate system
        self.x_n = geo_props['x-'] - J['x']  #Local Coordinates
        self.y_n = geo_props['y-'] - J['y'] #Local Coordinates
        self.x_p = geo_props['x+'] - J['x'] #Local Coordinates
        self.y_p = geo_props['y+'] - J['y'] #Local Coordinates
        self.x_J = J['x']                   #Global Coordinates
        self.y_J = J['y']                   #Global Coordinates
        self.material = material
        self.design = design
        
        # Projections of original wire length r
        self.r_1_0 = self.x_p - self.x_n
        self.r_2_0 = self.y_p - self.y_n
        self.length_r_0 = math.sqrt(self.r_1_0**2 + self.r_2_0**2)    
        
        # Initializing current wire length r
        self.r_1 = self.r_1_0
        self.r_2 = self.r_2_0
        self.length_r = self.length_r_0
        
        #calculate rigid body distances
        self.r_n = math.sqrt(self.x_n**2 + self.y_n**2)
        self.r_p = math.sqrt(self.x_p**2 + self.y_p**2)
        
        #define initial values for r, theta and force
        self.theta = 0
        self.F = 0.
        #Cross section area
        try:
            self.area = geo_props['area']
        except:
            self.area = area
        self.sigma = None
        # type of actuator
        try:
            self.actuator_type = geo_props['actuator_type']
        except:
            self.actuator_type = "wire"
        
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
        Calculate angle for given deformation epsilon via the newton method.
        """
        #TODO: include option for design B
        if self.design == "A":
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
    #        print self.eps, self.eps_0
            self.theta = newton(eq, theta_0, diff_eq, maxiter = 100)
            
            if abs(self.theta) > math.pi:
                self.theta = self.theta % (2.*math.pi)
        return self.theta
        
    def update(self, theta = None):
        """If vector length r or theta is changed, new coordinates are 
        calculated"""
        if theta != None:
            self.theta = theta
        else:
            #TODO: Add option for B
            if self.design == 'A':
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
        return self.F      
          
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
        import matplotlib.pyplot as plt
        if self.material == 'linear':
            colour = 'b'
        elif self.material == 'SMA':
            colour = 'r'
        plt.figure(1)
        plt.axes().set_aspect('equal')
        
        if self.actuator_type == "wire":
            plt.scatter([self.x_n + self.x_J, self.x_n + self.x_J + self.r_1], 
                        [self.y_n + self.y_J, self.y_n + self.y_J + self.r_2], 
                        c=colour)
            plt.scatter([self.x_J],[self.y_J], c = 'g')
            plt.plot([self.x_n + self.x_J, self.x_n + self.x_J + self.r_1], 
                     [self.y_n + self.y_J, self.y_n + self.y_J + self.r_2],
                     colour)
        #TODO: include spring plot

    def find_limits(self, y, theta_0 = 0):
        """The actuator has two major constraints:
            A - Because there is no physical sense of an actuator that has any
        part of it outside of the aircraft. We need to find the maximum
        theta and eps the actuator can have without this taking place.
            B - When r+ and r- are aligned, but + is between - and J, we 
            have the minimum length possible for the actuator. Below this,
            it is quite unrealistic
            The maximum and minimum theta is defined by the smallest of
            theta_A and theta_B
        """
 
        def diff_eq(theta):
            sin = math.sin(theta)
            cos = math.cos(theta)
            diff = -a*(-self.x_p*sin - self.y_p*cos) + self.x_p*cos - self.y_p*sin
            return diff
        
        def eq_theta_A(theta):
            sin = math.sin(theta)
            cos = math.cos(theta)
            return -a*(self.x_p*cos - self.y_p*sin - 0) + \
                    self.x_p*sin + self.y_p*cos - y['l']
        
        def eq_theta_B():
            A = 2. - math.sqrt(self.x_p**2 + self.y_p**2)/math.sqrt(self.x_n**2 + self.y_n**2)
            sin = A * (self.y_n - self.x_n*self.y_p/self.x_p)/(self.x_p + self.y_p**2/self.x_p)
            cos = (A*self.x_n + self.y_p*sin)/self.x_p
            return math.atan2(sin, cos)

        # Constraint B
        if self.r_n > self.r_p:
            self.max_theta_B = math.atan2(self.y_n*self.x_p - self.x_n*self.y_p,
                                          self.x_n*self.x_p + self.y_n*self.y_p)
            self.max_theta_B = np.sign(self.max_theta_B) * (abs(self.max_theta_B) % (2*math.pi))

            if self.max_theta_B > 0.:
                self.min_theta_B = self.max_theta_B
                self.max_theta_B = self.max_theta_B - 2*math.pi
            else:
                self.min_theta_B = self.max_theta_B + 2*math.pi

        else:
            self.max_theta_B = -math.pi/2.
            self.min_theta_B = math.pi/2.
            
        # Constraint A
        #Avoid division by zero for when x_n is the origin
        if abs(self.x_n) > 1e-4:
#            print 'comparison', self.r_p, abs(y['l'])
            if self.r_p >= abs(y['l']):
                a = (y['l'] - self.y_n)/(0. - self.x_n)
                self.max_theta_A = newton(eq_theta_A, theta_0, diff_eq, maxiter = 1000)
            else:
                self.max_theta_A = -math.pi/2.
            if self.r_p >= abs(y['u']):
                a = (y['u'] - self.y_n)/(0. - self.x_n)
                self.min_theta_A = newton(eq_theta_A, theta_0, diff_eq, maxiter = 1000)
            else:
                self.min_theta_A = math.pi/2.

        else:
            self.max_theta_A = -math.pi/2.
            self.min_theta_A = math.pi/2.

        self.max_theta_A = np.sign(self.max_theta_A) * (abs(self.max_theta_A) % (2*math.pi))
        self.min_theta_A = np.sign(self.min_theta_A) * (abs(self.min_theta_A) % (2*math.pi))

        self.max_theta = max(self.max_theta_A, self.max_theta_B)
        self.min_theta = min(self.min_theta_A, self.min_theta_B)         

        #In case of full transformation, we have the maximum eps        
        r_1 = self.x_p*math.cos(self.max_theta) - \
                   self.y_p*math.sin(self.max_theta) - self.x_n
        r_2 = self.y_p*math.cos(self.max_theta) + \
                   self.x_p*math.sin(self.max_theta) - self.y_n
        length_r = math.sqrt(r_1**2 + r_2**2)
        self.max_eps = length_r/self.zero_stress_length - 1. 
        

        r_1 = self.x_p*math.cos(self.min_theta) - \
                   self.y_p*math.sin(self.min_theta) - self.x_n
        r_2 = self.y_p*math.cos(self.min_theta) + \
                   self.x_p*math.sin(self.min_theta) - self.y_n
        length_r = math.sqrt(r_1**2 + r_2**2)
        self.min_eps = length_r/self.zero_stress_length - 1. 
    
    def check_crossing_joint(self, tol = 1e-3):
        """Does the actuator cross the joint? Should not happen"""
        
        #rotationed 
        x_p = self.x_p*math.cos(self.theta) - \
                   self.y_p*math.sin(self.theta)
        y_p = self.y_p*math.cos(self.theta) + \
                   self.x_p*math.sin(self.theta)
                   
        B = (y_p - self.y_n)/(x_p - self.x_n)
        y_at_J = self.y_n - B* self.x_n
        
        if abs(y_at_J) < tol:
            return True
        else:
            return False
#        print 'theta_A', self.max_theta_A, self.min_theta_A
#        print 'theta_B', self.max_theta_B, self.min_theta_B

        #To constraint the secant method, I need to know the global
        
#        def diff_eps(theta):
#            sin = math.sin(theta)
#            cos = math.cos(theta)
#            
#            diff = (2.*self.x_p*cos - 2.*self.y_p*sin)*(self.x_p*sin + \
#                    self.y_p*cos - self.y_n) - (2.*self.x_p*sin + \
#                    2.*self.y_p*cos)*(self.x_p*cos - self.x_n - self.y_p*sin) 
#            
#            return diff
#            
#        self.theta_max_eps = newton(diff_eps, self.eps_0)
#        r_1 = self.x_p*math.cos(self.theta) - \
#                   self.y_p*math.sin(self.theta) - self.x_n
#        r_2 = self.y_p*math.cos(self.theta) + \
#                   self.x_p*math.sin(self.theta) - self.y_n
#        length_r = math.sqrt(r_1**2 + r_2**2)
#        self.global_max_eps = length_r/self.zero_stress_length - 1.
        