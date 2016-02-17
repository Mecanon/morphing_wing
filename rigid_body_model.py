# -*- coding: utf-8 -*-
"""
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
    - material: 
    """
    #
    def __init__(self, F, A, J, material = 'elastic'):
        """
        Initiate class and it's basic atributes:
        - f_1 and f_2: are the x & y coordinates of the forwards end in the
          local coordinate system (origin at the joint)
        - a_1 and a_2: are the x & y coordinates of the afterwards end in the
          local coordinate system (origin at the joint)
        """
        #Storing inputs in local coordinate system
        self.f_1 = F['x'] - J['x']
        self.f_2 = F['y'] - J['y']
        self.a_1 = A['x'] - J['x']
        self.a_2 = A['y'] - J['y']
        self.material = material
        
        # Projections of original wire length r
        self.rx_0 = self.a_1 - self.f_1
        self.ry_0 = self.a_2 - self.f_2
                
    def angle(self, eta = 1.):
        """
        Calculate angle for given deformation eta.
        """
        A1 = eta*self.rx_0 + self.f_1
        A2 = eta*self.ry_0 + self.f_2
        cos = (A2 + A1*(self.a_1/self.a_2))/(self.a_2 +(self.a_1/self.a_2)*self.a_1)
        sin = (self.a_1*cos - A1) / self.a_2
        print sin, cos
        theta = math.atan2(sin,cos)
        return math.degrees(theta)
    
if __name__ == '__main__':
    chord = 1.
    F1 = {'x': -1*chord, 'y': 1*chord}
    A1 = {'x': 1*chord, 'y': 1*chord}
    F2 = {'x': 0.4*chord, 'y': -0.2*chord}
    A2 = {'x': 0.6*chord, 'y': 0.2*chord}
    J = {'x': 0.*chord, 'y': 0.*chord}
    
    spring = actuator(F1, A1, J, material = 'elasticity')
    
    print spring.angle()