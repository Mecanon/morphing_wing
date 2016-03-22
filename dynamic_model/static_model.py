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
import numpy as np
    
from AeroPy import calculate_flap_moment
from aero_module import air_properties
import airfoil_module as af
import xfoil_module as xf

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
        self.x_n= geo_props['x-'] - J['x']  #Local Coordinates
        self.y_n = geo_props['y-'] - J['y'] #Local Coordinates
        self.x_p = geo_props['x+'] - J['x'] #Local Coordinates
        self.y_p = geo_props['y+'] - J['y'] #Local Coordinates
        self.x_J = J['x']                   #Global Coordinates
        self.y_J = J['y']                   #Global Coordinates
        self.material = material
        
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
        self.area = geo_props['area']
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
        Calculate angle for given deformation epsilon via the newton method.
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
#        print self.eps, self.eps_0
        self.theta = newton(eq, theta_0, diff_eq, maxiter = 1000)
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
        import matplotlib.pyplot as plt
        if self.material == 'linear':
            colour = 'b'
        elif self.material == 'SMA':
            colour = 'r'
        plt.figure(1)
        plt.axes().set_aspect('equal')
        plt.scatter([self.x_n + self.x_J, self.x_n + self.x_J + self.r_1], 
                    [self.y_n + self.y_J, self.y_n + self.y_J + self.r_2], 
                    c=colour)
        plt.scatter([self.x_J],[self.y_J], c = 'g')
        plt.plot([self.x_n + self.x_J, self.x_n + self.x_J + self.r_1], 
                 [self.y_n + self.y_J, self.y_n + self.y_J + self.r_2],
                 colour)
   
    def find_limits(self, t = 0.12, c = 1., theta_0 = 0):
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
        y = af.Naca00XX(c, t, [self.x_J], return_dict = 'y')
        
        def diff_eq(theta):
            sin = math.sin(theta)
            cos = math.cos(theta)
            diff = -a*(-self.x_p*sin - self.y_p*cos) + self.x_p*cos - self.y_p*sin
            return diff
        
        def eq_theta_A(theta):
            sin = math.sin(theta)
            cos = math.cos(theta)
            return -a*(self.x_p*cos - self.y_p*sin - self.x_n) + \
                    self.x_p*sin + self.y_p*cos - self.y_n 
        
        def eq_theta_B():
            A = 2. - math.sqrt(self.x_p**2 + self.y_p**2)/math.sqrt(self.x_n**2 + self.y_n**2)
            sin = A * (self.y_n - self.x_n*self.y_p/self.x_p)/(self.x_p + self.y_p**2/self.x_p)
            cos = (A*self.x_n + self.y_p*sin)/self.x_p
            return math.atan2(sin, cos)

        # Constraint B
        self.max_theta_B = math.atan2(self.y_n*self.x_p - self.x_n*self.y_p,
                                      self.x_n*self.x_p - self.y_n*self.y_p)
        if abs(self.max_theta_B) > 2*math.pi:
            self.max_theta_B = self.max_theta_B % 2*math.pi
        if self.max_theta_B > 0.:
            self.min_theta_B = self.max_theta_B
            self.max_theta_B = self.max_theta_B - 2*math.pi
        else:
            self.min_theta_B = self.max_theta_B + 2*math.pi
        
        # Constraint A
        #Avoid division by zero for when x_n is the origin
        print 'comparison', self.r_p, y['l']
        if abs(self.x_n) > 1e-4:
            a = (y['l'] - self.y_n)/(0. - self.x_n)
            self.max_theta_A = newton(eq_theta_A, theta_0,  maxiter = 1000)
            a = (y['u'] - self.y_n)/(0. - self.x_n)
            self.min_theta_A = newton(eq_theta_A, theta_0,  maxiter = 1000)
        else:
            self.max_theta_A = -math.pi/4.
            self.min_theta_A = math.pi/4.

        if self.r_n >= self.r_p:
            self.max_theta = max(self.max_theta_A, self.max_theta_B)
            self.min_theta = min(self.min_theta_A, self.min_theta_B)
        else:
            self.max_theta = self.max_theta_A
            self.min_theta = self.min_theta_A           
#        print 'min', self.min_theta_A, self.min_theta_B

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
        
        print 'theta_A', self.max_theta_A, self.min_theta_A
        print 'theta_B', self.max_theta_B, self.min_theta_B
        
def flap(airfoil, chord, J, sma, linear, sigma_o, length_l, W, r_w, V,
         altitude, alpha, T_0, T_final, MVF_init, n, all_outputs = False,
         import_matlab = True, eng = None, aero_loads = True):
    """
    solve actuation problem for flap driven by an antagonistic mechanism
    using SMA and linear actuators
    :param J: dictionary with coordinates of Joint
    :param sigma_o: pre-stress, if not defined, it calculates the biggest
                    one respecting the constraints (max is H_max)
    :param T_0: initial temperature
    :param T_final: final temperature
    :param MVF_init: initial martensitic volume fraction
    :param n: number of steps in simulation"""

    
    from scipy.optimize import newton    

    if import_matlab:
        import matlab.engine
        #Start Matlab engine
        eng = matlab.engine.start_matlab()
        #Go to directory where matlab file is
        if import_matlab:
            eng.cd(' ..')
            eng.cd('SMA_temperature_strain_driven')
        else:
            eng.cd('SMA_temperature_strain_driven')
            
    def constitutive_model(T_0, T_final, MVF_init, i, n, eps, eps_t_0, sigma_0 = 0,
            eps_0 = 0, plot = 'True'):
        """Run SMA model
        
        - all inputs are scalars"""
        k = i+1
        if k == n:
            data = eng.OneD_SMA_Model(k, eps, T_0, T_final, MVF_init, 
                                      eps_t_0, sigma_0, eps_0, n, plot, nargout=5) 
        else:
            data = eng.OneD_SMA_Model(k, eps, T_0, T_final, MVF_init,
                                      eps_t_0, sigma_0, eps_0, n, 'False', nargout=5)
        return data

    def equilibrium(eps_s, s, l, T_0, T_final, MVF_init, sigma_0,
                   i, n, r_w, W, x = None, y = None, alpha = 0.,
                   q = 1., chord = 1., x_hinge = 0.25, aero_loads = aero_loads):
        """Calculates the moment equilibrium. Function used for the 
        secant method.
        """
        
        #calculate new theta for eps_s and update all the parameter
        #of the actuator class
        s.eps = eps_s
        print s.eps, s.eps_0, s.r_1, s.r_2, s.r_1_0, s.r_2_0, s.max_eps, s.min_eps
        print s.min_theta, s.max_theta
        s.calculate_theta()
        s.update()
        
        l.theta = s.theta
        l.update()

        #SMA (Constitutive equation: coupling via sigma)
        data = constitutive_model(T_0, T_final, MVF_init, i, n, eps_s,
                                  eps_t_0, sigma_0, s.eps_0, plot = 'False')
        
        s.sigma = data[0][i][0]
        s.calculate_force(source = 'sigma')
        tau_s = s.calculate_torque()
        
        #Linear (Geometric equation: coupling via theta)
        l.calculate_force()
        tau_l = l.calculate_torque()
        
        #weight (Geometric equation: coupling via theta)
        tau_w = - r_w*W*math.cos(l.theta)
        
        #aerodynamic (Panel method: coupling via theta)
        if aero_loads:
            # The deflection considered for the flap is positivite in
            # the clockwise, contrary to the dynamic system. Hence we need
            # to multiply it by -1.
            Cm = calculate_flap_moment(x, y, alpha, J['x'], - l.theta,
                                       unit_deflection = 'rad')
            tau_a = Cm*q*chord**2
        else:
            tau_a = 0.
            
        print 'tau', tau_s, tau_l, tau_w, tau_a, tau_s + tau_l + tau_w + tau_a
        return tau_s + tau_l + tau_w + tau_a
        
    def deformation_theta(theta = -math.pi/2., plot = False):
        """Return lists of deformation os SMA actuator per theta"""
        
        theta_list = np.linspace(0, theta)
        eps_s_list = []
        eps_l_list = []
        
        for theta in theta_list:
            s.theta = theta
            l.theta = theta
            
            s.update()
            l.update()
            
            l.calculate_force()
            
            eps_s_list.append(s.eps)
            eps_l_list.append(l.eps)
        if plot:
            import matplotlib.pyplot as plt
            plt.figure()    
            plt.plot(np.degrees(theta_list), eps_s_list, 'r', np.degrees(theta_list), eps_l_list, 'b')  
            plt.xlabel('$\\theta (degrees)$')
            plt.ylabel('$\epsilon$')
    
        return eps_s_list, theta_list

    def plot_flap(x, y, x_J, theta):
        """
        Plot flap with actuators. theta is clockwise positive.
        @Author: Endryws (modified by Pedro Leal)
        
        Created on Fri Mar 18 14:26:54 2016
        """
        import matplotlib.pyplot as plt
        plt.figure()
        x_dict, y_dict = af.separate_upper_lower(x, y)
        
        # Below I create the dictionarys used to pass to the function find_hinge
        upper = {'x': x_dict['upper'], 'y': y_dict['upper']} # x and y upper points
        lower = {'x': x_dict['lower'], 'y': y_dict['lower']} # x and y lower points
        hinge = af.find_hinge(x_J, upper, lower) 
        
        #=======================================================================
        # With the Joint (hinge) point, i can use the find flap function to
        # found the points of the flap in the airfoil.
        #=======================================================================
        
        data = {'x': x, 'y': y}
        static_data, flap_data = af.find_flap(data, hinge)
        R = hinge['y_upper']
        theta_list = np.linspace(3*math.pi/2, math.pi/2, 50)
        x_circle_list = hinge['x'] + R*np.cos(theta_list)
        y_circle_list = hinge['y'] + R*np.sin(theta_list)

        n_upper = len(flap_data['x'])/2
        
        # Ploting the flap in the original position
        plt.plot(flap_data['x'][:n_upper],flap_data['y'][:n_upper],'k--')
        plt.plot(flap_data['x'][n_upper:],flap_data['y'][n_upper:],'k--')         
        # Rotate and plot
        upper = {'x': np.concatenate((flap_data['x'][:n_upper], x_circle_list)),
                 'y': np.concatenate((flap_data['y'][:n_upper], y_circle_list))}
        lower = {'x':(flap_data['x'][n_upper:]),
                 'y':(flap_data['y'][n_upper:])}
                
        rotated_upper, rotated_lower = af.rotate(upper, lower, hinge, theta, 
                                                 unit_theta = 'rad')
        plt.plot(static_data['x'], static_data['y'],'k')
        
        plt.plot(rotated_upper['x'], rotated_upper['y'],'k')
        plt.plot(rotated_lower['x'], rotated_lower['y'],'k')
        plt.axes().set_aspect('equal')
        
        s.plot_actuator()
        l.plot_actuator()
#==============================================================================
# Material and flow properties
#==============================================================================
    #SMA properties
    E_M = 60E9
    E_A = 60E9
    sigma_crit = 140e6
    H_max = 0.04
    H_min = 0.
    k = 0.021e-6
    
    Air_props= air_properties(altitude, unit='feet')
    rho = Air_props['Density']
    q = 0.5*rho*V**2

#===========================================================================
# Generate airfoil (NACA0012)
#===========================================================================
    xf.call(airfoil, output='Coordinates')
    filename = xf.file_name(airfoil, output='Coordinates')
    Data = xf.output_reader(filename, output='Coordinates', header = ['x','y'])
    #The coordinates from Xfoil are normalized, hence we have to multiply
    #by the chord
    x = []
    y = []
    for i in range(len(Data['x'])):
        x.append( Data['x'][i]*chord )
        y.append( Data['y'][i]*chord )

#==============================================================================
# Initial conditions   
#==============================================================================
    #Initial transformation strain
    eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - sigma_crit)))
    #Define initial strain
    eps_0 = eps_t_0 + sigma_o/E_M
    #Linear actuator (s)
    l = actuator(linear, J, zero_stress_length = length_l, material = 'linear')
    #Sma actuator (l)
    s = actuator(sma, J, eps_0 = eps_0, material = 'SMA')
    
#===========================================================================
# Although commented this part can be used to calculate sigma_0 for full transformation
#============================================================================
#    def diff_eq_sigma(sigma_o):
#        return  1./E_M + (H_max - H_min)*k*math.exp(-k*(abs(sigma_o) - sigma_crit))
#        
#    def eq_sigma(sigma_o):
#        s.eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - sigma_crit)))
#        return - s.eps_0 + s.eps_t_0 + sigma_o/E_M
#    
#    if s.max_eps < eps_0:
#        s.eps_0 = s.max_eps
#        s.sigma = newton(eq_sigma, sigma_o, diff_eq_sigma)

    #Input initial stress   
    s.sigma = sigma_o
    s.calculate_force(source = 'sigma')
    s.eps_t_0 = eps_t_0
#    print s.eps_0, s.eps_t_0, s.sigma, s.max_theta,

    # TODO: For now K is calculated in function of the rest. Afterwards it will
    # be an input
    if alpha != 0.:
        raise Exception('The initial equilibirum equation only makes sense for alpha equal to zero!!')
    #aerodynamic (Panel method: coupling via theta)
    if aero_loads:
        # The deflection considered for the flap is positivite in
        # the clockwise, contrary to the dynamic system. Hence we need
        # to multiply it by -1.
        
        Cm = calculate_flap_moment(x, y, alpha, J['x'], - l.theta,
                                   unit_deflection = 'rad')
        tau_a = Cm*q*chord**2
    else:
        tau_a = 0.
    
    l.k = ((s.F/s.length_r)*(s.x_p*s.r_2 - s.y_p*s.r_1) + r_w*W + tau_a)/(l.eps*(l.y_p*l.r_1 - l.x_p*l.r_2))
#    print 'stiffness: ', l.k

    l.calculate_force(source = 'strain')
    
    s.theta = l.calculate_theta()
    l.theta = s.theta
    
    #Calculate initial torques
    s.calculate_torque()    
    l.calculate_torque()

#    s.plot_actuator()
#    l.plot_actuator()

    s.find_limits(t = 0.12*chord, c = chord, theta_0 = 0)
    l.find_limits(t = 0.12*chord, c = chord, theta_0 = 0)
#    print s.max_theta, s.min_theta
#    print l.max_theta, l.min_theta
    print 's: limits', s.max_theta, s.min_theta
    print 'l: limits', l.max_theta, l.min_theta
    s.theta = s.max_theta
    s.update()
    l.theta = s.max_theta
    l.update()
    
#    deformation_theta(theta = -math.pi/2., plot = True)
    plot_flap(x, y, J['x'], -s.theta)

    BREAK
##==============================================================================
# Matlab simulation
##==============================================================================         
    eps_s = eps_0
    eps_s_list = [eps_s]
    eps_l_list = [l.eps]
    theta_list = [s.theta]

    #Because of the constraint of the maximum deflection, it is possible that
    #the number of steps is smaller than n
    n_real = 1
    for i in range(1, n):
        eps_s = newton(equilibrium, x0 = eps_s, args = ((s, l, T_0, T_final, 
                       MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
                       q, chord, J['x'], True)), maxiter = 500, 
                       tol = 1.0e-6)

        s.eps = eps_s
        s.calculate_theta()
        if s.theta <= s.max_theta:
            break
        else:
            n_real +=1
        l.theta = s.theta
        l.update()
        eps_s_list.append(eps_s)
        eps_l_list.append(l.eps)
        theta_list.append(math.degrees(s.theta))
        print i, eps_s, eps_0, s.theta
    
    #Extra run with prescribed deformation (which has already been calculated)
    # to get all the properties
    for i in range(1, n_real):
        data = constitutive_model(T_0, T_final, MVF_init, i, n, 
                                  eps_s_list[i], eps_t_0, sigma_0 = sigma_o,
                                  eps_0 = eps_0, plot = 'False')
    delta_xi = 1. - data[1][-1][0]
    
    if all_outputs:
        sigma_list = []
        for i in range(len(data[0])):
            sigma_list.append(data[0][i][0])
        MVF_list = []
        for i in range(len(data[1])):
            MVF_list.append(data[1][i][0])    
        T_list = []
        for i in range(len(data[2])):
            T_list.append(data[2][i][0])    
        eps_t_list = []
        for i in range(len(data[3])):
            eps_t_list.append(data[3][i][0])
        return eps_s_list, theta_list, sigma_list, MVF_list, T_list, eps_t_list
    else:
        T = data[2][n_real][0]
        return delta_xi, s.theta, l.k, T

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
    sigma_o = 400e6
           
    airfoil = "naca0012"
    chord = 0.6175
    t = 0.12*chord

    J = {'x':0.25, 'y':0.}
    
    # need to transform normalized coordiantes in to global coordinates
    sma['y+'] = sma['y+']*thickness(sma['x+'], t, chord)/2.
    sma['y-'] = sma['y-']*thickness(sma['x-'], t, chord)/2.
    
    linear['y+'] = linear['y+']*thickness(linear['x+'], t, chord)/2.
    linear['y-'] =  linear['y-']*thickness(linear['x-'], t, chord)/2.
    
    #Adding the area key to the dictionaries
    sma['area'] = math.pi*0.00025**2
    linear['area'] = 0.001
    
    # Design constants

    #original bias spring length
    length_l = 0.06 #
    

    #arm length to center of gravity
    r_w = 0.15
    
    #Aicraft weight (mass times gravity)
    W = 0.06*9.8
    alpha = 0.
    V = 10 #m/s
    altitude = 10000. #feet
    

    
    ## Temperature
    T_0 = 220.15
    T_final = 400.15
     
    #Initial martensitic volume fraction
    MVF_init = 1.
    
    # Number of steps
    n = 200
    
    delta_xi, theta, k, T = flap(airfoil, chord, J, sma, linear, sigma_o, 
                           length_l, W, r_w, V, altitude, alpha, T_0, 
                           T_final, MVF_init, n, all_outputs = False,
                           import_matlab = import_matlab, eng=eng)
    
    return {'delta_xi': delta_xi, 'theta': theta, 'k': k, 'T':T}
    
if __name__ == '__main__':
    J = {'x':0.25, 'y':0.}
    # Position coordinates from holes. y coordinates are a fraction of thickness/2.
    sma = {'x-': J['x'], 'y-': -0.02*2/0.0441137488474, 'x+': 0.1225 + J['x'],
           'y+': 0.0135*2/0.0345972364185}
    linear = {'x-': J['x'], 'y-': 0.032*2/0.0441137488474, 'x+': 0.146 + J['x'], 
              'y+': -0.0135*2/0.0321083851839}

    #SMA Pre-stress
    sigma_o = 400e6
    data = run({'sma':sma, 'linear':linear, 'sigma_o':sigma_o})
    print 'delta_xi', data['delta_xi'], 'theta: ', data['theta']
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
    
        
