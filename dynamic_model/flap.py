# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 17:27:40 2016

@author: Pedro Leal
"""
import math
from scipy.optimize import brentq, minimize_scalar, fixed_point
import numpy as np
    
from AeroPy import calculate_flap_moment
from aero_module import air_properties
import airfoil_module as af
import xfoil_module as xf

from actuator import actuator

def flap(airfoil, chord, J, sma, linear, sigma_o, W, r_w, V,
         altitude, alpha, T_0, T_final, MVF_init, n, all_outputs = False,
         import_matlab = True, eng = None, aero_loads = True, 
         cycling = 'False', n_cycles = 0):
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
    
    from scipy.interpolate import interp1d
    import pickle
    import os.path
 
    if import_matlab and eng == None:
        import matlab.engine
        #Start Matlab engine
        eng = matlab.engine.start_matlab()
        #Go to directory where matlab file is
        if import_matlab:
            eng.cd(' ..')
            eng.cd('SMA_temperature_strain_driven')
        else:
            eng.cd('SMA_temperature_strain_driven')
            
    def constitutive_model(T, MVF_init, i, n, eps, eps_t_0, sigma_0 = 0,
            eps_0 = 0, plot = 'True'):
        """Run SMA model
        
        - all inputs are scalars"""
        k = i+1
        if k == n:
            data = eng.OneD_SMA_Model(k, eps, T, MVF_init, 
                                      eps_t_0, sigma_0, eps_0, n, plot,
                                      nargout=6) 
        else:
            data = eng.OneD_SMA_Model(k, eps, T, MVF_init,
                                      eps_t_0, sigma_0, eps_0, n, 'False',
                                      nargout=6)
        return data

    def equilibrium(eps_s, s, l, T, MVF_init, sigma_0,
                   i, n, r_w, W, x = None, y = None, alpha = 0.,
                   q = 1., chord = 1., x_hinge = 0.25, aero_loads = aero_loads,
                   return_abs = False):
        """Calculates the moment equilibrium. Function used for the 
        secant method.
        """
        #calculate new theta for eps_s and update all the parameter
        #of the actuator class
        s.eps = eps_s
#        print s.eps, s.eps_0, s.r_1, s.r_2, s.r_1_0, s.r_2_0, s.max_eps, s.min_eps
#        print s.min_theta, s.max_theta
#        try:
#        print eps_s
        s.calculate_theta(theta_0 = s.theta)
#        except:
#            s.theta = max(s.max_theta, l.max_theta)
#            s.update()
#            l.theta = max_theta
#            l.update()
#            plot_flap(x, y, J['x'], -s.theta)
#            raise Exception("Inverse solution for theta did not converge")
        s.update()
        
        l.theta = s.theta
        l.update()
        #SMA (Constitutive equation: coupling via sigma)
        data = constitutive_model(T, MVF_init, i, n, eps_s,
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
#            print alpha, x_hinge, l.theta
            if abs(l.theta) > math.pi/2.:
                tau_a = - np.sign(l.theta)*4
            else:
                Cm = Cm_function(l.theta)
                tau_a = Cm*q*chord**2
#            Cm = calculate_flap_moment(x, y, alpha, x_hinge, - l.theta,
#                                       unit_deflection = 'rad')
        else:
            tau_a = 0.
            
#        print 'tau', tau_s, tau_l, tau_w, tau_a, tau_s + tau_l + tau_w + tau_a
        f = open('data', 'a')
        f.write('\t Inner loop \t'+ str(T[i]) + '\t' + str( eps_s) + '\t' + \
                str( tau_s + tau_l + tau_w + tau_a) + '\t' + str( tau_s) + '\t' + str(tau_l) + '\t' + str(tau_w) + '\t' + str(tau_a) + '\t' + \
                str(l.theta)  + '\n')
        f.write('\t ' + str(i) + '\t' + str(n) +'\n')
        f.write('\t ' + str(l.F) + '\t' + str(s.F) + '\n')
        f.close()
        if return_abs:
            return abs(tau_s + tau_l + tau_w + tau_a)
        else:
            return tau_s + tau_l + tau_w + tau_a

    def eps_s_fixed_point(eps_s, s, l, T, MVF_init, sigma_0,
                   i, n, r_w, W, x = None, y = None, alpha = 0.,
                   q = 1., chord = 1., x_hinge = 0.25, aero_loads = aero_loads,
                   return_abs = False):
        """Calculates the SMA strain. Function used for the 
        fixed position iteration method. Restriction d epsilon_s / dt < 1.
        """
        #calculate new theta for eps_s and update all the parameter
        #of the actuator class
        s.eps = eps_s
#        print s.eps, s.eps_0, s.r_1, s.r_2, s.r_1_0, s.r_2_0, s.max_eps, s.min_eps
#        print s.min_theta, s.max_theta
#        try:
#        print eps_s
        s.calculate_theta(theta_0 = s.theta)
#        except:
#            s.theta = max(s.max_theta, l.max_theta)
#            s.update()
#            l.theta = max_theta
#            l.update()
#            plot_flap(x, y, J['x'], -s.theta)
#            raise Exception("Inverse solution for theta did not converge")
        s.update()
        
        l.theta = s.theta
        l.update()
        #SMA (Constitutive equation: coupling via sigma)
        data = constitutive_model(T, MVF_init, i, n, eps_s,
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
#            print alpha, x_hinge, l.theta
            if abs(l.theta) > math.pi/2.:
                tau_a = - np.sign(l.theta)*4
            else:
                Cm = Cm_function(l.theta)
                tau_a = Cm*q*chord**2
#            Cm = calculate_flap_moment(x, y, alpha, x_hinge, - l.theta,
#                                       unit_deflection = 'rad')
        else:
            tau_a = 0.
        
        B = tau_s/s.sigma
        eps_t = data[3][i][0] 
        E = data[5][i][0] 
        eps_s = eps_t - (tau_l + tau_w + tau_a)/(B*E)
        
#        print 'tau', tau_s, tau_l, tau_w, tau_a, tau_s + tau_l + tau_w + tau_a
        f = open('data', 'a')
        f.write('\t Inner loop \t'+ str(T[i]) + '\t' + str( eps_s) + '\t' + \
                str( tau_s + tau_l + tau_w + tau_a) + '\t' + str( tau_s) + '\t' + str(tau_l) + '\t' + str(tau_w) + '\t' + str(tau_a) + '\t' + \
                str(l.theta)  + '\n')
        f.write('\t ' + str(i) + '\t' + str(n) +'\n')
        f.write('\t ' + str(l.F) + '\t' + str(s.F) + '\n')
        f.close()        
        return eps_s
            
    def deformation_theta(theta = -math.pi/2., plot = False):
        """Return lists of deformation os SMA actuator per theta"""
        
        theta_list = np.linspace(-theta, theta)
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

    def plot_flap(x, y, x_J, y_J = None, theta = 0):
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
        
        if y_J != None:
            for i in range(y_J):
                plt.scatter(x_J , y_J[i])
        plt.grid()
        plt.locator_params(axis = 'y', nbins=6)
        plt.xlabel('${}_{I}x + x_J$')
        plt.ylabel('${}_{I}y$')
        border = 0.05
        plt.xlim(-border, 1+border)
        plt.savefig( str(np.floor(100*abs(theta))) + "_configuration.png")
        plt.close()

#==============================================================================
# Material and flow properties
#==============================================================================
    #SMA properties
    E_M = 60E9
    E_A = 60E9
    sigma_crit = 140e6
    H_max = 0.047
    H_min = 0.
    k = 0.021e-6
    
    Air_props= air_properties(altitude, unit='feet')
    rho = Air_props['Density']
    q = 0.5*rho*V**2

    #Spring properties(Squared or closed)
    C = 10. #Spring index
    Nt = 17. #Total number of springs
    safety_factor = 1.2
    A = 2211e6*10**(0.435) #Mpa.mm**0.125
    
#===========================================================================
# Generate Temperature arrays
#===========================================================================
    if n_cycles == 0:
        T = np.linspace(T_0, T_final, n)
    else:
        T = np.append(np.linspace(T_0, T_final, n),
                      np.linspace(T_final, T_0, n))
        for i in range(n_cycles-1):
            T = np.append(T,T)
    T = list(T)
#===========================================================================
# Generate airfoil (NACA0012)
#===========================================================================
    xf.call(airfoil, output='Coordinates')
    filename = xf.file_name(airfoil, output='Coordinates')
    Data = xf.output_reader(filename, output='Coordinates',
                            header = ['x','y'])
    #The coordinates from Xfoil are normalized, hence we have to multiply
    #by the chord
    x = []
    y = []
    for i in range(len(Data['x'])):
        x.append( Data['x'][i]*chord )
        y.append( Data['y'][i]*chord )
    
    filename = airfoil + '_' + str(int(100*J['x'])) + '_' + \
               str(int(100*chord)) + '.p'
    
    #Generate aerodynamic moment data. If already exists, load it
    if os.path.isfile(filename):
        with open(filename, "rb") as f:
            Cm_list = pickle.load( f )
            theta_list = pickle.load( f )
    else:
        theta_list = np.linspace(-math.pi/2., math.pi/2., 100)
        Cm_list = []
        for theta in theta_list:
            Cm = calculate_flap_moment(x, y, alpha, J['x'], - theta,
                                       unit_deflection = 'rad')
            Cm_list.append(Cm)
        with open(filename, "wb") as f:
            pickle.dump( Cm_list, f )
            pickle.dump( theta_list, f )
        
    Cm_function = interp1d(theta_list, Cm_list)
#===========================================================================
# SMA actuator  
#===========================================================================
    #Initial transformation strain
    eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(sigma_o) - \
                                                     sigma_crit)))
    #Define initial strain
    eps_0 = eps_t_0 + sigma_o/E_M
    #Sma actuator (s)
    s = actuator(sma, J, eps_0 = eps_0, material = 'SMA')

    #Check if crossing joint. If True do nothing
    if s.check_crossing_joint(tol = 0.01):
        if all_outputs:
            return 0., s.theta, 0., 200.
        else:
            return s.theta
    else:
        #Input initial stress   
        s.sigma = sigma_o
        s.calculate_force(source = 'sigma')
        s.eps_t_0 = eps_t_0
    #    print s.eps_0, s.eps_t_0, s.sigma, s.max_theta,
    
        if alpha != 0.:
            raise Exception('The initial equilibirum equation only ' + \
                            'makes sense for alpha equal to zero!!')
    
        s.update()     
        #Calculate initial torques
        s.calculate_torque()   
        
        tau_w = - r_w*W 

        if aero_loads:
            # The deflection considered for the flap is positivite in
            # the clockwise, contrary to the dynamic system. Hence we need
            # to multiply it by -1.
            Cm = Cm_function(s.theta)
            tau_a = Cm*q*chord**2
        else:
            tau_a = 0.

#===========================================================================
# Linear actuator
#===========================================================================
        #Linear actuator (l)
        l = actuator(linear, J, material = 'linear')
        l.theta = s.theta
        
        #Check if crossing joint. If True do nothing
        if l.check_crossing_joint(tol = 0.01):
            if all_outputs:
                return 0., s.theta, 0., 200.
            else:
                return s.theta
        else:        
            #l.k = - (s.torque + tau_w + tau_a)/(l.eps*(l.y_p*l.r_1 - l.x_p*l.r_2))
            l.F = - l.length_r_0*(s.torque + tau_w + \
                    tau_a)/(l.y_p*l.r_1 - l.x_p*l.r_2)
            #Spring components
            if l.F < 0:
                l.d = (safety_factor*(2*C+1)*abs(1.5*l.F)/(0.1125*A*math.pi))**(1./1.855)
            else:
                l.d = (safety_factor*(C*(4*C+2)/(4*C-3))*abs(1.5*l.F)/(0.05625*A*math.pi))**(1./1.855)

            l.D = C*l.d
            
            if l.d < 0.000838:
                G = 82.7e9
                E = 203.4e9
            elif l.d < 0.0016:
                G = 81.7e9
                E = 200.0e9
            elif l.d < 0.00318:
                G = 81.0e9
                E = 196.6e9
            else:
                G = 80.0e9
                E = 193.0e9

            if l.F < 0:
                Na = Nt - 2 #Active number of springs
            else:
                Nt = (E*G*l.d*l.length_r_0 + (1.-2.*C)*E*G*l.d**2 - \
                     8.*l.F*G*C**3)/(E*G*l.d**2 + 8.*E*l.F*C**3)
                Na = Nt + G/E #Equivalent active number of springs
                print "Ns", Nt, Na
                
            l.k = l.d**4*G/(8*l.D**3*Na)
            if l.F < 0.:
                l.zero_stress_length = -l.F/l.k + l.length_r_0
                l.solid = (Nt + 1)*l.d
            else:
                l.zero_stress_length = (2*C - 1 + Nt)*l.d
                print "different zero lengths", (2*C - 1 + Nt)*l.d, -l.F/l.k + l.length_r_0
                l.solid = l.zero_stress_length
            print 'before'    
            print l.F,  l.k, l.d, G, l.D, Na, Nt, l.zero_stress_length, safety_factor, l.length_r_0, l.length_r_0/l.zero_stress_length
            print 'strain', l.length_r_0/l.zero_stress_length -1.         
            l.update()
            l.eps_0 = l.eps
            l.calculate_force(source = 'strain')
            l.calculate_torque()
            l.calculate_theta()
            print "after"
            print l.F,  l.k, l.d, G, l.D, Na, Nt, l.zero_stress_length, safety_factor, l.length_r_0
            print 'strain', l.eps          
            print "torques"            
            print l.torque, l.theta,  s.torque + tau_w + tau_a, s.torque, tau_w, tau_a, l.torque
            t = 0.12
            
            y_J = af.Naca00XX(chord, t, [J['x']], return_dict = 'y')
            
            s.find_limits(y_J, theta_0 = 0)
            l.find_limits(y_J, theta_0 = 0)
    
            print 's: limits', s.max_theta, s.min_theta
            print s.max_theta_A, s.max_theta_B
            print 'l: limits', l.max_theta, l.min_theta
            print l.max_theta_A, l.max_theta_B
            
            max_theta = max(s.max_theta, l.max_theta)
            
    #        #The following code is good for plotting the airfoil and etc
    #        s.theta = s.max_theta_B
    #        s.update()
    #        l.theta = s.max_theta_B
    #        l.update()
    #    #    
    ##        deformation_theta(theta = -math.pi/2., plot = True)
    #        plot_flap(x, y, J['x'], theta = -s.theta)
    #        plt.scatter(J['x'] , y_J['l'])
    #        plt.scatter(J['x'] , y_J['u'])
        ##==============================================================================
        # Matlab simulation
        ##==============================================================================         
            eps_s = eps_0
            eps_s_list = [eps_s]
            eps_l_list = [l.eps]
            theta_list = [s.theta]
            F_l_list =[l.calculate_force()]
            L_s_list = [s.length_r]
            #Because of the constraint of the maximum deflection, it is possible that
            #the number of steps is smaller than n
            n_real = 1
            
            #Create new data file and erase everything inside
            f = open('data', 'w')
            f.close()
            
            if n_cycles == 0:
                iterations = n
            else:
                iterations = 2*n*n_cycles
            a_counter = 0
            #original number of steps without any restrictions
            n_o = n

            eps_s_prev = [eps_0, eps_0]
            for i in range(1, iterations):
                #because n_real can be different of n, it is posssible that
                #i is greater than the actual number of iterations. Need
                #to break
                if i == iterations:
                    break
#                equilibrium_0 = equilibrium(eps_s, s, l, T, MVF_init, sigma_o,
#                           i, n, r_w, W, x, y, alpha, q, chord, J['x'], True)
#                f = open('data', 'a')
#                f.write(str(equilibrium_0) + '\t' + str(eps_0) + '\n')
#                f.close()
                print "T length: ", len(T)
#                if not abs(equilibrium_0) < 1e-8:
                delta_eps = eps_s_prev[1] - eps_s_prev[0]
                eps_s = fixed_point(eps_s_fixed_point, eps_s + delta_eps, args=((s, l, T, 
                                       MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
                                       q, chord, J['x'], True)), xtol = 1e-8)
                eps_s_prev[0] = eps_s_prev[1]
                eps_s_prev[1] = eps_s
#                    if 1.1*eps_s> eps_0:
#                        try:
#                            eps_s = brentq(equilibrium, eps_s*0.9, eps_0, args=((s, l, T, 
#                                           MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                           q, chord, J['x'], True)), rtol = 1e-8)
#                        except:
#                            eps_before = (0.999)*eps_s
#                            for j in range(10):
#                                try:
#                                    eps_s = brentq(equilibrium, eps_before, eps_0, args=((s, l, T, 
#                                                   MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                                   q, chord, J['x'], True)), rtol = 1e-8)
#                                    break
#                                except:
#                                    eps_before = eps_before - 0.001*eps_s
#    
#                            if j == 9:
#                                f = open('data', 'a')
#                                f.write('bounded1')
#                                f.close()
#                                OptimizeResult = minimize_scalar(equilibrium, (0.9*eps_0, eps_0), bounds = (0.9*eps_0, eps_0), args=((s, l, T, 
#                                               MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                               q, chord, J['x'], True, True)), method = 'bounded', options = {'xatol' : 1e-08})
#                                eps_s = OptimizeResult.x
#    #                            print OptimizeResult.fun
#    
#                    else:
#                        if n != n_real:
#                            try:
#                                eps_s = brentq(equilibrium, eps_s*0.9, 1.1*eps_s, args=((s, l, T, 
#                                               MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                               q, chord, J['x'], True)), rtol = 1e-8)
#                            except:
#                                eps_before = (0.999)*eps_s
#                                for j in range(10):
#                                    try:
#                                        eps_s = brentq(equilibrium, eps_before, eps_s, args=((s, l, T, 
#                                                       MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                                       q, chord, J['x'], True)), rtol = 1e-8)
#                                        break
#                                    except:
#                                        eps_before = eps_before - 0.002*eps_s
#         
#                                if j==9:
#                                    f = open('data', 'a')
#                                    f.write('bounded2')
#                                    f.close()
#                                    OptimizeResult = minimize_scalar(equilibrium, (0.5*eps_s, 1.1*eps_s), bounds = (0.*eps_s, 1.1*eps_s), args=((s, l, T, 
#                                                   MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                                   q, chord, J['x'], True, True)), method = 'bounded', options = {'xatol' : 1e-08})
#                                    eps_s = OptimizeResult.x
#                        else:
#                            try:
#                                eps_s = brentq(equilibrium, eps_s*0.9, 1.2*eps_s, args=((s, l, T, 
#                                               MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                               q, chord, J['x'], True)), rtol = 1e-8)
#                            except:
#                                eps_after = (1.001)*eps_s
#                                for j in range(10):
#                                    try:
#                                        eps_s = brentq(equilibrium, eps_s, eps_after, args=((s, l, T, 
#                                                       MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                                       q, chord, J['x'], True)), rtol = 1e-8)
#                                        break
#                                    except:
#                                        eps_before = eps_before + 0.002*eps_s
#         
#                                if j==9:
#                                    f = open('data', 'a')
#                                    f.write('bounded2')
#                                    f.close()
#                                    OptimizeResult = minimize_scalar(equilibrium, (0.9*eps_s, eps_0), bounds = (0.*eps_s, eps_0), args=((s, l, T, 
#                                                   MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                                                   q, chord, J['x'], True, True)), method = 'bounded', options = {'xatol' : 1e-08})
#                                    eps_s = OptimizeResult.x                                
#        #                        eps_s = newton(equilibrium, x0 = eps_s, args = ((s, l, T_0, T_final, 
#        #                                       MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#        #                                       q, chord, J['x'], True)), maxiter = 500, 
#        #                                       tol = 1.0e-6)
                s.eps = eps_s
                s.calculate_theta(theta_0 = s.theta)
                s.update()
                
                f = open('data', 'a')
                f.write('Outer loop \t'+ str(i)  + '\t'+ str(eps_s) + '\t'+ str(s.theta)+ '\n')
                f.close()
                
                #stop if actuator crosses joints, exceds maximum theta and theta is positively increasing
                if s.theta <= max_theta or s.check_crossing_joint(tol = 0.001) or l.check_crossing_joint(tol = 0.001) or s.theta > 0.01 or l.length_r <= l.solid:
                    if n_cycles == 0:
                        iterations = n_real
                        T = T[:n_real]
                        break
                    else:
                        a_counter += 1
                        n = n_real #+ 10
                        iterations = 2*n_real*n_cycles
                        T_cycle = T[:n_real]  + T[:n_real][::-1] #+ 10*[T[n_real-1]]
                        T = T_cycle
                        for k in range(n_cycles-1):
                            T += T_cycle
                        print 'iterations: ', iterations, n_real, n_cycles, a_counter, len(T)
                else:
                    print 'n_o: ', n_o, ', n: ', n, 'n_real', n_real
                    if n != n_real:
                        n_real +=1
                    if n_o == n_real:
                        #WRITE APPEND T HERE
                        n = n_real #+ 10
                        if n_cycles == 0:
                            iterations = n_real
                        else:
                            iterations = 2*n_real*n_cycles
                            T_cycle = T[:n_real]  + T[:n_real][::-1] #+ 10*[T[n_real-1]]
                            T = T_cycle
                            for k in range(n_cycles-1):
                                T += T_cycle
                        print 'length of T: ', len(T)
                    l.theta = s.theta
                    l.update()
                    
                    eps_s_list.append(eps_s)
                    eps_l_list.append(l.eps)
                    theta_list.append(math.degrees(s.theta))
                    F_l_list.append(l.calculate_force())
                    L_s_list.append(s.length_r)
        #        print i, eps_s, eps_0, s.theta
    
            plot_flap(x, y, J['x'], theta= -s.theta)
            
            if all_outputs:
                if n_real == 1:
                    return 0., s.theta, 0., 200.
                else:
                    #Extra run with prescribed deformation (which has already been calculated)
                    # to get all the properties]
                    #TODO: include back n_real
                    for i in range(1, iterations):
                        data = constitutive_model(T, MVF_init, i, iterations, 
                                                  eps_s_list[i], eps_t_0,
                                                  sigma_0 = sigma_o,
                                                  eps_0 = eps_0, 
                                                  plot = 'False')
                    delta_xi = 1. - data[1][-1][0]
                    
                    sigma_list = []
                    for i in range(len(data[0])):
                        sigma_list.append(data[0][i][0])
                    MVF_list = []
                    for i in range(len(data[1])):
                        MVF_list.append(data[1][i][0])    
                    eps_t_list = []
                    for i in range(len(data[3])):
                        eps_t_list.append(data[3][i][0])
                    return eps_s_list, eps_l_list, theta_list, sigma_list, MVF_list, T, eps_t_list, theta_list, F_l_list, l.k, L_s_list
            else:
                print "k", l.k
                return s.theta, l.k