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

def flap(airfoil, chord, J, sma, linear, spring, W, r_w, V,
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
#        print 'eps_s: ', eps_s
        s.calculate_theta(theta_0 = s.theta)
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
        tau_a = 0.
        
        B = tau_s/s.sigma
        eps_t = data[3][i][0] 
        E = data[5][i][0] 
        

        eps_s = eps_t - (tau_l + tau_w + tau_a)/(B*E)
        
#        print 'tau', tau_s, tau_l, tau_w, tau_a, tau_s + tau_l + tau_w + tau_a
        f = open('data', 'a')
        f.write('\t Inner loop 2\t'+ str(T[i]) + '\t' + str( eps_s) + '\t' + \
                str( tau_s + tau_l + tau_w + tau_a) + '\t' + str( tau_s) + '\t' + str(tau_l) + '\t' + str(tau_w) + '\t' + str(tau_a) + '\t' + \
                str(l.theta)  + '\n')
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
        
        l.plot_actuator()
        s.plot_actuator()
        
        
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
    E_M = 8.8888e+10
    E_A = 3.7427e+10
    sigma_crit = 0
    H_max = 0.0550
    H_min = 0.0387
    k = 4.6849e-09
    
    Air_props= air_properties(altitude, unit='feet')
    rho = Air_props['Density']
    q = 0.5*rho*V**2
    
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
#    if os.path.isfile(filename):
#        with open(filename, "rb") as f:
#            Cm_list = pickle.load( f )
#            theta_list = pickle.load( f )
#    else:
#        theta_list = np.linspace(-math.pi/2., math.pi/2., 100)
#        Cm_list = []
#        for theta in theta_list:
#            Cm = calculate_flap_moment(x, y, alpha, J['x'], - theta,
#                                       unit_deflection = 'rad')
#            Cm_list.append(Cm)
#        with open(filename, "wb") as f:
#            pickle.dump( Cm_list, f )
#            pickle.dump( theta_list, f )
#        
#    Cm_function = interp1d(theta_list, Cm_list)

#===========================================================================
# Linear actuator
#===========================================================================
    #Linear actuator (l)
    l = actuator(linear, J, material = 'linear')
    l.k = spring['k']
    #Check if crossing joint. If True do nothing
    if l.check_crossing_joint(tol = 0.005):
        print "problem with linear"
        plot_flap(x, y, J['x'], theta= 0)
        if all_outputs:
            return 0., 0, 0., 200.
        else:
            return 0, 9999.
    else:        

        l.zero_stress_length = spring['L_0']
        l.solid = l.zero_stress_length
  
        l.update()
        l.eps_0 = l.eps
        l.calculate_force(source = 'strain')
        l.calculate_torque()
        l.calculate_theta()
#===========================================================================
# SMA actuator  
#===========================================================================
        tau_w = - r_w*W 
        
        #Sma actuator (s)
        s = actuator(sma, J, material = 'SMA')
        
        #Initial force
        s.F = - s.length_r_0*(l.torque + tau_w)/(s.y_p*s.r_1 - s.x_p*s.r_2)
        
        #Initial stress
        s.sigma = s.F/s.area
        sigma_o = s.sigma
        
        #Initial transformation strain
        eps_t_0 = H_min + (H_max - H_min)*(1. - math.exp(-k*(abs(s.sigma) - \
                                                         sigma_crit)))
        s.eps_t_0 = eps_t_0
        
        #Define initial strain
        eps_0 = s.eps_t_0 + s.sigma/E_M
        s.eps_0 = eps_0
        s.zero_stress_length = s.length_r_0/(1. + s.eps_0)
        s.update()
        
        #Calculate initial torques
        s.calculate_torque()  
        print s.eps, s.F, l.F, s.torque, l.torque, tau_w, s.sigma, sigma_o
        #Check if crossing joint. If True do nothing
        if s.check_crossing_joint(tol = 0.005):
            print "problem with SMA"
            plot_flap(x, y, J['x'], theta= 0)
            if all_outputs:
                return 0., s.theta, 0., 200.
            else:
                return s.theta, 9999.
        else:
    
            t = 0.12
            
            y_J = af.Naca00XX(chord, t, [J['x']], return_dict = 'y')
            
            s.find_limits(y_J, theta_0 = 0)
            l.find_limits(y_J, theta_0 = 0)
    
#            print 's: limits', s.max_theta, s.min_theta
#            print s.max_theta_A, s.max_theta_B
#            print 'l: limits', l.max_theta, l.min_theta
#            print l.max_theta_A, l.max_theta_B
#            plot_flap(x, y, J['x'], theta= 0.)            
            max_theta = max(s.max_theta, l.max_theta)
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
            
#            plot_flap(x, y, J['x'], theta= -s.theta)
            
            #Create new data file and erase everything inside
            f = open('data', 'w')
            f.close()
            
            if n_cycles == 0:
                iterations = n
            else:
                iterations = 2*n*n_cycles
            #original number of steps without any restrictions
            n_o = n
            print s.eps, s.F, l.F, s.torque, l.torque, tau_w
            eps_s_prev = [eps_0, eps_0]
            for i in range(1, iterations):
                #because n_real can be different of n, it is posssible that
                #i is greater than the actual number of iterations. Need
                #to break
                if i == iterations:
                    break
                
                #Linear prediction 
                delta_eps = eps_s_prev[1] - eps_s_prev[0]
                #Correction with fixed point iteration
                try:
                    eps_s = fixed_point(eps_s_fixed_point, eps_s + delta_eps, args=((s, l, T, 
                                           MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
                                           q, chord, J['x'], True)), xtol = 1e-8)
                except:
#                    eps_s = fixed_point(eps_s_fixed_point, eps_s, args=((s, l, T, 
#                           MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
#                           q, chord, J['x'], True)), xtol = 1e-8)
                    eps_s_start = eps_s + delta_eps
                    for j in range(10):
                        try:
                            eps_s = fixed_point(eps_s_fixed_point, eps_s_start, args=((s, l, T, 
                                                   MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
                                                   q, chord, J['x'], True)), xtol = 1e-8)
                            break
                        except:
                            eps_s_start  = eps_s_start*float(j)/10.
                eps_s_prev[0] = eps_s_prev[1]
                eps_s_prev[1] = eps_s

                s.eps = eps_s
                s.calculate_theta(theta_0 = s.theta)
                s.update()
                
                f = open('data', 'a')
                f.write('Outer loop \t'+ str(i)  + '\t'+ str(eps_s) + '\t'+ str(s.theta)+ '\n')
                f.close()
#                print 'new theta: ', s.theta
                #stop if actuator crosses joints, exceds maximum theta and theta is positively increasing
                if s.theta <= max_theta or s.check_crossing_joint(tol = 0.001) or l.check_crossing_joint(tol = 0.001) or s.theta > 0.01 or l.length_r <= l.solid:
                    if n_cycles == 0:
                        iterations = n_real
                        T = T[:n_real]
                        break
                    else:
                        n = n_real #+ 10
                        iterations = 2*n_real*n_cycles
                        T_cycle = T[:n_real]  + T[:n_real][::-1] #+ 10*[T[n_real-1]]
                        T = T_cycle
                        for k in range(n_cycles-1):
                            T += T_cycle
                        #Correction with fixed point iteration
                            eps_s = fixed_point(eps_s_fixed_point, eps_s_prev[1], args=((s, l, T, 
                                                   MVF_init, sigma_o, i, n, r_w, W, x, y, alpha, 
                                                   q, chord, J['x'], True)), xtol = 1e-8)
                        eps_s_prev[0] = eps_s_prev[1]
                        eps_s_prev[1] = eps_s
                        
                        s.eps = eps_s
                        s.calculate_theta(theta_0 = s.theta)
                        s.update()
                        
                else:
                    if n != n_real:
                        n_real +=1
                    if n_o == n_real:
                        n = n_real #+ 10
                        if n_cycles == 0:
                            iterations = n_real
                        else:
                            iterations = 2*n_real*n_cycles
                            T_cycle = T[:n_real]  + T[:n_real][::-1] #+ 10*[T[n_real-1]]
                            T = T_cycle
                            for k in range(n_cycles-1):
                                T += T_cycle
                l.theta = s.theta
                l.update()
                
                eps_s_list.append(eps_s)
                eps_l_list.append(l.eps)
                theta_list.append(s.theta)
                F_l_list.append(l.calculate_force())
                L_s_list.append(s.length_r)
    
#            print 'final theta: ', s.theta
            if all_outputs:
                s.theta = theta_list[-1]
                l.theta = s.theta
                s.update()
                l.update()
                plot_flap(x, y, J['x'], theta= -s.theta)
                
                s.theta = 0.
                l.theta = s.theta
                s.update()
                l.update()
                plot_flap(x, y, J['x'], theta= -s.theta)
                
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
#                print theta_list
                return theta_list[-1], l.k
