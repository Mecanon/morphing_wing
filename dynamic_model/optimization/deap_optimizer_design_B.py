# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 17:40:21 2016

@author: Pedro Leal
"""

#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.

import array
import random
import numpy as np
import os
import sys
import math
from scipy.interpolate import interp1d
import pickle

from deap import base
from deap import benchmarks
from deap import creator
from deap import tools
import matlab.engine

#adding path to static model
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)

from static_model_B import run_multiobjective
#from power_usage import power

#==============================================================================
# Power calculation
#==============================================================================

def power(delta_t, sigma, T, xi, eps_s, L_s, output = "all"):
    """
    Calculate work, power and current.
    
    - output: defines what is the function output (Power or all)
    """
    sigma_o = 100e6
    r = 0.000381/2.
    d = 2*r
    T_o = 273.15 + 30
    
    alpha = 0.           #set to zero on purpose
    c = 837.36              #invented
    rho = 6450.
    
    #Transformation strain properties
    H_max = 0.0550
    H_min = 0.0387
    sigma_crit = 0
    k = 4.6849e-09
    
    rho_E_M = 0.8e-6         #Dynalloy
    rho_E_A = 1.0e-6         #Dynalloy
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
    CP_list = [1.0038, 1.0049, 1.0063, 1.0082, 1.0106, 1.0135, 1.0206]
    T_list = [275., 300., 325., 350., 375., 400., 450.]
    Cp_f = interp1d(T_list, CP_list)
    # Air conductivity
    k_list = [2.428e-5, 2.624e-5, 2.816e-5, 3.003e-5, 3.186e-5, 3.365e-5, 3.710e-5]
    k_f = interp1d(T_list, k_list)
    
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
        
        T_avg = (T[i] + T[i-1])/2.
        Cp_air = Cp_f(T_avg)
        k_air = k_f(T_avg)
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
        P = math.pi*r**2*L_s[i]*((T[i]*alpha*delta_sigma + \
            rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) )/delta_t + \
            2.*(h/r)*(T[i] - T_o))
        
        P_list.append(P)
        
        if output == 'all':
            I = r*math.pi*math.sqrt((r/rho_E)*((r/delta_t)*((T[i]*alpha*delta_sigma + \
                rho*c*delta_T + delta_xi*(-pi_t + rho_delta_s0*T[i]) ) + \
                2.*h*(T[i] - T_o))))
        
            
            dW = math.pi*r**2*L_s[0]*0.5*(sigma[i]+sigma[i-1])*delta_eps
            
            I_list.append(I)
            W_list.append(dW)
        
    Total_power = 0
    for i in range(len(P_list)-1):
        Total_power += delta_t*(P_list[i] + P_list[i+1])/2.
    if output == 'all':
        return I_list, P_list, W_list, Total_power
    elif output == "power":
        return Total_power
        
#==============================================================================
# Objective function
#==============================================================================
def objfunc(x):
    x_J = .75
    length_steel = 0.05

    #SMA Pre-stress
    sigma_o = 100e6
    
    inputs = {'sma':{'x-': x_J - length_steel - x[0], 'y-':-x[2],
                     'x+': x_J - length_steel, 'y+':-x[2],
                     'pulley_position':'down'},
              'linear':{'x-':x_J - length_steel - x[1], 'y-':x[2],
                        'x+':x_J - length_steel, 'y+':x[2],
                       'actuator_type': 'wire',
                       'pulley_position':'up'},
              'sigma_o':sigma_o, 'R':x[2], 'T_f': x[3]}

    DataFile = open('opt_data.txt','a')
    for x_i in x:
        DataFile.write( '\t %.5f' % (x_i) )
    DataFile.close()
            
    theta, sigma, T, MVF, eps_s, L_s = run_multiobjective(inputs = inputs, parameters = [eng])
    
    theta = theta[-1]
    
    delta_t = 0.05
    
    P= power(delta_t, sigma, T, MVF, eps_s, L_s, output = "power")
    
    DataFile = open('opt_data.txt','a')
    DataFile.write( '\t %.5f \t %.5f' % (theta, P) )
    DataFile.write('\n')
    DataFile.close()
    
    return theta, P
#==============================================================================
# Start Matlab engine
#==============================================================================
eng = matlab.engine.start_matlab()
#Go to directory where matlab file is
eng.cd('..')
eng.cd('SMA_temperature_strain_driven')

#==============================================================================
# DEAP algorithm   
#==============================================================================
creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0))
creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)

toolbox = base.Toolbox()

chord = 1.
x_hinge = 0.75
safety = 0.005*chord
# Problem definition
BOUND_LOW = [0.1, 0.1, 0.001, 273.15+30.]

BOUND_UP = [0.6, 0.6, 0.03, 273.15+140.]

NDIM = 9

def uniform(low, up, size=None):
    try:
        return [random.uniform(a, b) for a, b in zip(low, up)]
    except TypeError:
        return [random.uniform(a, b) for a, b in zip([low] * size, [up] * size)]

toolbox.register("attr_float", uniform, BOUND_LOW, BOUND_UP, NDIM)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", objfunc)
toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0)
toolbox.register("mutate", tools.mutPolynomialBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0, indpb=1.0/NDIM)
toolbox.register("select", tools.selNSGA2)

def main(seed=None):
    random.seed(seed)

    # Number of generations
    NGEN = 50
    # Population size (has to be a multiple of 4)
    MU = 40
    # Mating probability
    CXPB = 0.9

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    # stats.register("avg", np.mean, axis=0)
    # stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    
    logbook = tools.Logbook()
    logbook.header = "gen", "evals", "std", "min", "avg", "max"
    
    pop = toolbox.population(n=MU)

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in pop if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    # This is just to assign the crowding distance to the individuals
    # no actual selection is done
    pop = toolbox.select(pop, len(pop))
    
    record = stats.compile(pop)
    logbook.record(gen=0, evals=len(invalid_ind), **record)
    print(logbook.stream)

    # Begin the generational process
    for gen in range(1, NGEN):
        # Vary the population
        offspring = tools.selTournamentDCD(pop, len(pop))
        offspring = [toolbox.clone(ind) for ind in offspring]
        
        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
            if random.random() <= CXPB:
                toolbox.mate(ind1, ind2)
            
            toolbox.mutate(ind1)
            toolbox.mutate(ind2)
            del ind1.fitness.values, ind2.fitness.values
        
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Select the next generation population
        pop = toolbox.select(pop + offspring, MU)
        record = stats.compile(pop)
        logbook.record(gen=gen, evals=len(invalid_ind), **record)
        print(logbook.stream)

#    print("Final population hypervolume is %f" % hypervolume(pop, [11.0, 11.0]))

    return pop, logbook
        
if __name__ == "__main__":

    DataFile = open('opt_data.txt','w')
    key_list = ['l_s', 'l_l', 'R', 'T_f']
    output_list = ['theta', 'power']
    for key in key_list + output_list:
        DataFile.write(key + '\t')
    DataFile.write('\n')
    DataFile.close()
   
    pop, stats = main()
    pop.sort(key=lambda x: x.fitness.values)
    
    print(stats)
    
    import matplotlib.pyplot as plt
    
    front = np.array([ind.fitness.values for ind in pop])
    
    pickle.dump( front, open( "front.p", "wb" ) )
    pickle.dump( pop, open( "pop.p", "wb" ) )
    pickle.dump( stats, open( "stats.p", "wb" ) )
    plt.scatter(np.rad2deg(front[:,0]), front[:,1], c="b")
    plt.axis("tight")
    plt.grid()
    plt.xlabel("Deflection angle (${}^{\circ}$)")
    plt.ylabel("Heating load (J)")
    plt.show()