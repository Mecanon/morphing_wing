# -*- coding: utf-8 -*-
"""
Created on Wed Jun 08 15:06:17 2016

@author: Pedro Leal
"""
import matplotlib.pyplot as plt

from scipy import polyfit, polyval
import numpy as np

from xfoil_module import output_reader


raw_data = output_reader("flexinol_monotonic_loading_martensite.csv", separator=",", 
                     rows_to_skip=4, header = ['Time', 'Extension',
                                               'Load', 
                                               "Strain",  "Stress"],)
#                                               'Load', 
#                                               "Strain", "Temperature", "Stress"],)
#Ignore initial data
for i in range(len(raw_data['Time'])):
    if raw_data['Stress'][i] > 100.  :
        break

young_data = {}
for key in raw_data:
    young_data[key] = raw_data[key][:i+1]

data = raw_data

#==============================================================================
# Fitting Young
#==============================================================================
(a, b)=polyfit(young_data["Strain"], young_data["Stress"], 1)
print "Young Modulus: ", a, b

fitted_strain = np.linspace(young_data["Strain"][0], 0.015)
fitted_stress = polyval([a,b], fitted_strain)
#==============================================================================
# Plotting
#==============================================================================
#plt.figure()
#plt.plot(data["Temperature"],data["Strain"])
#plt.xlabel("Temperature (C)")
#plt.ylabel("Strain (m/m)")
#plt.grid()

plt.figure()
plt.plot(data["Strain"], data["Stress"], lw = 2)
plt.plot(fitted_strain, fitted_stress, '--', color = "0.75", lw =2)
plt.xlabel("Strain (m/m)")
plt.ylabel("Stress (MPa)")
plt.grid()

plt.figure()
plt.plot(data["Temperature"], data["Stress"])
#plt.plot(smoothed_T, smoothed_sigma, 'g')
#plt.plot(T_interp, sigma_interp, 'r')
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")
plt.grid()

plt.figure()
plt.plot( np.array(data["Time"]) - data["Time"][0],data["Temperature"])
#plt.plot(xx, smoothed_T, 'g')
#plt.plot(xx - xx[0], T_interp, 'r')
plt.xlabel("Time (t)")
plt.ylabel("Temperature (C)")
plt.grid()