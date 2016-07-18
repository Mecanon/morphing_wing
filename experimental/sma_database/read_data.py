# -*- coding: utf-8 -*-
"""
Created on Wed Jun 08 15:06:17 2016

@author: Pedro Leal
"""
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.interpolate as interpolate
import numpy as np

from xfoil_module import output_reader

stress = "50MPa"
raw_data = output_reader("flexinol_isobaric_trained_" + stress + ".csv", separator=",", 
                     rows_to_skip=4, header = ['Time', 'Extension',
                                               'Load', "Temperature",
                                               "Strain", "Stress"],)
#first_data = output_reader("flexinol_untrained_cyclic_temperature_0p.csv", separator=",", 
#                     rows_to_skip=4, header = ['Time', 'Extension',
#                                               'Load',  "Strain",
#                                               "Stress", "Temperature"],)

##Ignore final data
#for i in range(len(first_data['Time'])):
#    if first_data['Time'][i] == 2908.75300: #120745.90000: 5958.89500
#        print 'ho'
#        break
#
#old_data = first_data
#data = {}
#for key in first_data:
#    data[key] = first_data[key][:i+1]
#
#first_data = data
#
## Ignore inital data
#for i in range(len(first_data['Time'])):
#    if first_data['Time'][i] == 1130.15300:
#        print 'hi'
#        break
#
#data = {}
#for key in first_data:
#    data[key] = first_data[key][i:]

#data = first_data
          
#==============================================================================
# Filtering data
#==============================================================================
#xx = np.linspace(data['Time'][0], data['Time'][-1],500)
#
#f = interpolate.interp1d(data['Time'],data['Strain'])
#eps_interp = f(xx)
#window = signal.gaussian(1, 1)
#smoothed_eps = signal.convolve(eps_interp, window/window.sum(), mode = 'same')

#f = interpolate.interp1d(data['Time'],data['Temperature'])
#T_interp = f(xx)
#window = signal.gaussian(1, 1)
#smoothed_T = signal.convolve(T_interp, window/window.sum(), mode = 'same')

#f = interpolate.interp1d(data['Time'],data['Stress'])
#sigma_interp = f(xx)
#window = signal.gaussian(1, 1)
#smoothed_sigma = signal.convolve(sigma_interp, window/window.sum(), mode = 'same')

#==============================================================================
# Formatting data
#==============================================================================
#for i in range(len(data["Stress"])):
#    if data["Stress"] > 100.:
#        break
#eps_0 = data["Strain"][i]
#
#eps_min = min(eps_interp)
#eps_interp = eps_interp - eps_min

#==============================================================================
# Plotting
#==============================================================================
plt.figure()
plt.plot(data["Temperature"],data["Strain"], label = 'First training')
#plt.plot(second_data["Temperature"],second_data["Strain"], label = 'Second training')
#plt.scatter(second_data["Temperature"][0],second_data["Strain"][0])
plt.xlabel("Temperature (C)")
plt.ylabel("Strain (m/m)")
plt.grid()
#plt.legend()

#plt.figure()
#plt.plot(data["Temperature"],data["Strain"])
#plt.xlabel("Temperature (C)")
#plt.ylabel("Strain (m/m)")
#plt.grid()
#
#plt.figure()
#plt.plot(first_data["Strain"], first_data["Stress"], label = 'First training')
#plt.plot(second_data["Strain"], second_data["Stress"], label = 'Second training')
#plt.xlabel("Strain (m/m)")
#plt.ylabel("Stress (MPa)")
#plt.grid()
#plt.legend()
#
plt.figure()
plt.plot(data["Temperature"], data["Stress"], label = 'First training')
#plt.plot(second_data["Temperature"], second_data["Stress"], label = 'Second training')
#plt.plot(smoothed_T, smoothed_sigma, 'g')
#plt.plot(T_interp, sigma_interp, 'r')
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")
plt.grid()
#plt.legend()
#
#plt.figure()
#plt.plot( np.array(first_data["Time"]) - data["Time"][0], first_data["Temperature"], label = 'First training')
#plt.plot( np.array(second_data["Time"]) - first_data["Time"][0], second_data["Temperature"], label = 'Second training')
##plt.plot(xx, smoothed_T, 'g')
##plt.plot(xx - xx[0], T_interp, 'r')
#plt.xlabel("Time (t)")
#plt.ylabel("Temperature (C)")
#plt.grid()
#plt.legend()
#
#plt.figure()
#plt.plot(first_data["Time"], first_data["Strain"], label = 'First training')
#plt.plot(second_data["Time"], second_data["Strain"], label = 'Second training')
##plt.plot(xx, smoothed_eps, 'g')
##plt.plot(xx - xx[0], eps_interp, 'r')
#plt.xlabel("Time (t)")
#plt.ylabel("Strain")
#plt.grid()
#plt.legend(loc=4)
#plt.figure()
##plt.plot(data["Temperature"],data["Strain"])
#plt.plot( smoothed_T, smoothed_eps, 'g')
#plt.plot(T_interp, eps_interp, 'r')
#plt.grid()
#plt.xlabel("Temperature")
#plt.ylabel("Strain")
#
##==============================================================================
## Stroing data
##==============================================================================
#data = np.array([T_interp, eps_interp, sigma_interp])
#try:
#    np.savetxt("filtered_data_"+ stress+".txt", data.T,fmt='%.18f')
#except:
#    print "No output file"