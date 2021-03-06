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

#stress = "50MPa"
#raw_data = output_reader("flexinol_isobaric_trained_" + stress + ".csv", separator=",", 
#                     rows_to_skip=4, header = ['Time', 'Extension',
#                                               'Load', "Temperature",
#                                               "Strain", "Stress"],)
first_data = output_reader("flexinol_training_first.csv", separator=",", 
                     rows_to_skip=4, header = ['Time', 'Extension',
                                               'Load',  "Temperature",
                                               "Strain", "Stress"],)

second_data = output_reader("flexinol_training_second.csv", separator=",", 
                     rows_to_skip=4, header = ['Time', 'Extension',
                                               'Load',  "Temperature",
                                               "Strain", "Stress"],)
#==============================================================================
# CORRECTION FACTOR BECAUSE OF INCORRENT SAMPLE LENGTH AT INSTRON
#==============================================================================
first_data["Strain"] = list(np.array(first_data["Strain"])*( 0.045/0.10769))

second_data["Strain"] = list(np.array(second_data["Strain"])*( 0.045/0.10769))

#==============================================================================
# Plot first and last with same scale
#==============================================================================
plt.figure()
plt.plot(first_data["Temperature"][11600:29078], first_data["Strain"][11600:29078])
plt.xlabel("Temperature ($^{\circ}C$)", fontsize = 16)
plt.ylabel("$\epsilon$", fontsize = 26)

plt.grid()
plt.ylim(-0.045,0.02)
plt.xlim(30,140)
plt.figure()
plt.plot(second_data["Temperature"][589527:615819], second_data["Strain"][589527:615819])
plt.xlabel("Temperature ($^{\circ}C$)", fontsize = 16)
plt.ylabel("$\epsilon$", fontsize = 26)
plt.grid()
plt.ylim(-0.045,0.02)
plt.xlim(30,140)

#==============================================================================
# 
#==============================================================================
#Ignore initial data
for i in range(len(first_data['Time'])):
    if first_data['Stress'][i] > 100.  :
        break
    
#Ignore final data
for i in range(len(first_data['Time'])):
    if first_data['Time'][i] == 54470.46: #120745.90000: 5958.89500
        print 'ho'
        break

old_data = first_data
data = {}
for key in first_data:
    data[key] = first_data[key][:i+1]

first_data = data

# Ignore inital data
for i in range(len(second_data['Time'])):
    if second_data['Time'][i] == 1601.507  :
        print 'hi'
        break

data = {}
for key in second_data:
    data[key] = second_data[key][i:]

second_data= data

# Translating second data to match first data
second_data['Strain'] = list(np.array(second_data['Strain']) + \
                        (first_data['Strain'][-1] - second_data['Strain'][0]))
second_data['Time'] = list(np.array(second_data['Time']) + \
                        (first_data['Time'][-1] - second_data['Time'][0]))
# Joining both datas
data = {}
for key in second_data:
    data[key] = first_data[key] + second_data[key]              
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
plt.plot(first_data["Temperature"],first_data["Strain"], label = 'First part')
plt.plot(second_data["Temperature"],second_data["Strain"], label = 'Second part')
plt.scatter(first_data["Temperature"][-1],first_data["Strain"][-1])
#plt.scatter(second_data["Temperature"][0],second_data["Strain"][0])
plt.xlabel("Temperature (C)")
plt.ylabel("Strain (m/m)")
plt.grid()
plt.legend()

plt.figure()
plt.plot(data["Temperature"],data["Strain"])
plt.xlabel("Temperature (C)")
plt.ylabel("Strain (m/m)")
plt.grid()

plt.figure()
plt.plot(first_data["Strain"], first_data["Stress"], label = 'First part')
plt.plot(second_data["Strain"], second_data["Stress"], label = 'Second part')
plt.xlabel("Strain (m/m)")
plt.ylabel("Stress (MPa)")
plt.grid()
plt.legend()

plt.figure()
plt.plot(first_data["Temperature"], first_data["Stress"], label = 'First part')
plt.plot(second_data["Temperature"], second_data["Stress"], label = 'Second part')
#plt.plot(smoothed_T, smoothed_sigma, 'g')
#plt.plot(T_interp, sigma_interp, 'r')
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")
plt.grid()
plt.legend()

plt.figure()
plt.plot( np.array(first_data["Time"]) - data["Time"][0], first_data["Temperature"], label = 'First part')
plt.plot( np.array(second_data["Time"]) - first_data["Time"][0], second_data["Temperature"], label = 'Second part')
#plt.plot(xx, smoothed_T, 'g')
#plt.plot(xx - xx[0], T_interp, 'r')
plt.xlabel("Time (t)")
plt.ylabel("Temperature (C)")
plt.grid()
plt.legend()

plt.figure()
plt.plot(first_data["Time"], first_data["Strain"], label = 'First part')
plt.plot(second_data["Time"], second_data["Strain"], label = 'Second part')
#plt.plot(xx, smoothed_eps, 'g')
#plt.plot(xx - xx[0], eps_interp, 'r')
plt.xlabel("Time (t)")
plt.ylabel("Strain")
plt.grid()
plt.legend(loc=4)
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