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

raw_data = output_reader("flexinol_isobaric_training_172MPa_2.csv", separator=",", 
                     rows_to_skip=4, header = ['Time', 'Extension',
                                               'Load', "Temperature",
                                               "Strain", "Stress"],)

#Ignore initial data
for i in range(len(raw_data['Time'])):
    if raw_data['Time'][i] == 111471.60000:
        break

data = {}
for key in raw_data:
    data[key] = raw_data[key][i:]

data = raw_data
#Ignore final data
for i in range(len(data['Time'])):
    if data['Time'][i] == 120745.90000 : #120745.90000: 5958.89500
        break

old_data = data
data = {}
for key in old_data:
    data[key] = old_data[key][:i+1]

#==============================================================================
# Filtering data
#==============================================================================
xx = np.linspace(data['Time'][0], data['Time'][-1],500)
print data['Time'][20000],data['Time'][-20000],len(data['Time'])
f = interpolate.interp1d(data['Time'],data['Strain'])
eps_interp = f(xx)
window = signal.gaussian(10, 1000)
smoothed_eps = signal.convolve(eps_interp, window/window.sum(), mode = 'same')

f = interpolate.interp1d(data['Time'],data['Temperature'])
T_interp = f(xx)
window = signal.gaussian(10, 1000)
smoothed_T = signal.convolve(T_interp, window/window.sum(), mode = 'same')

f = interpolate.interp1d(data['Time'],data['Stress'])
sigma_interp = f(xx)
window = signal.gaussian(10, 1000)
smoothed_sigma = signal.convolve(sigma_interp, window/window.sum(), mode = 'same')

#==============================================================================
# Formatting data
#==============================================================================
for i in range(len(data["Stress"])):
    if data["Stress"] > 100.:
        break
eps_0 = data["Strain"][i]

eps_min = min(eps_interp)
eps_interp = eps_interp - eps_min + eps_0

#==============================================================================
# Plotting
#==============================================================================
plt.figure()
plt.plot(data["Temperature"],data["Strain"])
plt.xlabel("Temperature (C)")
plt.ylabel("Strain (m/m)")
plt.grid()

plt.figure()
plt.plot(data["Strain"], data["Stress"])
plt.xlabel("Strain (m/m)")
plt.ylabel("Stress (MPa)")

plt.figure()
plt.plot(data["Temperature"], data["Stress"])
plt.plot(smoothed_T, smoothed_sigma, 'g')
#plt.plot(T_interp, sigma_interp, 'r')
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")

plt.figure()
plt.plot( data["Time"],data["Temperature"])
plt.plot(xx, smoothed_T, 'g')
plt.plot(xx, T_interp, 'r')
plt.xlabel("Time")
plt.ylabel("Temperature (C)")

plt.figure()
plt.plot(data["Time"], data["Strain"])
plt.plot(xx, smoothed_eps, 'g')
plt.plot(xx, eps_interp, 'r')
plt.xlabel("Time")
plt.ylabel("Strain")

plt.figure()
#plt.plot(data["Temperature"],data["Strain"])
#plt.plot( smoothed_T, smoothed_eps, 'g')
plt.plot(T_interp, eps_interp, 'r')
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Strain")

#==============================================================================
# Stroing data
#==============================================================================
data = np.array([T_interp, eps_interp, sigma_interp])
np.savetxt("filtered_data_500MPa.txt", data.T,fmt='%.18f')