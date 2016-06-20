# -*- coding: utf-8 -*-
"""
Created on Wed Jun 08 15:06:17 2016

@author: Pedro Leal
"""
import matplotlib.pyplot as plt

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
    
#Ignore final data
for i in range(len(data['Time'])):
    if data['Time'][i] ==  120745.90000:
        break

old_data = data
data = {}
for key in old_data:
    data[key] = old_data[key][:i+1]

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
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")