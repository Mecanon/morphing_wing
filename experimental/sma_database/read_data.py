# -*- coding: utf-8 -*-
"""
Created on Wed Jun 08 15:06:17 2016

@author: Pedro Leal
"""
import matplotlib.pyplot as plt

from xfoil_module import output_reader

raw_data = output_reader("flexinol_untrained_isobaric_172MPa.csv", separator=",", 
                     rows_to_skip=3, header = ['Time', 'Extension',
                                               'Load', "Temperature",
                                               "Strain", "Stress"],)

#Use data only after T>80
for i in range(len(raw_data['Temperature'])):
    if raw_data['Temperature'][i]> 40.:
        break

data = {}
for key in raw_data:
    data[key] = raw_data[key][i:]

plt.figure()
plt.plot(data["Strain"], data["Temperature"])
plt.xlabel("Strain")
plt.ylabel("Temperature (C)")
plt.grid()

plt.figure()
plt.plot(data["Strain"], data["Stress"])
plt.xlabel("Strain")
plt.ylabel("Stress (MPa)")

plt.figure()
plt.plot(data["Temperature"], data["Stress"])
plt.xlabel("Temperature (C)")
plt.ylabel("Stress (MPa)")