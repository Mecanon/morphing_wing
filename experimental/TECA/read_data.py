# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:33:31 2016

@author: Pedro Leal
"""
import datetime
import time
import numpy as np
from scipy import polyfit, polyval
import matplotlib.pyplot as plt

from xfoil_module import output_reader

# Number of point to ignore at each voltage
N = 20
# Import data from IR camera
filename = "untrained_wire_temperature_run2.txt"
Data_temperature = output_reader(filename, separator='\t', output=None, rows_to_skip=13,
                     header=['Date', 'Time', 'Miliseconds', 'Relative time', 
                     'Temperature'], column_types = [str, str, int, float,
                     float, float])

# Convert time to total seconds
for i in range(len(Data_temperature['Time'])):
    time_i = time.strptime(Data_temperature['Time'][i].split(',')[0],
                           '%H:%M:%S')
    Data_temperature['Time'][i] = datetime.timedelta(hours=time_i.tm_hour, 
                                  minutes=time_i.tm_min, 
                                  seconds=time_i.tm_sec).total_seconds()

# Import data from Ir camera
filename = "untrained_wire_voltage_current_run2.txt"
Data_eletric = output_reader(filename, separator='\t', output=None, rows_to_skip=1,
                     header=['Date', 'Time',  'Voltage', 
                     'Current'], column_types = [str, str,  float,
                     float])

# Convert time to total seconds
for i in range(len(Data_eletric ['Time'])):
    time_i = time.strptime(Data_eletric ['Time'][i].split(',')[0],
                           '%H:%M:%S')
    Data_eletric['Time'][i] = datetime.timedelta(hours=time_i.tm_hour, 
                                  minutes=time_i.tm_min, 
                                  seconds=time_i.tm_sec).total_seconds()
    
# Fix problem that voltage "drops in the middle"

for i in range(1800, len(Data_eletric["Voltage"])):
    Data_eletric["Voltage"][i] += 3000.

# There is a voltage offset
for i in range(len(Data_eletric["Voltage"])):
    Data_eletric["Voltage"][i] -= 1210.07069871

# Since arduino time is included in IR time, to make both lists just
start_time = Data_eletric["Time"][0]
end_time = Data_eletric["Time"][-1]

start_index = Data_temperature['Time'].index(start_time)
end_index = len(Data_temperature['Time']) - Data_temperature['Time'][::-1].index(end_time) - 1

plt.figure()
plt.plot(Data_temperature['Time'], Data_temperature['Temperature'])
plt.plot(Data_eletric['Time'], Data_eletric['Current'])

new_Data = {}
for key in Data_temperature:
    new_Data[key] = Data_temperature[key][start_index:end_index+1]
Data_temperature = new_Data

# Average values for temperature so that both datas have same length
all_Data= {'Voltage': Data_eletric['Voltage'], 'Current': Data_eletric['Current'], 'Temperature':[]}
temperature_readings = []
current_time = 0
for i in range(len(Data_temperature['Temperature'])):
    if Data_temperature['Time'][i] != current_time:
        if current_time != 0:
            all_Data['Temperature'].append(sum(temperature_readings)/len(temperature_readings))
        current_time = Data_temperature['Time'][i]
        temperature_readings = [Data_temperature['Temperature'][i]]
    else:
        temperature_readings.append(Data_temperature['Temperature'][i])          
print len(all_Data['Voltage']), len(all_Data['Current']), len(all_Data['Temperature'])

# Ignore 20 first measurement for each voltage
current_voltage = 0
counter = 1
filtered_Data = {}
current_readings = []
temperature_readings = []
filtered_Data = {'Voltage':[], 'Current':[], 'Temperature':[]}
for i in range(len(all_Data['Current'])):
    if Data_eletric['Voltage'][i] != current_voltage:
        if current_voltage != 0:
            filtered_Data['Voltage'].append(current_voltage)
            filtered_Data['Current'].append(sum(current_readings)/len(current_readings))
            filtered_Data['Temperature'].append(sum(temperature_readings)/len(temperature_readings))
        counter = 1
        current_voltage = all_Data['Voltage'][i]
        current_readings = []
        temperature_readings = []
    else:
        if counter >N:
            current_readings.append(all_Data['Current'][i])
            temperature_readings.append(all_Data['Temperature'][i])   
        counter += 1

# Fit a line for initial data
(a, b)=polyfit(filtered_Data["Voltage"][3:7], filtered_Data["Current"][3:7], 1)
print "Initial resistivity, A0, V0: ", a, b, -b/a

fitted_voltage= np.linspace(0, 500)
fitted_current = polyval([a,b], fitted_voltage)

# Plot some results
plt.figure()
plt.plot(fitted_voltage, fitted_current)
plt.scatter(filtered_Data["Voltage"], filtered_Data["Current"])
plt.grid()

plt.figure()
plt.plot(filtered_Data["Temperature"],np.array(filtered_Data["Voltage"])/np.array(filtered_Data["Current"])/0.05)
plt.grid()

plt.figure()
plt.plot(filtered_Data["Current"], filtered_Data["Temperature"])
plt.grid()