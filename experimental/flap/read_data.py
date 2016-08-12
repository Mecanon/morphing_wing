# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:33:31 2016

@author: Pedro Leal
"""
import math
import datetime
import time
import pickle
import numpy as np
from scipy import polyfit, polyval
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from xfoil_module import output_reader


#Data = {'theta': theta, 'eps_s': eps_s, 'eps_l': eps_l, 
#                'sigma': sigma, 'xi': MVF, 'T': T, 'eps_t': eps_t,
#                'F_l': F_l, 'k': k, 'L_s':L_s}
Data = pickle.load( open( "data.p", "rb" ) )
T_num = np.array(Data['T']) - 273.15 
theta_num = np.rad2deg(np.array(Data['theta']))

# Wire properties
rho = 3.55041728247e-06
A = math.pi*(0.000381/2.)**2
L = math.sqrt((0.168 + 0.018)**2 +  (0.06 + 0.01)**2)

# Number of point to ignore at each voltage
N = 2
# Import data from IR camera
filename = "voltage_angle_run4.txt"
Data_arduino = output_reader(filename, separator='\t', output=None, rows_to_skip=1,
                     header=['Date', 'Time', 'Voltage', 'delta t', 'Z'], column_types = [str, str, float,
                     float, float])

# Convert time to total seconds
for i in range(len(Data_arduino['Time'])):
    time_i = time.strptime(Data_arduino['Time'][i].split(',')[0],
                           '%H:%M:%S')
    Data_arduino['Time'][i] = datetime.timedelta(hours=time_i.tm_hour, 
                                  minutes=time_i.tm_min, 
                                  seconds=time_i.tm_sec).total_seconds()
    
# Import data from IR camera
filename = "temperature_run4.txt"
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

# Since arduino time is included in IR time, to make both lists just
start_time = Data_arduino["Time"][0]
end_time = Data_arduino["Time"][-1]

start_index = Data_temperature['Time'].index(start_time)
end_index = len(Data_temperature['Time']) - Data_temperature['Time'][::-1].index(end_time) - 1

# Eliminate temperature data for where there is no arduino data
new_Data = {}
for key in Data_temperature:
    new_Data[key] = Data_temperature[key][start_index:end_index+1]
Data_temperature = new_Data

#==============================================================================
# Plot raw data
#==============================================================================

fig, ax1 = plt.subplots()
ax1.plot(np.array(Data_temperature['Time']) - min(Data_temperature['Time']), Data_temperature['Temperature'], 'b-')
ax1.set_xlabel('Time (s)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel("Temperature ($^{\circ}C$)", color='b')
ax1.set_xlim([min(Data_temperature['Time']), max(Data_temperature['Time'])])
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(np.array(Data_arduino['Time']) - min(Data_arduino['Time']), Data_arduino['Voltage'], 'r')
ax2.set_ylabel("Voltage (mV)", color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()


fig, ax1 = plt.subplots()
ax1.plot(np.array(Data_temperature['Time']) - min(Data_temperature['Time']), Data_temperature['Temperature'], 'b-')
ax1.set_xlabel('Time (s)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel("Temperature ($^{\circ}C$)", color='b')
ax1.set_xlim([min(Data_temperature['Time']), max(Data_temperature['Time'])])
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(np.array(Data_arduino['Time']) - min(Data_arduino['Time']), Data_arduino['Z'], 'r')
ax2.set_ylabel("Deflection ($^{\circ}$)", color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()



#==============================================================================
# # Average values for temperature so that both datas have same length
#==============================================================================
averaged_temperature_data= {'Time': [], 'Temperature':[]}

temperature_readings = []
current_time = 0
for i in range(len(Data_temperature['Temperature'])):
    if Data_temperature['Time'][i] != current_time:
        if current_time != 0:
            averaged_temperature_data['Temperature'].append(sum(temperature_readings)/len(temperature_readings))
        
        current_time = Data_temperature['Time'][i]
        averaged_temperature_data['Time'].append(current_time)
        temperature_readings = [Data_temperature['Temperature'][i]]
    else:
        temperature_readings.append(Data_temperature['Temperature'][i])          
# If last set of data not used, store it
if len(averaged_temperature_data['Time']) != len(averaged_temperature_data['Temperature']):
    averaged_temperature_data['Temperature'].append(sum(temperature_readings)/len(temperature_readings))   

#==============================================================================
# # Average values for voltage and deflection for same time
#==============================================================================
averaged_arduino_data= { 'Time':[], 'Voltage': [], 'Deflection': []}

voltage_readings = []
deflection_readings = []
current_time = 0
for i in range(len(Data_arduino['Z'])):
    if Data_arduino['Time'][i] != current_time:
        if current_time != 0:
            averaged_arduino_data['Deflection'].append(sum(deflection_readings)/len(deflection_readings))
            averaged_arduino_data['Voltage'].append(sum(voltage_readings)/len(voltage_readings))
        current_time = Data_arduino['Time'][i]
        averaged_arduino_data['Time'].append(current_time)
        deflection_readings = [Data_arduino['Z'][i]]
        voltage_readings = [Data_arduino['Voltage'][i]]
    else:
        deflection_readings.append(Data_arduino['Z'][i])
        voltage_readings.append(Data_arduino['Voltage'][i])
# If last set of data not used, store it
if len(averaged_arduino_data['Deflection']) != len(averaged_arduino_data['Time']):
    averaged_arduino_data['Deflection'].append(sum(deflection_readings)/len(deflection_readings))
    averaged_arduino_data['Voltage'].append(sum(voltage_readings)/len(voltage_readings))   
         
#==============================================================================
# Same length dataset
#==============================================================================
all_Data = {'Voltage': [], 'Deflection': [], 'Temperature': [], 'Time': []}

all_Data['Time'] = averaged_arduino_data['Time']
all_Data['Voltage'] = averaged_arduino_data['Voltage']
all_Data['Deflection'] = averaged_arduino_data['Deflection']
f = interp1d(averaged_temperature_data['Time'], averaged_temperature_data['Temperature'])
all_Data['Temperature'] = f(np.array(averaged_arduino_data['Time']))

#print len(all_Data['Voltage']), len(all_Data['Deflection']), len(all_Data['Temperature'])   

plt.figure()
plt.scatter(min(all_Data['Temperature']) + np.array(all_Data['Temperature'] - min(all_Data['Temperature']))*1.5,all_Data['Deflection'], label = 'Post-processed experiment')
plt.plot(T_num, theta_num, label = 'Mathematical model')
plt.xlabel("Temperature ($^{\circ}C$)")
plt.ylabel("Flap deflection ($^{\circ}$)")
plt.grid()
plt.legend(loc = 'best')
#==============================================================================
# Ignore N first measurement for each voltage
#==============================================================================
current_voltage = 0
counter = 1
filtered_Data = {}

time = []
deflection = []
voltage = []
temperature = []

time_readings = []
deflection_readings = []
temperature_readings = []
            
for i in range(len(all_Data['Deflection'])):
    print current_voltage, all_Data['Voltage'][i] 
    if all_Data['Voltage'][i] != current_voltage:
        if current_voltage != 0:
            print sum(time_readings), len(time_readings), time_readings
            time.append(sum(time_readings)/len(time_readings))
            deflection.append(sum(deflection_readings)/len(deflection_readings))
            temperature.append(sum(temperature_readings)/len(temperature_readings))
            voltage.append(current_voltage)   
        counter = 1
        current_voltage = all_Data['Voltage'][i]
    else:
        if counter >N:
            time_readings.append(all_Data['Time'][i]) 
            deflection_readings.append(all_Data['Deflection'][i])
            temperature_readings.append(all_Data['Temperature'][i]) 
        counter += 1
# If last set of data not used, store it
if current_voltage != voltage[-1]:
    voltage.append(current_voltage)
    time.append(sum(time_readings)/len(time_readings))
    deflection.append(sum(deflection_readings)/len(deflection_readings))
    temperature.append(sum(temperature_readings)/len(temperature_readings))   
    
print len(voltage), len(deflection), len(temperature)  

#==============================================================================
# Plotting
#==============================================================================
time = np.array(time) - time[0]

fig, ax1 = plt.subplots()
ax1.plot(time, deflection, 'b-')
ax1.set_xlabel('Time (s)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel("Flap deflection (${}^{\circ}$)", color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(time, voltage, 'r')
ax2.set_ylabel("Voltage (mV)", color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.grid()
plt.show()

plt.figure()
plt.scatter(voltage, deflection)
plt.xlabel("Voltage (mV)")
plt.ylabel("Flap deflection (${}^{\circ}$)")
plt.grid()

plt.figure()
plt.scatter(voltage, temperature)
plt.xlabel("Voltage (mV)")
plt.ylabel("Temperature (${}^{\circ}C$)")
plt.grid()

plt.figure()
plt.scatter(temperature, deflection)
plt.xlabel("Temperature (${}^{\circ}$C)")
plt.ylabel("Flap deflection (${}^{\circ}$)")
plt.grid()

current = A*np.array(voltage)/(rho*L)

plt.figure()
plt.scatter(current, temperature)
plt.xlabel("Current (mA)")
plt.ylabel("Temperature (${}^{\circ}$C)")
plt.grid()

#fig, ax1 = plt.subplots()
#ax1.plot(time, deflection, 'b-')
#ax1.set_xlabel('Time (s)')
## Make the y-axis label and tick labels match the line color.
#ax1.set_ylabel("Flap deflection (${}^{\circ}$)", color='b')
#for tl in ax1.get_yticklabels():
#    tl.set_color('b')
#
#
#ax2 = ax1.twinx()
#ax2.plot(time, voltage, 'r')
#ax2.set_ylabel("Voltage (mV)", color='r')
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')
#plt.grid()
#plt.show()