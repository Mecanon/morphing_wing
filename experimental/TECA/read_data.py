# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:33:31 2016

@author: Pedro Leal
"""
import datetime
import time
import numpy as np
from scipy import polyfit, polyval
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from xfoil_module import output_reader

# Number of point to ignore at each voltage
N = 20
# Wire length
L = 0.19
# Wire radius
r = 0.000381/2.
# Wire cross section area
A = np.pi*r**2
# Start austenite temperature
A_s = 66.87
# Finish austenite temperature
A_f = 89.68

strain_correction = True

# Import data from IR camera
filename = "untrained_wire_temperature_run5.txt"
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

# Import data from Arduino
filename = "untrained_wire_voltage_current_run5.txt"
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
    Data_eletric["Voltage"][i] -= 1210.07069871 + 474.710578241 -949.421156482 - 504.768862708

# Since arduino time is included in IR time, to make both lists just
start_time = Data_eletric["Time"][0]
end_time = Data_eletric["Time"][-1]

start_index = Data_temperature['Time'].index(start_time)
end_index = len(Data_temperature['Time']) - Data_temperature['Time'][::-1].index(end_time) - 1

plt.figure()
plt.plot(Data_temperature['Time'], Data_temperature['Temperature'])
plt.plot(Data_eletric['Time'], Data_eletric['Current'])
plt.plot(Data_eletric['Time'], Data_eletric['Voltage'])
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

# Applying strain correction, L is a function of temperature.
if strain_correction:
    strain_data = output_reader("filtered_data_50MPa.txt", separator=",",
                               header = ['Temperature', "Strain", "Stress"],)
    f = interp1d(strain_data['Temperature'], strain_data["Strain"])
    
    for i in range(len(strain_data['Temperature'])):
        if filtered_Data['Temperature'][i] > strain_data['Temperature'][0]:
            break
    temperature_interp = filtered_Data['Temperature'][i:] 
    voltage_interp = filtered_Data['Voltage'][i:]
    current_interp = filtered_Data['Current'][i:]
    strain_interp = f(temperature_interp)
    delta_eps = max(strain_data['Strain']) - min(strain_data['Strain'])
current = filtered_Data["Current"]
voltage = filtered_Data["Voltage"]
temperature = filtered_Data["Temperature"]

# Find voltage where transformation occurs
f = interp1d(temperature, voltage)
V_As = f(A_s)
V_Af = f(A_f)

# Linear regression for Martensite
for i in range(len(temperature)):
    if temperature[i] > A_s:
        break
martensite_current = current[:i]
martensite_voltage = voltage[:i]
(a_M, b_M)=polyfit(martensite_voltage, martensite_current, 1)
print "Martensite resistivity  and resistnace, A0, V0: ", (1/a_M)*(A/L), (1/a_M), b_M, -b_M/a_M

# Linear regression for Austenite
for i in range(len(temperature)):
    if temperature[i] > A_f:
        break
austenite_current = current[i:]
austenite_voltage = voltage[i:]
(a_A, b_A)=polyfit(austenite_voltage, austenite_current, 1)
print "Austenite resistivity and resistance", (1/a_A)*(A/(L*(1-delta_eps))), (1/a_A), a_A

# Linear regression for all
(R, R_0)=polyfit(filtered_Data["Voltage"], filtered_Data["Current"], 1)
print "Linear regressed electrical resistivity", (1/R)*(A/L), (1/R)*(A/0.04)

#fitted_voltage= np.linspace(0, 500)
#fitted_current = polyval([a,b], fitted_voltage)

# Plot some results
plt.figure()
#plt.plot(fitted_voltage, fitted_current)

plt.axvline(V_As, color='r')
plt.axvline(V_Af, color='r')
plt.plot(filtered_Data["Voltage"], polyval((R, R_0), filtered_Data["Voltage"]),
         lw = 2, label = "Full linear regresion", c="0.5")
plt.plot(martensite_voltage, polyval((a_M, b_M), martensite_voltage), 'b',
         lw = 2, label = "Martensite linear regresion")
plt.plot(austenite_voltage, polyval((a_A, b_A), austenite_voltage),'g',
         lw = 2, label = "Austenite linear regresion")
plt.scatter(filtered_Data["Voltage"], filtered_Data["Current"],
            label = "Experimental data")
plt.xlabel("Voltage (mV)")
plt.ylabel("Current (mA)")
plt.grid()
plt.legend(loc='best')

plt.figure()
plt.plot(filtered_Data["Temperature"],
         np.array(filtered_Data["Voltage"])/np.array(filtered_Data["Current"]))
plt.axvline(A_s, color='r')
plt.axvline(A_f, color='r')
plt.xlabel(r"Temperature (${}^{\circ}$C)")
plt.ylabel(r"Electrical Resistance ($\Omega$)")
plt.grid()

plt.figure()
plt.plot(temperature_interp,
         np.array(voltage_interp)/np.array(current_interp)*(A/(L*(1-strain_interp)))*10**7)
plt.xlabel(r"Temperature (${}^{\circ}$C)")
plt.axvline(A_s, color='r')
plt.axvline(A_f, color='r')
plt.ylabel(r"Electrical Resistivity($\mu \Omega \times$ cm)")
plt.grid()

plt.figure()
plt.plot(filtered_Data["Current"], filtered_Data["Temperature"])
plt.axhline(A_s, color='r')
plt.axhline(A_f, color='r')
plt.xlabel("Current (mA)")
plt.ylabel(r"Temperature (${}^{\circ}$C)")
plt.grid()

