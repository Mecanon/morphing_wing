# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:53:47 2016

@author: Pedro Leal
"""
import numpy as np

from xfoil_module import output_reader

stress = "50MPa"
filename = "filtered_data_" + stress + ".txt"
max_172 = 0.10769
Data = output_reader(filename, separator='\t', output=None, rows_to_skip=1,
                     header=['Temperature', 'Strain', 'Stress'])

T = np.array(Data['Temperature'])
eps = np.array(Data['Strain'])*(0.045/max_172)
sigma = Data['Stress']

data = np.array([T, eps, sigma])

np.savetxt("filtered_data_"+ stress+".txt", data.T,fmt='%.18f')

print max(Data['Strain'])