# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:31:45 2016

@author: Pedro Leal
"""

import matlab.engine

#Start Matlab engine
eng = matlab.engine.start_matlab()
#Go to directory where matlab file is
eng.cd('SMA_temperature_strain_driven')
#Run SMA model
data = eng.OneD_SMA_Model('True')

