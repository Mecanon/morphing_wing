# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 15:46:41 2016

@author: Pedro Leal
"""

import matlab.engine
#Start Matlab engine
eng = matlab.engine.start_matlab()

a_list = [1,2,3,4]
back_list = eng.test(a_list, nargout = 0)