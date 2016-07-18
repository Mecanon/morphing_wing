# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:39:31 2016

@author: Eduardo Tancredo & Pedro Leal
"""
import math
import numpy as np
#import matplotlib.pyplot as plt

n_1 = 0.1919
n_2 = 0.1823
n_3 = 0.1623
n_4 = 0.2188

Ms_list = np.array([ 68.74, 75.71, 82.33, 84.77, 88.27]) 
Mf_list = np.array([ 57.74, 65.39, 71.29, 74.07, 77.88]) 
As_list = np.array([ 78.47, 83.82, 88.81, 91.38, 94.78]) 
Af_list = np.array([ 88.75, 95.02, 102.15, 105.12, 108.85]) 

calibrated_Ms = np.array([ 88.44,  112.6]) 
calibrated_Mf = np.array([ 24.33,  48.49]) 
calibrated_As = np.array([ 45.21,  66.6]) 
calibrated_Af = np.array([ 113.7,  135.09])
 
def smoothed_temperatures(Ms_list , Mf_list , As_list , Af_list, 
                          print_output = True):
    """Convert tangent line results into smoothed line results"""
#   smoothed lists
    MsSmooth_list = [ ]
    MfSmooth_list = [ ]
    AsSmooth_list = [ ]
    AfSmooth_list = [ ]
    
#   creating the smoothed lists
    for i in range(len(Ms_list)):
        
        MsSmooth_list.append(  ( Ms_list[i]/2  *  (    ( 1+   ( 2**(-n_1))*(n_1+1)   +    (2**(-n_2))*(n_2-1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2))))     +    Mf_list[i]/2  *  (    ( -1+   ( 2**(-n_1))*(n_1-1)   +    (2**(-n_2))*(n_2+1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2))))
        MfSmooth_list.append(  ( Ms_list[i]/2  *  (    ( -1+   ( 2**(-n_1))*(n_1+1)   +    (2**(-n_2))*(n_2-1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2))))     +    Mf_list[i]/2  *  (    ( 1+   ( 2**(-n_1))*(n_1-1)   +    (2**(-n_2))*(n_2+1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2))))
        AsSmooth_list.append(  ( As_list[i]/2  *  (    ( 1+   ( 2**(-n_3))*(n_3-1)   +    (2**(-n_4))*(n_4+1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4))))     +    Af_list[i]/2  *  (    ( -1+   ( 2**(-n_3))*(n_3+1)   +    (2**(-n_4))*(n_4-1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4))))
        AfSmooth_list.append(  ( As_list[i]/2  *  (    ( -1+   ( 2**(-n_3))*(n_3-1)   +    (2**(-n_4))*(n_4+1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4))))     +    Af_list[i]/2  *  (    ( 1+   ( 2**(-n_3))*(n_3+1)   +    (2**(-n_4))*(n_4-1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4))))

    if print_output:
        print "smoothed_Ms: ", Ms_list
        print "smoothed_Mf: ", Mf_list
        print "smoothed_As: ", As_list
        print "smoothed_Af: ", Af_list
        
    
    return MsSmooth_list, MfSmooth_list, AsSmooth_list, AfSmooth_list

def tangent_temperatures(MsSmooth_list, MfSmooth_list, AsSmooth_list, 
                         AfSmooth_list, print_output = True):
    """Convert smoothed line results into tangent line results"""
#   lists created on the reverse problem
    mst=[]
    mft=[]
    ast=[]
    aft=[]
#   creating the reverse lists
    for i in range(len(MsSmooth_list)):
        k1 = (( 1+   ( 2**(-n_1))*(n_1+1)   +    (2**(-n_2))*(n_2-1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2)))/2
        k2 = (( -1+   ( 2**(-n_1))*(n_1-1)   +    (2**(-n_2))*(n_2+1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2)))/2
        k3 = (( -1+   ( 2**(-n_1))*(n_1+1)   +    (2**(-n_2))*(n_2-1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2)))/2 
        k4 = (( 1+   ( 2**(-n_1))*(n_1-1)   +    (2**(-n_2))*(n_2+1))  /  ( n_1 * (2**-n_1)  +  n_2*(2**-n_2)))/2
        k5 = (( 1+   ( 2**(-n_3))*(n_3-1)   +    (2**(-n_4))*(n_4+1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4)))/2
        k6 = (( -1+   ( 2**(-n_3))*(n_3+1)   +    (2**(-n_4))*(n_4-1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4)))/2
        k7 = (( -1+   ( 2**(-n_3))*(n_3-1)   +    (2**(-n_4))*(n_4+1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4)))/2
        k8 = (( 1+   ( 2**(-n_3))*(n_3+1)   +    (2**(-n_4))*(n_4-1))  /  ( n_3 * (2**-n_3)  +  n_4*(2**-n_4)))/2
    
    
        q = np.array([[k1,k3], [k2,k4]])
        l = np.array([[k5,k6], [k7,k8]])
        m = np.array([MsSmooth_list[i],MfSmooth_list[i]])
        a = np.array([AsSmooth_list[i],AfSmooth_list[i]])
    
        x = np.linalg.solve(q, m)
        mst.append(x[0])
        mft.append(x[1])
        y = np.linalg.solve(l, a)
        ast.append(y[0])
        aft.append(y[1])

    if print_output:
        print "tangent_Ms: ", mst
        print "tangent_Mf: ", mft
        print "tangent_As: ", ast
        print "tangent_Af: ", aft
        
    return mst, mft, ast, aft
   
#smoothed_temperatures(Ms_list , Mf_list , As_list , Af_list )
tangent_temperatures(calibrated_Ms , calibrated_Mf , calibrated_As , calibrated_Af )