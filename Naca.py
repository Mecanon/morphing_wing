# -*- coding: utf-8 -*-
"""
Created on Wed Feb 03 12:50:52 2016

@author: Endryws
"""
def Naca(c, t, x_list):
    yt = []
    yt_inside = []
    for x in x_list:
        xc = x/c # Is just for increase speed and facilitate future changes.
        a1 = 5*t*c
        t1 = 0.2969*(math.sqrt(xc))
        t2 = -0.1260*xc
        t3 = -0.3516*(xc**2)
        t4 = 0.2843*(xc**3)
        t5 = -0.1015*(xc**4)
        y = (a1*(t1+t2+t3+t4+t5))
        y_inside = y*p
        yt.append(y)
        yt_inside.append(y_inside)
    return yt,yt_inside
 

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import math

    c = 1 #c is the chord length
    t = 0.12 # is the maximum thickness as a fraction of the chord 
             #(t*100 = the last two nunbers of the NACA denomination)
    
    number_of_points = 200
    
    p = 0.5 #this variabel is here to force the points fall inside the model.
        
    xt = np.linspace(.001,c,number_of_points) # This function receives inicial
                                          # point,finalpoin,number of therms 
                                         #between inicial and final point
    
    yt,yt_inside = Naca(c,t,xt)
    
    fig = plt.figure()

    fig.add_subplot(111)
    plt.scatter(xt,yt)
    plt.scatter(xt,yt_inside,c='r' )

    yt_neg = [ -y for y in yt] # is just for pick the negative numbers
    yt_inside_neg = [ -y for y in yt_inside]
    
    plt.scatter(xt,yt_neg)
    plt.scatter(xt,yt_inside_neg,c='r')
    
    plt.grid()
    plt.xlabel('xt')
    plt.ylabel('yt')