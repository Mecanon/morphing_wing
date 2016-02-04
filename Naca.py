# -*- coding: utf-8 -*-
"""
Created on Wed Feb 03 12:50:52 2016

@author: Endryws
"""
c = 1 #c is the chord length
t = 0.012 # is the maximum thickness as a fraction of the chord 
         #(t*100 = the last two nunbers of the NACA denomination)

number_of_points = 200


p = 0.3 #I can´t remember how we are puting this p,
        # but I remeber Pedro say to put.

def Naca(c,t,xt):
    yt = []
    for x in xt:
        chimera = x/c # Is just for increase speed and facilitate future changes.
        a1 = 5*t*c
        t1 = 0.2969*(math.sqrt(chimera))
        t2 = -0.1260*chimera
        t3 = -0.3516*(chimera**2)
        t4 = 0.2843*(chimera**3)
        t5 = -0.1015*(chimera**4)
        y = (a1*(t1+t2+t3+t4+t5))*p
        yt.append(y)
    return yt
 

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    
    xt = np.linspace(.001,c,number_of_points) # This function receives inicial
                                          # point,finalpoin,number of therms 
                                         #between inicial and final point
    
    yt = Naca(c,t,xt)
    
    fig = plt.figure()

    fig.add_subplot(111)
    plt.scatter(xt,yt)

    yt_neg = [ -y for y in yt] # is just for pick the negative numbers
    plt.scatter(xt,yt_neg)
    plt.grid()
    plt.xlabel('xt')
    plt.ylabel('yt')