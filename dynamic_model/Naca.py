# -*- coding: utf-8 -*-

def Naca_00XX(c, t, x_list):
    """
    Generates a simetric NACA airfoil.
    Inputs:
    :param c: chord
    :param t: max thickness as a fraction of the chord.
    :param x_list: points between 0 and c (chord lenght) 
     
    Returns dictionary with keys 'u'pper and 'l'ower
    
    The Naca function can be found in: https://en.wikipedia.org/wiki/NACA_airfoil  

    Created on Wed Feb 03 12:50:52 2016
    
    @author: Endryws and Pedro Leal
    """
    y_upper = []
    y_lower = []
    for x in x_list:
        xc= x/c # Is just for increase speed and facilitate future changes.
        a1 = 5*t*c
        t1 = 0.2969*(math.sqrt(xc))
        t2 = -0.1260*xc
        t3 = -0.3516*(xc**2)
        t4 = 0.2843*(xc**3)
        t5 = -0.1015*(xc**4)
        y = (a1*(t1+t2+t3+t4+t5))
        y_upper.append(y)
        y_lower.append(y*(-1)) # is just for pick the y axis
                                       # negative numbers
    y = {'u': y_upper, 'l':y_lower}
    return y
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import math

    from airfoil_module import create_x
    
    c = 1. #c is the chord length
    t = 0.12 # is the maximum thickness as a fraction of the chord 
             #(t*100 = the last two nunbers of the NACA denomination)
    
    number_of_points = 200
    
    p = 0.5 #this variabel is here to force the points fall inside the model.
        
    x_list = create_x(c, n = number_of_points, distribution = 'polar')  # This function receives 
                                                   # inicial point,finalpoint,
                                                   #number of therms between 
                                                   # inicial and final point
    y = Naca_00XX(c,t,x_list)
    
    
    fig = plt.figure()

    fig.add_subplot(111)
    plt.scatter(x_list,y['u'])
    plt.scatter(x_list,y['l'])

    
    plt.grid()
    plt.xlabel('x values')
    plt.ylabel('y values')