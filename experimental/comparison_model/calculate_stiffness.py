# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:42:47 2016

@author: Pedro Leal
"""
#Total number of springs
Nt = 67. 

safety_factor = 5.
A = 2211e6*10**(0.435) #Mpa.mm**0.125

#wire diameter
d = 0.00075
# Spring diameter
D = 0.00725
# Spring index
C = D/d
# zero length
L_o = 0.05880

if d < 0.000838:
    G = 82.7e9
    E = 203.4e9
elif d < 0.0016:
    G = 81.7e9
    E = 200.0e9
elif d < 0.00318:
    G = 81.0e9
    E = 196.6e9
else:
    G = 80.0e9
    E = 193.0e9

Na = Nt + G/E #Equivalent active number of springs
                
k = G*d**4/(8*Na*D**3)

print "spring coefficient:", k