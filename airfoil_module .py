# -*- coding: utf-8 -*-
import math 

def CST(x,deltasz,A0s,A1s,c):
    """
    Based on the paper ""Fundamental" Parametric Geometry Representations for
    Aircraft Component Shapes" from Brenda M. Kulfan and John E. Bussoletti. 
    The code uses a 1st order Bernstein Polynomial for the “Class Function” / 
    “Shape Function” airfoil representation.
    
    The inputs are:

        - x:list of points along the chord, from TE and the LE, or vice-versa.
          The code works both ways.
    
    For deltaz, A0s and A1s: first element in the list is related to the upper
    surface and the second to the lower surface. There are two because the CST
    method treats the airfoil surfaces as two different surfaces (upper and 
    lower)
    
        - deltasz: list of thicknesses on the TE.
        
        - A0s: list of A0 coefficients, which are design parameters.
        
        - A1s: list of A1 coefficients, which are design parameters.
        
        - c: chord
    
    The outputs are:
        - y_u: list of y position of the upper surface according to the x 
        points.
        
        - y_l: same as y_u, but for the lower surface.
        
    Created on Sun Jan 19 16:36:55 2014

    @author: Pedro Leal
    """
    def K(r,n):
        K= math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K        
    def S(r,n,psi):
        S=K(r,n)*(psi**r)*(1-psi)**(n-r)
        return S
        
    def C(N1,N2,psi):
        C=((psi)**N1)*((1-psi)**N2)
        return C
     
    psi=x/c;
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Treating inputs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Giving names to everything
    A0_u=A0s[0]
    A0_l=A0s[1]
    A1_u=A1s[0]
    A1_l=A1s[1]
    deltaz_u=deltasz[0]
    deltaz_l=deltasz[1]
    
    # Adimensionalizing
    psi=x/c;    
    deltaeta_u=deltaz_u/c
    deltaeta_l=deltaz_l/c # Will end up contributing negatively
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Class Function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # The Coefficients for an airfoil with a rounded leading edge and a sharp
    # trailing edge are N1=0.5 and N2=1.0.
    N1=0.5;
    N2=1.0;
    n=1
    C=C(N1,N2,psi);
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Shape Function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    S_u=A0_u*S(0,n,psi)+A1_u*S(1,n,psi)
    S_l=-(A0_l*S(0,n,psi)+A1_l*S(1,n,psi))
    
    # Need to make sure that the upper surface is higher than the lower surface                
    for i in range(0,len(S_u)):
        if S_u[i]>S_l[i]:
            break
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Airfoil Shape (eta=z/c)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    # Airfoil Shape (eta=z/c)
    eta_u= C*S_u+psi*deltaeta_u;
    eta_l= C*S_l-psi*deltaeta_l;
    
    # Giving back the dimensions
    y_u=c*eta_u
    y_l=c*eta_l
    return y_u,y_l

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    
    x = np.linspace(0,1)
    deltasz = [0.001,0.001]
    A0s = [0.16,0.3]
    A1s = [0.4,0.2]
    c = 1.
    y_u, y_l = CST(x, deltasz, A0s, A1s, c)
    
    plt.plot(x,y_u,'b')
    plt.plot(x,y_l,'b')
    
    A0s = [0.16,0.3]
    A1s = [0.4,0.2]
    c = 1.
    y_u, y_l = CST(x, deltasz, A0s, A1s, c)
    
    plt.plot(x,y_u)
    plt.plot(x,y_l)