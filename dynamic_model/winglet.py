# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:31 2016

@author: endryws
"""


from actuator import actuator


if __name__ == '__main__':
    
#=======================================================================
# Variaveis do projeto    
#=======================================================================
    
    posicoes = {'xl_n': -1, 'xl_p': -0.5, 'xs_n': -1, 'xs_p': -0.5}
    
    radius = 0.2
    
    k = 0.5
    
    deformacao_inicial = 0.2
    
#=======================================================================
# Criando uma instância de actuator
#=======================================================================

    wing = actuator(posicoes,radius,k,epss_0 = deformacao_inicial)

#=======================================================================
# Criando um desenho inicial do modelo
#=======================================================================

    wing.imprime_modelo('r')
    
#=======================================================================
# tendo em mãos um ângulo, calculo a deformacao que foi gerada.
# calculo a força gerada e o torque gerados a partir deste angulo que 
# foi girado.
#=======================================================================
    angulo_de_rotacao = 90 
    
    material_l = 'linear'
    
    material_SMA= 'SMA'
    
    wing.update(angulo_de_rotacao)
    
    wing.calcula_forca(material_l)
    
    wing.calcula_forca(material_SMA)
    
    wing.calcula_torque()

    

#=======================================================================
# Redesenha o modelo
#=======================================================================

    wing.imprime_modelo('b')
    
#=======================================================================
# 
#=======================================================================