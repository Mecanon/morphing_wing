# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:31 2016

@author: endryws
"""

import pylab as plt
import math

class Winglet(object):
    def __init__(self, posicoes, R, k, area=None, epss_0=None,
                 tamanho_sem_deformacao_s=None):
                     
     """                
     - posicoes e um dicionario com as posicoes da mola linear e da Sma    
     o n significa o ponto fixo e o p, a outra ponta da mola, este ponto
     varia de acordo com a deformacao epsilon sofrida pela mesma
     - epss_0 é a deformação inicial(épsilon) da mola SMA; se não for 
     definida ela é igual a zero.
     - k e o coeficiente elástico da mola linear.
     """
     
    # Guardando as posicoes para o sistema de coordenadas
     self.xl_n = posicoes['xl_n']
     self.xl_p = posicoes['xl_p']
     self.xs_n = posicoes['xs_n']
     self.xs_p = posicoes['xs_p']
        
    # Encontrando os tamanhos iniciais das molas
        
     self.rl_0 = self.xl_p - self.xl_n  # Jogada para receber 
                                        # tamanho > 0
     self.rs_0 = self.xs_p - self.xs_n  # Jogada para receber 
                                        # tamanho > 0
    
    # Iniciando o tamanho das molas com os valores iniciais, os 
    # definidos acima.
    # Estes valores serão alterados, quando a polia girar.
    
     self.rl = self.rl_0
     self.rs = self.rs_0
        
    #Definindo os valores iniciais para o raio, o angulo e a força.
     self.R = R        
     self.theta_0 = 0
     self.theta = self.theta_0
     self.F = 0.
    
    # Definindo informacoes importantes para o cálculo da forca exercida
    # pela SMA
    
     self.area = area
     self.sigma = None
     self.k = k
    
    # tamanho_sem_deformacao_s é o tamanho original da SMA.
    # No caso de tamanho_sem_deformacao_s nao ser definido, mas a tensão
    # inicial é conhecida, calcula a variavel tamanho_sem_deformacao.
    # Se tamanho_sem_deformacao for definido, calcula a deformacao
    # inicial (eps_0).
    
     if tamanho_sem_deformacao_s == None and epss_0 != None:
         
         self.epss_0 = epss_0
         self.epss = self.epss_0
            
    # Deve ser fornecida a deformacao sofrida pela SMA.
            
         self.tamanho_sem_deformacao_s = self.rs_0/(1 + self.epss)
            
    # Neste caso, como ambas as molas estão conectadas, a compressão da
    # SMA significa um elongamento da linear e vice-versa; por isso
    # utilizo (1 - deformação da SMA).
            
         self.tamanho_sem_deformacao_l = self.rl_0/(1 - self.epss)
            
     elif tamanho_sem_deformacao_s != None and epss_0 == None:
            
         self.tamanho_sem_deformacao_s = self.rs_0
         self.tamanho_sem_deformacao_l = self.rl_0
            
         self.epss_0 = (self.rs / self.tamanho_sem_deformacao_s) - 1
         self.epss = self.epss_0
            
    def update(self,theta):
        
    # math.radius devolve o angulo entregue em graus, em radianos, e 
    # pode retornar valores positivos ou negativos.
    # Por conta disto, delta_r é positivo para angulos positivos, e 
    # negativo para angulos negativos.
    # theta > 0 para angulos que vão de x+ para y+.
        
        self.theta = theta
        
        delta_r = math.radians(self.theta)*self.R
         
        self.rs = self.xl_p - self.xl_n - delta_r
        self.rl = self.xs_p - self.xs_n + delta_r
        
        self.epss = (self.rs / self.tamanho_sem_deformacao_s) - 1
        self.epsl = (self.rl / self.tamanho_sem_deformacao_l) - 1
        
    def calcula_forca(self,material, tipo='tensao'):
        
        if tipo == 'tensao':
            if material == 'linear':
                self.F = self.k*(self.epsl*self.rl)
            elif material == 'SMA':
                print "Colocar o modelo constitutivo para o SMA"         
                
    # Calcula a força de tensão através da secao transversal.
                
        elif tipo == 'sigma':
            self.F = self.area * self.sigma      
    
    def calcula_torque(self):
    # Calcula o torque dada a forca do atuador.
    # Não me parece necessária, mas mantenho com as devidas alterações; 
    # para futuramente unir os codigos.
    # Como ambas as forças Fs e Fl tem um braco de alavanca R.
        self.torque = self.F*self.R
        
        return self.torque
    
    def imprime_modelo(self,cor):

    # Cor é uma variável que é utilizada para diferenciar cada uma 
    # das vezes que é impresso o modelo.
    # k = preto ; r = vermelho ; b = azul ; y = amarelo              
    
        sup = self.R + self.R / 4
        inf = self.R - self.R / 4
        tamanho_l = self.rl
        tamanho_s = self.rs
        anel_l = tamanho_l / 10
        anel_s = tamanho_s / 10
        
        cir1 = plt.Circle((0,0), radius= self.R, alpha =.7, fc='k')
        cir2 = plt.Circle((0,0), radius=(self.R/2), alpha =.5, fc='w')
                                                        
        ax = plt.axes(aspect=1) # Cria eixos vazios (aspect=1 tem a ver
                                # com a escala)

        ax.add_patch(cir1)                    
        ax.add_patch(cir2) 

        plt.plot([self.xl_n,self.xl_n + anel_l,self.xl_n + 2*anel_l,
                  self.xl_n + 3*anel_l, self.xl_n + 4*anel_l, 
                  self.xl_n + 5*anel_l, self.xl_n + 6*anel_l,
                  self.xl_n + 7*anel_l, self.xl_n + 8*anel_l, 
                  self.xl_n + 9*anel_l, self.xl_n + 10*anel_l,0], 
                  [-self.R, -sup, -inf, -sup, -inf, -sup, -inf, -sup,
                   -inf, -sup, -self.R, -self.R], cor )

        plt.plot([self.xs_n, self.xs_n + anel_s, self.xs_n + 2*anel_s,
                  self.xs_n + 3*anel_s, self.xs_n + 4*anel_s, 
                  self.xs_n + 5*anel_s, self.xs_n + 6*anel_s,
                  self.xs_n + 7*anel_s, self.xs_n + 8*anel_s, 
                  self.xs_n + 9*anel_s, self.xs_n + 10*anel_s,0], 
                  [self.R,sup,inf,sup,inf,sup,inf,sup,
                   inf,sup,self.R,self.R], cor )
          
        plt.plot([0, self.R + 0.1],[0,0],'b')
        
        if self.theta == 0:
            pass
        else:
            plt.plot([0, self.R*math.cos(math.radians(self.theta))],
                      [0, self.R*math.sin(math.radians(self.theta))],
            'r')
        
    
    #Definindo a operação imprimir / print (utilizada para teste)
    
    def __str__(self, ponto):
       return self.posicoes[ponto]
       
########################################################################
# Aqui comeca a utilizacao da classe para o winglet.
########################################################################
 
if __name__ == '__main__':
    
#=======================================================================
# Variaveis do projeto    
#=======================================================================
    
    posicoes = {'xl_n': -1, 'xl_p': -0.5, 'xs_n': -1, 'xs_p': -0.5}
    
    radius = 0.2
    
    k = 0.5
    
    deformacao_inicial = 0.2
    
#=======================================================================
# Criando uma instância de Winglet
#=======================================================================

    wing = Winglet(posicoes,radius,k,epss_0 = deformacao_inicial)

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