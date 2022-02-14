# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 16:40:45 2022

@author: Alberto Zigante
"""

import numpy as np 
import matplotlib.pyplot as plt


# Object function
def plot_data(pop):
    
    plt.figure()
    plt.grid(which='major', axis='both')
    for i in range(pop):
            
        name = 'data/a'+str(i)+'.log' # .txt for XFLR5 outputs
        acquisition = np.loadtxt(name, skiprows=12) #skiprows=11 for XFLR5 outputs, skiprows=12 for XFOIL outputs
        Cl      = acquisition[:,1]
        Cd      = acquisition[:,2]
    
        End = np.zeros(len(Cl))
        for i in range(len(Cl)-1):
                End[i] = Cl[i]**1.5/Cd[i]
       
        plt.plot(Cl, End)
        plt.axis([0.9, 2.0, 100, 160])
        #plt.grid(which='major', axis='both')
        plt.xlabel('Cl')
        plt.ylabel('Endurance')

    
    for i in range (4):
        name = 'reference_data/a_ref_'+str(i)+'.txt' # .txt for XFLR5 outputs
        acquisition = np.loadtxt(name, skiprows=11) #skiprows=11 for XFLR5 outputs, skiprows=12 for XFOIL outputs
        Cl      = acquisition[:,1]
        Cd      = acquisition[:,2]
    
        End = np.zeros(len(Cl))
        for i in range(len(Cl)-1):
                End[i] = Cl[i]**1.5/Cd[i]
        
        plt.plot(Cl, End, 'k--')
        plt.axis([0.9, 2.0, 100, 160])
        plt.xlabel('Cl')
        plt.ylabel('Endurance')

         
    return
