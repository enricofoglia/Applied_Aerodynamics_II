# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 17:05:16 2022

@author: Alberto Zigante
"""


import matplotlib.pyplot as plt


# Object function
def plot_iterations(max_iter, pop, fitness):
        
    #PLOT of the fitnes per iteration
    plt.figure()
        
    for it in range(max_iter):
        for p in range(pop):
            plt.plot(it+1, fitness[it*pop+p], 'ro', markersize=1)
    plt.xlim((0, max_iter+1))
    plt.grid(which='major', axis='both')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    return


#This function doesn't show with different size the number of airfoils of the same type,
#partially because we don't have anymore same airfoils, and because the old function was wrong