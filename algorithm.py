# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 21:55:46 2022

@author: Asus
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 11:20:02 2022

@author: Enrico Foglia

- Every airfoil is represented as a string of binary values representing the 
percentage of merging of the initial profiles.
As a first iteration percentages will go from 0 to 100 (integers), thus requiring 
a string length of 7 for each percentage.
Given the example of a candidate library on n=4 airfoils A, B, C, D, the encoding
represents the percentage of A, then that of AB (previous merge), then that of ABC
(previuos merge) 
"""

# TODO: some plotting of the evolution of airfoils


import numpy as np
import heapq
import random
import matplotlib.pyplot as plt
from obj_fun import obj_fun
from PYfoil import *

class Airfoil: 
    def __init__(self, genome):
        self.genome = genome
    def set_fitness(self, fitness):
        self.fitness = fitness

def initialize(pop, library, n, hypP, constraints): # population, library (a directory), number of airfoils, hyperparameters, constraints 
    gen = [] # creating the generation 0
    for ind in range(pop):
        airfoil_list = []
        for i in range(n-1):
            airfoil_list.append(random.randint(0,127))
        gen.append(Airfoil(encode(airfoil_list)))
        name = 'a'+str(ind)
        Merge(library, decode(gen[-1].genome), name)
        Polar(name, iterations = 100)
        fitness = obj_fun(name, hypP, constraints)
        gen[-1].set_fitness(fitness)
    return gen
    

def reproduction(gen, pop, fitness_scaling): # 16/02/2022: modified to account for scaling
    '''
    Selects 2 parent airfoils using a "weighted roulette" approach
    The function returns the indexes of the two airfoils selected
    '''
    pop = len(gen)
    w = [gen[i].fitness for i in range(pop)]    # probabilities of survival
    if fitness_scaling != 0: # linear fitness scaling
        a, b = pre_scaling(gen, fitness_scaling)
        w = [a*w[i] + b for i in range(pop)]
    indexes = np.linspace(0, pop-1, num=pop)    # list with all the indexes: 0, 1, 2, 3, ..., pop-1
    index1, index2 = random.choices(indexes, weights=w, k = 2)
    return int(index1), int(index2)

def crossover(gen, index1, index2, pc):
    '''
    Performs crossovers at a rate defined by pc
    The function returns the modified strings
    '''
    def cross(s1, s2):
        cross_locus = random.randint(1,length-1)
        return s1[0:cross_locus] + s2[cross_locus:], s2[0:cross_locus] + s1[cross_locus:]
    length = len(gen[index1].genome)
    if random.random() <= pc:
        return cross(gen[index1].genome, gen[index2].genome)#Return the two modified strings
    return gen[index1].genome, gen[index2].genome           #Return the original strings

def elitism(gen):
    '''
    Select the two best airfoils based on their fitness values and return their genome
    '''
    pop = len(gen)
    w   = [gen[i].fitness for i in range(pop)]
    index1, index2 = heapq.nlargest(2, range(pop), w.__getitem__)
    # print(f'Best index = {index1}, {index2}')
    return gen[index1].genome, gen[index2].genome
    
def mutation(string, pm):
    '''
    Performs a mutation at a random location
    '''
    if random.random()<=pm:
        length   = len(string)
        location = random.randint(0,length-1)#
        print('Mutation position:: ', location)
        if string[location] == '0':
            string = string[:location]+'1'+string[location+1:]
        else:
            string = string[:location]+'0'+string[location+1:]
    return string

def pre_scaling(gen, C): # 16/02/2022: added 
    '''
    Generates the coefficients a and b to perform a linear fitness scaling such 
    that:
        f' = a*f + b
    respecting the constraints:
        f'_avg = f_avg
        f'_max = C*f_avg
        f'_min > 0
    '''
    pop = len(gen)
    f_avg = 0 # average fitness
    f_max = 0 # maximum fitness
    f_min = 1E6 # minimum fitness
    for i in range(pop):
        f_avg += gen[i].fitness/pop
        if gen[i].fitness > f_max:
            f_max = gen[i].fitness
        if gen[i].fitness < f_min:
            f_min = gen[i].fitness
    if f_min > (C*f_avg - f_max)/(C-1): # respecting non negativity
        delta = f_max - f_min
        a = (C-1)*f_avg/delta
        b = f_avg*(f_max - C*f_avg)/delta
    else: # max possible scaling
        delta = f_avg - f_min
        a = f_avg/delta
        b = -f_min*f_avg/delta
    return a,b 
    
def plot_data(pop,n): # 16/02/2022 modified to account for different library sizes
    '''
    Plot the Endurance over Cl for all the airfoils
    '''
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
        plt.plot(Cl, End, color = 'tab:blue')
        plt.axis([0.9, 2.0, 100, 160])
        plt.xlabel('Cl')
        plt.ylabel('Endurance')

    for i in range (n-1): # the minus one is temporary for the lack of one polar xoxo
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


def plot_iterations(max_iter, pop, fitness):
    '''
    PLOT of the fitnes function at every iteration
    '''
    plt.figure()
    mean = np.zeros(max_iter)   
    for it in range(max_iter):
        for p in range(pop):
            plt.plot(it+1, fitness[it*pop+p], 'ro', markersize=3)
            mean[it] = mean[it] + fitness[it*pop+p]
        mean[it] = mean[it]/pop
    plt.plot(range(1,max_iter+1), mean, 'ko--', markersize=3)
    plt.xlim((0, max_iter+1))
    plt.grid(which='major', axis='both')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    return



def ga(n, pop, pc, pm, library, hypP, constraints, max_iter, fitness_scaling = 0):
    # INITIALIZATION
    if pop%2 != 0:
        raise ValueError('The population size must be an even number')
    # out = open('out.txt', 'w')
    gen = initialize(pop, library, n, hypP, constraints)
        
    # Algorithm
    new_gen      = gen.copy() 
    fitness_list = np.zeros(max_iter*pop)
    
    for it in range(max_iter):
        gen     = new_gen.copy()
        # for i in range(pop):
        #     print(decode(gen[i].genome))
        new_gen = []
         
        for i in range(int((pop-2)/2)): # genetic operations
            index1, index2 = reproduction(gen, pop, fitness_scaling)# index1 and index2 are the airfoils selected
            # print(f'AfterReproduction: {index1}')
            # print(f'AfterReproduction: {index2}')
            S_child1, S_child2 = crossover(gen, index1, index2, pc)# Return the string with the modified genomes

            # print(f'AfterCrossover: {S_child1}')
            # print(f'AfterCrossover: {S_child2}')
            S_child1 = mutation(S_child1, pm)# Return the string with the modified genomes
            S_child2 = mutation(S_child2, pm)
            # print(f'AfterMutation:  {S_child1}')
            # print(f'AfterMutation:  {S_child2}')
            # pp1 = decode(S_child1)
            # pp2 = decode(S_child2)
            # print(f'Decode: {pp1}')
            # print(f'Decode: {pp2}')
            new_gen.append(Airfoil(S_child1))# New elements, given its string (genome)
            new_gen.append(Airfoil(S_child2))
            # print('\n')

        
        # elitism
        best = elitism(gen)                 # Return the strings of the best two airfoils
        new_gen.append(Airfoil(best[0]))    # New elements, given its string (genome)
        new_gen.append(Airfoil(best[1]))
        
        print(f'\n New generation: {it+1}')        
        for ind in range(pop):
            percentage = decode(new_gen[ind].genome)
            # print(f'Percentages: {percentage}')
            new_name = 'a'+str(ind)
            Merge(library, percentage, new_name)
            Polar(new_name, iterations = 100)
            new_gen[ind].fitness = obj_fun(new_name, hypP, constraints)
            fitness_list[(it)*pop+ind] = new_gen[ind].fitness    
            # print(new_gen[ind].fitness)
        
        # out.write(f'Iteration {it+1}\n')
        # fprint_results(out, new_gen)
    
    #PLOTS
    plot_iterations(max_iter, pop, fitness_list)
    plot_data(pop,n)
    # out.close()
    return new_gen


def encode(airfoil_list):
    '''
    Given a list of integers from 0 to 100 generates the related binary string
    '''
    airfoil_str = ''
    for elem in airfoil_list:
        airfoil_str += format(elem, '07b')
    return airfoil_str    
    
def decode(airfoil_str):
    '''
    Translates the airfoil representation from the binary string to a list of integers.
    '''
    n = int(len(airfoil_str)/7)
    split = [airfoil_str[i:i+7] for i in range(0, len(airfoil_str), 7)]
    decoded = []
    for i in range(n):
        decoded.append(format(int(split[i],2)/127,'.3f'))
    return decoded

def fprint_results(out, gen):
    for i in range(len(gen)):
        out.write(f'Airfoil {gen[i].genome}: {gen[i].fitness}\n')
        
if __name__ == '__main__':  # this runs only if this script is the main, thus allowing to import it in other scripts without issues
    import os
    import sys
    sys.path.insert(1, '/Users/enricofoglia/Documents/python')
    from send_message import send_message
    n           = len(os.listdir('library')) # number of airfoils in the library
    pop         = 50   # number of airfoils per generation
    pc          = 0.5  # probability of crossover
    pm          = 0.05  # probablity of mutation
    max_iter    = 50   # maximum number of generations
    library     = 'library'
    hypP        = np.array([0.1, 1, 1, 1, 1]) #Hyperparameters: [End_max, Cl_max, Delta_alpha, End_Cl_133, End_integral]
    constraints = (1.2, 3)  
    try:
        new = ga(n, pop, pc, pm, library, hypP, constraints, max_iter, fitness_scaling=1.2)
        send_message('Ehi boi, il codice ha finito di runnare, vai a darci un occhio!')
    except:
        send_message('/!\ Ohi bro, qualcosa è andato storto!\nVai a controllare asap!')

           
