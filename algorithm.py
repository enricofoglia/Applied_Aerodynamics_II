#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 11:20:02 2022

@author: Enrico Foglia

- Every airfoil is represented as a string of binary values representing the 
percentage of merging of the initial profiles.
As a first iteration percentages will go from 0 to 100 (integers), thus requiring 
a string length of 7 for each percentage (from '0000000' to '1100100').
Given the example of a candidate library on n=4 airfoils A, B, C, D, the encoding
represents the percentage of A, then that of AB (previous merge), then that of ABC
(previuos merge) 
"""

import numpy as np
import random
from obj_fun import obj_fun
from Pyfoil import *

class Airfoil: 
    def __init__(self, genome, fitness):
        self.genome = genome
        self.fitness = fitness

def initialize(n,ind):
    '''
    Generates the initial airfoils strings.

    Parameters
    ----------
    n : int
        Number of airfoils in the initial library.
    ind : int
        Index of the airfoil.

    Returns
    -------
    airfoil : string
        Binary string representing the 'pure' airfoil.

    '''
    airfoil = ''
    for i in range(n-1):
        if i == ind:
            airfoil += '1'*7 # binary for 127
        else:
            airfoil += '0'*7
    return airfoil

def reproduction(gen):
    '''
    Selects 2 parent airfoils using a "weighted roulette" approach
    '''
    pop = len(gen)
    w = [gen[i].fitness for i in range(pop)] # probabilities of survival
    parent1, parent2 = random.choices(gen, weights=w, k = 2)
    return parent1, parent2

def crossover(parent1, parent2, pc):
    '''
    Performs crossovers at a rate defined by pc
    '''
    def cross(s1, s2):
        cross_locus = random.randint(1,length-1)
        return s1[0:cross_locus] + s2[cross_locus:], s2[0:cross_locus] + s1[cross_locus:]
    length = len(parent1.genome)
    if random.random() <= pc:
        parent1.genome, parent2.genome = cross(parent1.genome, parent2.genome)
    return parent1, parent2
        
def mutation(individual, pm):
    '''
    Performs a mutation at a random location
    '''
    if random.random()<=pm:
        length = len(individual[0])
        location = random.randint(0,length)
        if individual.genome[location] == '0':
            individual.genome[location] = '1'
        else:
            individual.genome[location] = '0'
    return individual

def fitness(gen):
    '''
    Will calculate the fitness of an airfoil and return it.
    In practice, given the airfoil string, will create the .dat using the merge
    funciton of PYfoil file and have the obj_fun analize it
    Will also calculate the probability of survival.
    
    Parameters
    ----------
    gen : list (of lists)
        list of binary strings of the airfoils and their respective fitness values.

    Returns
    -------
    None.

    '''
    pass

def ga(n, pop, pc, pm, max_iter, library, hypP, constraints):
    # INITIALIZATION
    if pop%2 != 0:
        raise ValueError('The population size must be an even number')
        
    gen = []
    w = random.choices(range(500), k = 10) # example weights
    for ind in range(pop):
        airfoil_list = []
        for i in range(n-1):
            airfoil_list.append(random.randint(0,127))
        gen.append(Airfoil(encode(airfoil_list), w[ind]))
    # Algorithm
    new_gen = gen
    for it in range(max_iter):
        gen = new_gen
        new_gen = []
        for i in range(int(pop/2)): # genetic operations
            parent1, parent2 = reproduction(gen)
            child1, child2 = crossover(parent1, parent2, pc)
            child1 = mutation(child1, pm)
            child2 = mutation(child2, pm)
            new_gen.append(child1)
            new_gen.append(child2)
        for ind in range(pop):
            percentage = decode(new_gen[ind].genome)
            new_name = 'a'+str(ind)
            Merge(library, percentage, new_name)
            Polar(new_name, iterations = 500)
            # new_gen[ind].fitness = obj_fun(new_name, hypP, constraints)
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
    n = int(len(airfoil_str)/7 + 1 )
    split = [airfoil_str[i:i+7] for i in range(0, len(airfoil_str), 7)]
    decoded = []
    for i in range(n-1):
        decoded.append(format(int(split[i],2)/127,'.3f'))
    return decoded


if __name__ == '__main__': # this runs only if this script is the main, thus allowing to import it in other scripts without issues
    n = 2 # number of airfoils in the library
    pop = 4 # number of airfoils per generation
    pc = 0.5 # probability of crossover
    pm = 0.001 # probablity of mutation
    max_iter = 2 # maximum number of generations
    library = 'library'
    hypP    = np.multiply([1, 1, 1, 1, 1], 10)      #Hyperparameters, to be set manually (trial and error)
    constraints = (1.2, 6.0)  
    new = ga(n, pop, pc, pm, max_iter, library, hypP, constraints)







