#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 18:55:16 2022

@author: enricofoglia
"""

import numpy as np
from matplotlib import pyplot as plt

def read_airfoil(path):
    x, y = np.loadtxt(path, unpack = True, skiprows=1)
    return x, y    

def plot_airfoil(x,y):
    fig, ax = plt.subplots()
    ax.plot(x,y)
    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.grid(True)
    ax.set_ylim((-0.1, 0.2))
    ax.set_aspect('equal',adjustable='box') 
    
    
path = '/Users/enricofoglia/Documents/ISAE-SUPAERO1/Applied Aerodynamics II/airfoils/GOGETA.dat'
GOGx,GOGy = read_airfoil(path) 
plot_airfoil(GOGx,GOGy)

path = '/Users/enricofoglia/Documents/ISAE-SUPAERO1/Applied Aerodynamics II/airfoils/Keemera.dat'
KEEx,KEEy = read_airfoil(path) 
plot_airfoil(KEEx,KEEy)