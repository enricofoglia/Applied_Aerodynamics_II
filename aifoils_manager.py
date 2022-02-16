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


    
# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/generation/a0.dat'
# x,y = read_airfoil(path) 
# plot_airfoil(x,y)

# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/generation/a1.dat'
# x1,y1 = read_airfoil(path) 
# plot_airfoil(x,y)

# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/generation/a2.dat'
# x,y = read_airfoil(path) 
# plot_airfoil(x,y)

# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/generation/a3.dat'
# x,y = read_airfoil(path) 
# plot_airfoil(x,y)

# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/generation/a6.dat'
# x_,y_ = read_airfoil(path) 
# plot_airfoil(x_,y_)

# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/library/AH79100C.dat'
# xe,ye = read_airfoil(path) 
# plot_airfoil(xe,ye)


# path = '/Users/enricofoglia/Documents/git/Applied_Aerodynamics_II/library/e216.dat'
# xd,yd = read_airfoil(path) 
# plot_airfoil(xd,yd)

