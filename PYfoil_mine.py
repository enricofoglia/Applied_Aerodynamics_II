#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 22:04:02 2022

@author: Enrico Foglia

This script allows to call Xfoil directly to perform different types of 
analysis. In particular the first version allows to calculate the polar for a 
given range of alpha and to calculate the Cp distribution and Boundary Layer 
quantities for a given alpha. Adapted vesions of the plotting routines of 
XFLR5dataReader are given if plotting is desired.
"""

import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt

#%% Functions to run Xfoil

def calculate_polar(airfoil,
                    Re = 350e3,
                    Ncrit = 9,
                    Xtr_top = 1, Xtr_bot = 1,
                    iterations = 100,
                    a_min = -5, 
                    a_max = 15, 
                    delta_a = 0.25):
    # write input file to run xfoil
    with open('input_file.in', 'w') as f:
        f.write(f'load {airfoil}.dat\n')
        f.write('pane\n')
        f.write(f'oper\nv {Re}\n')
        f.write(f'vpar\nN {Ncrit}\nXtr {Xtr_top} {Xtr_bot}\n\n')
        f.write(f'iter {iterations}\n')
        f.write('pacc\npolar.txt\n\n')
        f.write(f'aseq {a_min} {a_max} {delta_a}\n\n')
        f.write('quit\n')
    # delete polar file if already existing
    if os.path.exists('polar.txt'):
        os.remove('polar.txt')
    # run xfoil
    subprocess.run('Xfoil/bin/xfoil < input_file.in', shell = True)
    # save polar as np.array
    polar_data = np.loadtxt('polar.txt', skiprows=12)
    return polar_data

def calculate_BL(airfoil,
                 Re = 350000,
                 Ncrit = 9,
                 Xtr_top = 1, Xtr_bot = 1,
                 iterations = 100,
                 alpha = 5):
    with open('input_file.in', 'w') as f:
        f.write(f'load {airfoil}.dat\n')
        f.write('pane\n')
        f.write(f'oper\nv {Re}\n')
        f.write(f'vpar\nN {Ncrit}\nXtr {Xtr_top} {Xtr_bot}\n \n')
        f.write(f'iter {iterations}\n')
        f.write(f'a {alpha}\n')
        f.write('dump BL.txt\n\n')
        f.write('quit\n')
    # delete polar file if already existing
    if os.path.exists('BL.txt'):
        os.remove('BL.txt')
    # run xfoil
    subprocess.run('Xfoil/bin/xfoil < input_file.in', shell = True)
    # save polar as np.array
    BL_data = np.loadtxt('BL.txt', skiprows=1, max_rows=160)
    # TODO: acquire N and CD data
    return BL_data
#%% Miscellanea of useful functions

def max3(a,b,c):
    return np.maximum(a,np.maximum(b,c))

def min3(a,b,c):
    return np.minimum(a,np.minimum(b,c))

def split_top_bottom(x):
    ind = 0
    while x[ind] >= x[ind+1]:
        ind += 1
    return ind+1

#%% Plotting functions

def polar_graphics(data, params = (), data_min = np.array([]), data_max = np.array([])):
    
    alpha, CL, CD, CDp, Cm, xtr_top, xtr_bottom, _, _ = np.split(data, 9, axis = 1) # modified w.r.t XFLR5dataReader
    # CL and Cm
    fig1, ax1 = plt.subplots()
    color1 = 'tab:blue'
    ax1.plot(alpha, CL, color = color1)
    ax1.set_xlabel(r'$\alpha$')
    ax1.set_ylabel(r'$C_L$', color = color1)
    ax1.set_title('Lift and moment coefficients')
    ax1.grid(True)
    ax_1 = ax1.twinx()
    color2 = 'tab:red'
    ax_1.plot(alpha, Cm, color = color2)
    ax_1.set_ylabel(r'$C_m$', color = color2)
    ax_1.set_ylim((-0.4,0.2))
    fig1.tight_layout()
    
    # polar
    fig2, ax2 = plt.subplots()
    ax2.plot(CD, CL)
    ax2.set_xlabel(r'$C_D$')
    ax2.set_ylabel(r'$C_L$')
    ax2.set_title('Polar')
    ax2.grid(True)

    # Efficiency 
    fig3, ax3 = plt.subplots()
    ax3.plot(alpha, CL**1.5/CD)
    ax3.set_xlabel(r'$\alpha$')
    ax3.set_ylabel(r'$C_L^{1.5}/C_D$')
    ax3.set_title('Efficiency parameter')
    ax3.grid(True)

    # transition
    fig4, ax4 = plt.subplots()
    ax4.plot(alpha, xtr_top, color = color1)
    ax4.plot(alpha, xtr_bottom, color = color2)
    ax4.set_xlabel(r'$\alpha$')
    ax4.set_ylabel(r'$x_{tr}/c$')
    ax4.set_title('Transition position')
    ax4.grid(True)
    ax4.legend(['top', 'bottom'])
    
    if params != (): # 24/01/2022:possibility of not showing title
        name, xtrf_top, xtrf_bottom, M, Re, Ncrit = params
        fig1.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
        fig2.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
        fig3.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
        fig4.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    
    if data_min.any() != False: # 24/01/2022: possibility to plot the errorbars for different Ncrit
        alpha_max, CL_max, CD_max, CDp_max, Cm_max, xtr_top_max, xtr_bottom_max, Cpmin_max, Chinge_max, XCp_max = np.split(data_max, 10, axis = 1)
        alpha_min, CL_min, CD_min, CDp_min, Cm_min, xtr_top_min, xtr_bottom_min, Cpmin_min, Chinge_min, XCp_min = np.split(data_min, 10, axis = 1)
        ax1.fill_between(alpha[:,0], min3(CL[:,0],CL_min[:,0],CL_max[:,0]), 
                         max3(CL[:,0],CL_min[:,0],CL_max[:,0]), alpha = 0.3, color = color1)
        ax_1.fill_between(alpha[:,0], max3(Cm[:,0],Cm_min[:,0],Cm_max[:,0]), 
                          min3(CL[:,0],Cm_min[:,0],Cm_max[:,0]), alpha = 0.3, color = color2)
        ax3.fill_between(alpha[:,0], min3(CL[:,0]**1.5/CD[:,0],CL_min[:,0]**1.5/CD_min[:,0],CL_max[:,0]**1.5/CD_max[:,0]), 
                         max3(CL[:,0]**1.5/CD[:,0],CL_min[:,0]**1.5/CD_min[:,0],CL_max[:,0]**1.5/CD_max[:,0]), alpha = 0.3)
        ax4.fill_between(alpha[:,0], min3(xtr_top[:,0], xtr_top_min[:,0], xtr_top_max[:,0]), 
                         max3(xtr_top[:,0], xtr_top_min[:,0], xtr_top_max[:,0]), color = color1, alpha = 0.3)
        ax4.fill_between(alpha[:,0],min3(xtr_bottom[:,0], xtr_bottom_min[:,0], xtr_bottom_max[:,0]), 
                         max3(xtr_bottom[:,0], xtr_bottom_min[:,0], xtr_bottom_max[:,0]), color = color2, alpha = 0.3)
        # the polar graph does not come out well, it is omitted
        
def BL_graphics(data, params = ()):
    s, x, y, Ue, D_star, theta, Cf, H, H_star, P, m, K = np.split(data, 12, axis = 1)
    color1 = 'tab:blue'
    color2 = 'tab:red'
    ind = split_top_bottom(x[x<=1])
    # diplacement and momentum thickness, shape parameter
    fig1, axs = plt.subplots(3,1,constrained_layout=True)
    axs[0].plot(x[x<=1][ind-1:], D_star[x<=1][ind-1:], color = color2)
    axs[0].plot(x[:ind], D_star[:ind], color = color1) 
    axs[1].plot(x[x<=1][ind-1:], theta[x<=1][ind-1:], color = color2)
    axs[1].plot(x[:ind], theta[:ind], color = color1)
    axs[2].plot(x[x<=1][ind-1:],H[x<=1][ind-1:], color = color2)
    axs[2].plot(x[:ind],H[:ind], color = color1)

    for ax in axs:
        ax.set_xlabel('x/c')
        ax.legend(['bottom', 'top'])
        ax.grid(True)
    axs[0].set_ylabel(r'$\delta^*$')
    axs[1].set_ylabel(r'$\vartheta$')
    axs[2].set_ylabel(r'$H_k$')
    axs[0].set_title('Displacement thikness')
    axs[1].set_title('Momentum thickness')
    axs[2].set_title('Kinematic shape parameter')
    
    # Skin friction coefficients
    fig2, axs = plt.subplots()
    axs.plot(x[x<=1][ind-1:],Cf[x<=1][ind-1:], color = color2)
    axs.plot(x[:ind],Cf[:ind], color = color1)
    ax.set_xlabel('x/c')
    axs.set_ylabel(r'$C_f$')
    axs.set_title('Skin friction coefficient')
    axs.legend(['bottom', 'top'])
    axs.grid(True)
    
    # External velocity Ue
    fig3, ax = plt.subplots(1,1)
    ax.plot(x[x<=1][ind-1:],Ue[x<=1][ind-1:], color = color2)
    ax.plot(x[:ind],Ue[:ind], color = color1)
    ax.set_xlabel('x/c')
    ax.set_ylabel(r'$U_e$')
    ax.set_title('External velocity')
    ax.legend(['bottom', 'top'])
    ax.grid(True)
    
    if params != (): # 24/01/2022:possibility of not showing title
        name, alpha,Re,M,Ncrit = params
        fig1.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
        fig2.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
        fig3.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    

#%% Run file (test)

polar_data = calculate_polar('e387')
polar_graphics(polar_data)
BL_data = calculate_BL('e387')
BL_graphics(BL_data)