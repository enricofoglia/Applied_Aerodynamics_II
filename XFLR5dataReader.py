#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 13:34:13 2022

@author: enricofoglia
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def Op_graphics(data, params):
    x,Cpi,Cpv,Qi,Qv = np.split(data,5,axis = 1)
    name, alpha,Re,M,Ncrit = params
    figure1, axs = plt.subplots(nrows=1, ncols=2, constrained_layout=True)
    axs[0].plot(x,-Cpi)
    axs[0].set_xlabel('x/c')
    axs[0].set_ylabel('-Cp')
    axs[0].set_title('Inviscid pressure coefficient')
    axs[0].grid(True)
    axs[1].plot(x,-Cpi)
    axs[1].set_xlabel('x/c')
    axs[1].set_ylabel('-Cp')
    axs[1].set_title('Viscous pressure coefficient')
    axs[1].grid(True)
    figure1.suptitle(f'{name}: alpha = {alpha}; Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    figure2, axs = plt.subplots(nrows=1, ncols=2, constrained_layout=True)
    axs[0].plot(x,Qi)
    axs[0].set_xlabel('x/c')
    axs[0].set_ylabel('Q')
    axs[0].set_title('Inviscid velocity distribution')
    axs[0].grid(True)
    axs[1].plot(x,Qv)
    axs[1].set_xlabel('x/c')
    axs[1].set_ylabel('Q')
    axs[1].set_title('Viscous velocity distribution')
    axs[1].grid(True)
    figure2.suptitle(f'{name}: alpha ={alpha}; Re ={Re}; Mach ={M}; Ncrit ={Ncrit}')
    
def BL_graphics(data, params):
    x,Hk,Ue,Cf,Cd,A,D_star,theta,CTq = np.split(data, 9, axis = 1)
    name, alpha,Re,M,Ncrit = params
    # diplacement and momentum thickness, shape parameter
    fig1, axs = plt.subplots(3,1,constrained_layout=True)
    fig1.suptitle(f'{name}: alpha = {alpha}; Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    axs[0].plot(x[x<=1],D_star[x<=1])
    axs[1].plot(x[x<=1],theta[x<=1])
    axs[2].plot(x[x<=1],Hk[x<=1])
    for ax in axs:
        ax.set_xlabel('x/c')
        ax.grid(True)
    axs[0].set_ylabel(r'$\delta^*$')
    axs[1].set_ylabel(r'$\vartheta$')
    axs[2].set_ylabel(r'$H_k$')
    axs[0].set_title('Displacement thikness')
    axs[1].set_title('Momentum thickness')
    axs[2].set_title('Kinematic shape parameter')
    
    # Drag coefficients
    fig2, axs = plt.subplots(2,1,constrained_layout=True)
    fig2.suptitle(f'{name}: alpha = {alpha}; Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    axs[0].plot(x[x<=1],Cd[x<=1])
    axs[1].plot(x[x<=1],Cf[x<=1])
    for ax in axs:
        ax.set_xlabel('x/c')
        ax.grid(True)
    axs[0].set_ylabel(r'$C_d$')
    axs[1].set_ylabel(r'$C_f$')
    axs[0].set_title('Drag coefficient')
    axs[1].set_title('Skin friction coefficient')
    
    # External velocity Ue
    fig3, ax = plt.subplots(1,1)
    fig3.suptitle(f'{name}: alpha = {alpha}; Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    ax.plot(x[x<=1],Ue[x<=1])
    ax.set_xlabel('x/c')
    ax.set_ylabel(r'$U_e$')
    ax.set_title('External velocity')
    ax.grid(True)
    
    # Amplification of disturbances A/A0
    fig4, ax = plt.subplots(1,1)
    fig4.suptitle(f'{name}: alpha = {alpha}; Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')
    ax.plot(x[x<=1],A[x<=1])
    ax.set_xlabel('x/c')
    ax.set_ylabel(r'$A/A_0$')
    ax.set_title('Amplification of disturbances')
    ax.grid(True)
    
def polar_graphics(data, params):
    name, xtrf_top, xtrf_bottom, M, Re, Ncrit = params
    alpha, CL, CD, CDp, Cm, xtr_top, xtr_bottom, Cpmin, Chinge, XCp = np.split(data, 10, axis = 1)
    # CL and Cm
    fig1, ax = plt.subplots()
    color1 = 'royalblue'
    ax.plot(alpha, CL, color = color1)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$C_L$', color = color1)
    ax.set_title('Lift and moment coefficients')
    ax.grid(True)
    ax_ = ax.twinx()
    color2 = 'navy'
    ax_.plot(alpha, Cm, color = color2)
    ax_.set_ylabel(r'$C_m$', color = color2)
    ax_.set_ylim((-0.5,0))
    fig1.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')

    fig1.tight_layout()

    # polar
    fig2, ax = plt.subplots()
    ax.plot(CD, CL)
    ax.set_xlabel(r'$C_D$')
    ax.set_ylabel(r'$C_L$')
    ax.set_title('Polar')
    ax.grid(True)
    fig2.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')

    # Efficiency 
    fig3, ax = plt.subplots()
    ax.plot(alpha, CL**1.5/CD)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$C_L^{1.5}/C_D$')
    ax.set_title('Efficiency parameter')
    ax.grid(True)
    fig3.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')

    # transition
    fig4, ax = plt.subplots()
    ax.plot(alpha, xtr_top, color = color1)
    ax.plot(alpha, xtr_bottom, color = color2)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$x_{tr}/c$')
    ax.set_title('Transition position')
    ax.grid(True)
    ax.legend(['top', 'bottom'])
    fig4.suptitle(f'{name}:  Re = {Re}; Mach = {M}; Ncrit = {Ncrit}')

    
def read_OpPoint(path, graphics = True):
    df = pd.read_csv(path,nrows = 3, header=None)
    name = df.iloc[1,0]
    
    df = pd.read_csv(path,nrows = 1, header=None, skiprows=3)
    alpha = df.iloc[0,1]
    Re = df.iloc[0,3]
    M = df.iloc[0,5]
    Ncrit = df.iloc[0,7]
    
    df = pd.read_csv(path, skiprows=4)
    data = df.to_numpy()
    
    if graphics == True:
        Op_graphics(data,(name,alpha,Re,M,Ncrit))
    return data

def read_BL(path, graphics_top = True, graphics_bottom = False):
    df = pd.read_csv(path,nrows = 2, header=None)
    name = df.iloc[1,0]
    
    df = pd.read_csv(path,nrows = 1, header=None, skiprows=2)
    alpha = df.iloc[0,1]
    Re = df.iloc[0,3]
    M = df.iloc[0,5]
    Ncrit = df.iloc[0,7]
    
    # the solution used here should be improved
    df = pd.read_csv(path, skiprows=6) # there are some blank rows
    data = df.to_numpy()
    bottom_start = np.where(data=='Bottom Side')[0][0]
    # df_top = pd.read_csv(path,skiprows=6, nrows = bottom_start)
    # df_bottom = pd.read_csv(path,skiprows= 10 + bottom_start) # stupid blank lines
    df_top = df.loc[:bottom_start-1,:]
    df_bottom = df.loc[bottom_start+2:,:]
    data_top = df_top.to_numpy(dtype = np.float64)
    data_bottom = df_bottom.to_numpy(dtype = np.float64)
    if graphics_top == True:
        BL_graphics(data_top, (name,alpha,Re,M,Ncrit))
    if graphics_bottom == True:
        BL_graphics(data_bottom, (name,alpha,Re,M,Ncrit))        
    return data_top, data_bottom

def read_polar(path, graphics = True):
    df = pd.read_csv(path,nrows = 2, header=None)
    name = df.iloc[1,0].replace('Calculated polar for: ','')
    
    df = pd.read_csv(path, skiprows=5, nrows = 2, header=None)
    xtrf_top, xtrf_bottom = df.iloc[0,0].split()[2:5:2]
    M = df.iloc[1,0].split()[2]
    Re = ''.join(df.iloc[1,0].split()[5:8])
    Ncrit = df.iloc[1,0].split()[-1]
    
    df = pd.read_csv(path, skiprows=8)
    data = df.to_numpy(dtype = np.float64)
    if graphics == True:
        polar_graphics(data,(name, xtrf_top, xtrf_bottom, M, Re, Ncrit))
    return data
    
    
path = 'data/AH79100C_OP_data_a5.csv'
data = read_OpPoint(path, graphics = True)
pathBL = 'data/AH79100C_BL_data_a5.csv'
dataBL_top, dataBL_bottom = read_BL(pathBL, graphics_top=True)
pathPolar = 'data/AH79100C_polar.csv'
polar = read_polar(pathPolar)


