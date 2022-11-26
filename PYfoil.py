# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 16:23:10 2022

@author: Aldo Schioppa, Enrico Foglia
"""

import subprocess as sp
import os
import numpy as np
import matplotlib.pyplot as plt


def Merge(directory,
          percentage,
          newname,
          ):
    '''
    Calls XFOIL in order to merge all airfoils in the directory according to the
    percentages specified in 'percentage'. The new airfoil is saved in the directory 
    'geeration/' as a .dat file with the name specified in 'newname'. 
    '''
    airfoils = os.listdir(directory)
    with open('input_file.in', 'w') as f:
        f.write('\nload library/' + airfoils[0] + '\n')
        for i in range(len(percentage)):
            f.write('inte\n')
            f.write('C\n')
            f.write('F\n')
            f.write('library/' + airfoils[i+1] + '\n') # 
            interpolating_fraction = percentage[i]
            f.write(str(interpolating_fraction) + '\n')
            f.write('generation/' + newname + '\n')# airfoil nam, to upload
            f.write('PANE\n')
        f.write('save\n')
        f.write('generation/' + newname + '.dat\n')
        f.write('quit')
    if os.path.exists('generation/'+newname+'.dat'):
        os.remove('generation/'+newname+'.dat')
    sp.run('Xfoil/bin/xfoil < input_file.in', shell=True, capture_output = True)  # WANRNING: has to be adapted to the operating system

def create_input_polar(airfoil,
          iterations = 100,
          Re = 350e3,
          M = 0,
          x_tr_top = 1.0,
          x_tr_bot = 1.0,
          N_crit = 9.0,
          alpha_min = -2.0,
          alpha_max = 15.0,
          step = 0.5): # 16/02/2022: moved outside to be able to call it multiple times
    '''
    Write the input file for XFOIL for computing the polar of the airfoil.
    '''
    with open('input_file.in', 'w') as f:
        # f.write(' \n')
        # f.write('plop\n')
        # f.write('g\n')
        f.write(' \n')
        f.write('load generation/' + airfoil + '.dat\n') # First thing Xfoil is gonna ask the user is the profile.
        f.write('PANE\n')
        f.write('OPER\n')
        f.write('Vpar\n') # managing transition abscissa, N_crit, x_tr and further random parameters we don't really care of
        f.write('N ' + str(N_crit) + '\n')
        f.write('xtr\n') # typed 'xtr' in Vpar, Xfoil firstly ask to enter xtr_top/c.
        f.write(str(x_tr_top) + '\n')
        f.write(str(x_tr_bot) + '\n') # Observe, if you want them to be unvaried (both at 100% of the chord), let the function pass them as 1
        f.write(' \n')
        f.write(f'iter {iterations}\n') 
        f.write('visc ' + str(Re) + '\n')
        f.write('Mach\n')
        f.write(str(M) + '\n') # I split the command into to 2 steps, because apparently putting them together does not work out.
        f.write('PACC\n') # initialise 1 polar file for data accumulation (it goes without saying, for the input airfoil)
        f.write('data/' + airfoil + '.log\n')  # output file, .log stands for 'test document file'.
        f.write(' \n')          # after, Xfoil will ask you to enter a name for a 'Dump file', which I don't know what is so no dump file
        f.write(f'aseq {alpha_min} {alpha_max} {step}\n') # command for Airfoil, concerning AoA sequences
        f.write(' \n')     # escape OPER
        f.write('quit\n')  # exit

def Polar(airfoil,
          iterations = 100,
          Re = 350e3,
          M = 0,
          x_tr_top = 1.0,
          x_tr_bot = 1.0,
          N_crit = 9.0,
          alpha_min = -2.0,
          alpha_max = 15.0,
          step = 0.5
          ):
    '''
    Compute the polar of an airfoil and save the results as a .txt file 
    in the 'data/' directory. If the convergence fails in the prescribed 
    number of iterations, the routine is tried again three times increasing
    the number of iterations before aborting.
    '''
        # Deleate the polar if already existing
    if os.path.exists('data/'+airfoil + '.log'):
        os.remove('data/'+airfoil + '.log')
        # Run Xfoil  
    create_input_polar(airfoil = airfoil,
                 iterations = iterations,
                 Re = Re,
                 M = M, 
                 x_tr_top= x_tr_top,
                 x_tr_bot = x_tr_bot,
                 N_crit = N_crit,  
                 alpha_min = alpha_min,
                 alpha_max = alpha_max,
                 step = step)
    out_message = sp.run('Xfoil/bin/xfoil < input_file.in', shell=True, capture_output = True) # WANRNING: has to be adapted to the operating system
    error = 'Sequence halted since previous  4 points did not converge'
    CONV_FLAG = 1
    while error in str(out_message.stdout) and CONV_FLAG <= 3: # 16/02/2022: try to add convergence check
        create_input_polar(airfoil = airfoil,
                      Re = Re,
                      M = M, 
                      x_tr_top= x_tr_top,
                      x_tr_bot = x_tr_bot,
                      N_crit = N_crit, 
                      iterations = iterations + CONV_FLAG*100, #### Change 
                      alpha_min = 0, # should improve convergence
                      alpha_max = alpha_max - CONV_FLAG*2,
                      step = step)
        os.remove('data/'+airfoil + '.log')
        out_message = sp.run('Xfoil/bin/xfoil < input_file.in', shell=True, capture_output = True) # WANRNING: has to be adapted to the operating system
        print(f'Convergence failed {CONV_FLAG} times\n')
        CONV_FLAG += 1

        

def calculate_polar(airfoil,
                    Re = 350e3,
                    Ncrit = 9,
                    Xtr_top = 1, Xtr_bot = 1,
                    iterations = 500,
                    a_min = -5, 
                    a_max = 15, 
                    delta_a = 0.5):
    '''
    Run XFOIL with the input parameters in order to compute the polar and
    store the data in a numpy matrix.
    '''
    # write input file to run xfoil
    with open('input_file.in', 'w') as f:
        f.write(f'load {airfoil}\n')
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
    sp.run('Xfoil/bin/xfoil < input_file.in', shell = True)
    # save polar as np.array
    polar_data = np.loadtxt('polar.txt', skiprows=12)
    return polar_data

def calculate_BL(airfoil,
                 Re = 350000,
                 Ncrit = 9,
                 Xtr_top = 1, Xtr_bot = 1,
                 iterations = 100,
                 alpha = 5):
    '''
    Run XFOIL using the input parameters in order to compute all boundary 
    layer parameters and save them in a numpy matrix.
    '''
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
    sp.run('Xfoil/bin/xfoil < input_file.in', shell = True)
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
    '''
    Plot polar graphics for the data, passed as a Nx9 matrix (to be consistent with 
    the file format of the XFOIL output). 
    If parameters are passed in 'params', these can be displayed in the titles of the 
    figures. 
    If 'data_min' and 'data_max' are passed, the plots include error bars by shading the 
    area enclosed by these two curves.
    '''
    
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
    ax3.plot(CL, CL**1.5/CD)
    ax3.set_xlabel(r'$C_L$')
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
    


if __name__ == '__main__':
    '''
    The main here is only used for testing purpuses.
    '''
    directory = 'library'
    Merge(directory,[0.8234],'a0')
    Polar('a0',1)#,0.35e6,0,1,1,9,0,15,1) # airfoil, Re = 350e3, M = 0, x_tr_top = 1.0, x_tr_bot = 1.0, N_crit = 9.0, iterations = 100, alpha_min = -2.0, alpha_max = 15.0, step = 0.5

    
        
        

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
