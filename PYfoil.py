# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 16:23:10 2022

@author: Asus
"""

import subprocess as sp
import os
import numpy as np


def Merge(directory,
          percentage,
          newname,
          ):
    airfoils = os.listdir(directory)
    with open('input_file.in', 'w') as f:
        f.write('\nload library/' + airfoils[0] + '\n')
        for i in range(len(percentage)):
            f.write('inte\n')
            f.write('C\n')
            f.write('F\n')
            f.write('library/' + airfoils[i] + '\n') # 
            interpolating_fraction = percentage[i]
            f.write(str(interpolating_fraction) + '\n')
            f.write('generation/' + newname + '\n')# airfoil nam, to upload
            f.write('PANE\n')
            f.write('save\n')
            f.write('generation/' + newname + '.dat\n')
            f.write('quit')
    sp.run('Xfoil/bin/xfoil < input_file.in', shell=True)  # WANRNING: has to be adapted to the operating system

def Polar(airfoil,
          Re = 350e3,
          M = 0,
          x_tr_top = 1.0,
          x_tr_bot = 1.0,
          N_crit = 9.0,
          iterations = 100,
          alpha_min = 0.0,
          alpha_max = 10.0,
          step = 0.25
          ):
    
    with open('input_file.in', 'w') as f:
        #f.write(' \n')
        #f.write('plop\n')
        #f.write('g\n')
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
        # Deleate the polar if already existing
    if os.path.exists(airfoil + '.log'):
        os.remove(airfoil + '.log')
        # Run Xfoil   
    sp.run('Xfoil/bin/xfoil < input_file.in', shell=True) # WANRNING: has to be adapted to the operating system


if __name__ == '__main__':
    directory = r'C:\Users\Asus\Desktop\XFOIL\library' # directory's name where we store the airfoils
    Merge(directory,[0.5],'lol.dat')
    Polar('e220',0.35e6,0,1,1,9,100,0,15,1)

    
        
        

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    