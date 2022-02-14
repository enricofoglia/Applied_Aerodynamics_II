# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:38:57 2022

@author: Alberto Zigante
"""
import numpy as np 


# Object function
def obj_fun(name, hypP, constraints, Cl_minD = 1.1, Cl_maxD = 1.7, Cl_cruise = 1.33):
        
    
       
    #OPEN FILE & ACQUISITION
    #Initialization of the parameters used by the objective function
    Cl_max          = 0
    End_Cl_133      = 0
    Delta_alpha     = 0
    End_max         = 0
    End_integral    = 0
    
    #Open file
    name = 'data/'+name+'.log'#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% .txt for XFLR5 outputs
                                    #% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% .log for XFOIL outputs

    acquisition = np.loadtxt(name, skiprows=12) #Acquisition from xfoil output %%%%%%%%%%%% skiprows=11 for XFLR5 outputs
                                                                            #% %%%%%%%%%%%% skiprows=12 for XFOIL output
    alpha   = acquisition[:,0]
    Cl      = acquisition[:,1]
    Cd      = acquisition[:,2]
    
    #STALL CONDITION: max Cl and relative angle of incidence
    Cl_max = Cl[len(Cl)-1]
    Alpha_Clmax = alpha[len(alpha)-1]
    for i in range(len(Cl)-1):
        if Cl[i] > Cl_max :
            Cl_max      = Cl[i]
            Alpha_Clmax = alpha[i]

    #DESIGN CONDITION: Cl=1.33, here looking at the FIRST value greater than 1.33
    for i in range(len(Cl)-1):

        if (Cl[i]-Cl_cruise)>0:
            Alpha_Cl133 = alpha[i]
            End_Cl_133  = Cl[i]**1.5/Cd[i]
            break
    Delta_alpha = Alpha_Clmax-Alpha_Cl133
            
    #ENDURANCE: computation only in the range of Cl that is of our interest
    #Outside [Cl_minD; Cl_maxD] the endurance is not computed : zeros
    End = np.zeros(len(Cl))
    for i in range(len(Cl)-1):
            End[i] = Cl[i]**1.5/Cd[i]
    
    End_max = max(End)
    End[Cl<Cl_minD] = 0
    End[Cl>Cl_maxD] = 0
    End_integral = sum(End[i]*(Cl[i+1]-Cl[i]) for i in range(len(Cl)-1))

    data = [End_max-100, Cl_max-1, Delta_alpha, End_Cl_133-100, End_integral-50]
    #-------------------------------------------------------------------------
    
    # CONSTRAINT
    S_cl, S_alpha = constraints
    
    if Cl_max < S_cl*Cl_maxD or Delta_alpha < S_alpha:
        data = np.zeros(5)                      # The function of the airfoil will be zero
    #-----------------------------------------------------------------------------
     
    
    # OBJECTIVE FUNCTION
    #All the parameters are multiplied for their hyperparameters and are summed
    f = sum(hypP[i]*data[i] for i in range(5))   
    #-----------------------------------------------------------------------------
    if f<=0: # avoid negative values 
        f = 1

    return f

if __name__ == '__main__':
    hypP    = np.multiply([1, 1, 1, 1, 1], 10)      #Hyperparameters, to be set manually (trial and error)
    constraints = (1.2, 6.0)  
    f = obj_fun('e_dae', hypP, constraints)












