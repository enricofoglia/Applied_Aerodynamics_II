# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:38:57 2022

@author: Alberto Zigante
"""


# Object function
def obj_fun(nn, hypP, constraints, Cl_minD = 1.1, Cl_maxD = 1.7, Cl_cruise = 1.33): #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To call the script as a function
        
    import numpy as np 
    

        
    #OPEN FILES & LOADING POLARS
    #Initialization of the parameters used by the objective function
    Cl_max          = np.zeros(nn)
    End_Cl_133      = np.zeros(nn)
    Delta_alpha     = np.zeros(nn)
    End_max         = np.zeros(nn)
    End_integral    = np.zeros(nn)
    
    for airfoil in range(nn):
    
        name = 'data/a'+str(airfoil)+'.log'#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% .txt for XFLR5 outputs
                                            # % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% .log for XFOIL outputs
    
        
        DATAlog = np.loadtxt(name, skiprows=12) #Acquisition from xfoil output %%%%%%%%%%%% skiprows=11 for XFLR5 outputs
                                                                                # % %%%%%%% skiprows=12 for XFOIL output
        alpha   = DATAlog[:,0]
        Cl      = DATAlog[:,1]
        Cd      = DATAlog[:,2]
        
        #Stall condition: max Cl and relative angle of incidence
        Cl_max[airfoil] = Cl[len(Cl)-1]
        Alpha_Clmax = alpha[len(alpha)-1]
        for i in range(len(Cl)-1):
            if Cl[i] > Cl_max[airfoil]:
                Cl_max[airfoil] = Cl[i]
                Alpha_Clmax = alpha[i]
        
        #Design condition: Cl=1.33, here looking at the FIRST value greater than 1.33
        for i in range(len(Cl)-1):
            if (Cl[i]-Cl_cruise)>0:
                Alpha_Cl133 = alpha[i]
                End_Cl_133[airfoil] = Cl[i]**1.5/Cd[i]
                Delta_alpha[airfoil] = Alpha_Clmax-Alpha_Cl133
                break
        
        #Computation of the Endurance only in the range of Cl that is of our interest
        #Outside [Cl_minD; Cl_maxD] the endurance is not computed : zeros
        End = np.zeros(len(Cl))
        for i in range(len(Cl)-1):
            #if Cl[i] > Cl_minD or Cl[i] < Cl_maxD:
                End[i] = Cl[i]**1.5/Cd[i]
        
        End_max[airfoil]= max(End)
        End[Cl<Cl_minD] = 0
        End[Cl>Cl_maxD] = 0
        End_integral[airfoil] = sum(End[i]*(Cl[i+1]-Cl[i]) for i in range(len(Cl)-1))
    #-----------------------------------------------------------------------------
    
    
    
    #COLLECTING ALL THE DATA FOR THE OBJECTVIE FUNCTION
    #          |  End_max | Cl_max  |  Delta_alpha  |  End_Cl:1.33  |   End_integral
    #------------------------------------------------------------------------------
    #airfoil1  |          |         |               |               |
    #airfoil2  |          |         |               |               |
    #   ...    |          |         |               |               |
    
    data = np.zeros((nn, 5))
    data[:,0] = End_max
    data[:,1] = Cl_max
    data[:,2] = Delta_alpha
    data[:,3] = End_Cl_133
    data[:,4] = End_integral
        
    #worst = [min(row) for row in zip(*data)] #Array with the lowest values per each column, to normalize all the others
    #-----------------------------------------------------------------------------
    
    
    
    # CONSTRAINT
    S_cl, S_alpha = constraints #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%input of the function
    
    for i in range(nn) :
        if data[i,1] < S_cl*Cl_maxD or data[i,2] < S_alpha:
            data[i,:] = np.zeros(5)                      # The function of the airfoil will be zero
    #-----------------------------------------------------------------------------
    
    
    
    # TODO: NORMALIZATION
    
    #-----------------------------------------------------------------------------
    
    
    
    # OBJECTIVE FUNCTION
    f       = np.zeros(nn)                          #Array with the values of the objective function
    
    #For each polar all the parameters are multipled for their hyperparameters and are summed
    for i in range(nn) :
        f[i] = sum(hypP[j]*data[i,j] for j in range(5))
    #-----------------------------------------------------------------------------
    

    return f #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The function returns the objective function













