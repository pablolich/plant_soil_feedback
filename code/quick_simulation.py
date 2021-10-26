#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from functions import *
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''

    A = np.array(pd.read_csv('../data/A_pre.csv', delimiter = ' ', 
                 header = None))
    A_per = np.array(pd.read_csv('../data/A_post.csv', delimiter = ' ', 
                     header = None))
    mu = pd.read_csv('../data/mu.csv', delimiter = ' ', header = None)
    B = np.diag(np.array(mu).flatten())
    B_per = B
    n = len(A)
    n_per = n
    z0 = 2*n*[1/n] 
    z0_per = z0
    tol = 1e-9
    equilibrium = 0
    n_int = 0
    ext_ind = []
    p_ab_per = []
    t_vec = []
    p_l = []
    present = np.ones(n)
    t_max = 500
    while not equilibrium:
        #Integrate unperurbed system
        plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [0, t_max], 
                                                         z0_per, 
                                                         (A_per, B, n),
                                                         method = 'RK45', 
                                                         max_step = 1)
        #Delete extinctions 
        rem_plant, rem_soil = remove_extinctions(plant_ab_per[:, -1], 
                                                 soil_ab_per[:, -1], 
                                                 tol) 
        #Check if we have reached equilibrium
        equilibrium = check_equilibrium_bis(plant_ab_per[:, -1], 
                                            soil_ab_per[:, -1])
        #Update number of plants
        n_per = len(rem_plant)
        #Find extinct indices of plants and soils 
        ext_plant, ext_soil = find_extinct_indices(plant_ab_per[:, -1], 
                                                   rem_plant, 
                                                   soil_ab_per[:, -1],
                                                   rem_soil)
        #Set extinct species to 0
        plant_ab_per[plant_ab_per[:, -1] < tol, :] = 0
        #soil_ab_per[soil_ab_per[:, -1] < tol, :] = 0
        #Store plant abundance
        p_ab_per += list(plant_ab_per.flatten())
        #Store time
        try:
            t_vec += list(np.tile(t_per+t_vec[-1], n))
        except:
            t_vec += list(np.tile(t_per, n))
        #Store species labels
        p_l += list(np.repeat(np.arange(n), len(t_per)))
        #Set initial conditions to the final state of previous integration 
        z0_per = list(np.hstack([plant_ab_per[:, -1], soil_ab_per[:, -1]]))
        n_int += 1
    sol = solve_ivp(model, [0, t_max*n_int], z0, args = (A, B, n), 
                    dense_output = True, method = 'RK45')
    #Store solution
    df = pd.DataFrame({'p_ab_per':p_ab_per,
                       't':t_vec,
                       'p_l':p_l
                       }) 
    #Sort it according to label
    df = df.sort_values(by = ['p_l', 't'])
    #Get a continuous time vector
    ind = np.where(np.array(p_l) == 0)[0]
    t_cont = np.array(t_vec)[ind] 
    #Evaluate unperturbed solution in time vector of perturbed one.
    z = sol.sol(t_cont)
    plant_ab = z[0:n, :]
    #Add to df
    df['p_ab'] = plant_ab.flatten()
    df.to_csv('../data/zach_results.csv')


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

