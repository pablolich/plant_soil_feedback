#!/usr/bin/env python3

__appname__ = '[model.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
import progressbar
from functions import *
import pandas as pd

def main(argv):
    '''Main function'''
    #Number of plants (and soils)
    n_vec = np.array([2, 3, 5, 6])
    #Names of columns
    names = list(['n_sim', 'n_p', 'n_p_f'])
    #Number of simulations
    n_sim = 5000
    #Set maximum number of integrations before interruption and skip to next 
    #simulation
    n_int_max = 100
    #Preallocate dataframe
    df = pd.DataFrame(0, index = np.arange(len(n_vec)*n_sim), columns = names)
    #Fill known data
    df.loc[:, 'n_p'] = np.repeat(n_vec, n_sim)
    df.loc[:, 'n_sim'] = np.tile(np.arange(n_sim), len(n_vec))
    #Extinct threshold
    tol = 1e-9
    #Iterator for each initial species number (goes from 0 to 3)
    n_vec_it = 0
    #Iterate over different values of inital species 
    for n in n_vec:
        #Current simulation number
        n_act = 0
        #Run n_sim simulations for n starting plant and soil species
        while n_act < n_sim:
            #Sample random matrix A 
            A_cand = np.random.random(size = (n, n))
            #Ensure that A is not singular and that it is feasible
            singular = check_singularity(A_cand)
            while singular:
                A_cand = np.random.random(size = (n, n))
                singular = check_singularity(A_cand)
            #After ensuring non-singularity of A_cand, declare it as A
            A = A_cand
            #Sample matrix B
            B = np.diag(np.random.random(size = n))
            #Integrate system
            plant_ab, soil_ab, t = integrate_PSF(model, [1, 2000], 2*n*[1/n], 
                                                 (A, B, n))
            convergence = check_convergence(plant_ab[:, -1], soil_ab[:, -1],
                                            tol)
            #Check for convergence 
            if not convergence:
                #If integration doesn't converge, we skip it.
                continue
            #Delete extinctions (only those that happen simultaneously in both
            #vectors). Failing to do so would be a conceptual mistake, since
            #If a plant goes extinct, so does its soil.
            rem_plant, rem_soil = remove_extinctions(plant_ab[:, -1], 
                                                     soil_ab[:, -1], 
                                                     tol) 
            #Check if we have reached equilibrium
            equilibrium = check_equilibrium(rem_plant, rem_soil)
            #Update number of plants
            n_plant = len(rem_plant)
            ###################################################################
            #NOTE THAT THIS IS NOT ENTIRELY CORRECT, NEED TO ADD A SNIPPET
            #CHECKING FOR PARTNERSHIP/ZERO-SUM GAME
            ###################################################################
            #Initialize counter for number of integrations (integration cycles 
            #will stop if n_int is more than n_int_max)
            n_int = 0
            #Keep integrating system until (1) it diverges, (2) it 
            #converges to an equilibrium with 2 or less than 2 species, or 
            #(3) the number of integration cycles surpasses our limit. 
            while (not equilibrium) & (n_int <= n_int_max):
                #Find extinct indices of plants and soils 
                ext_plant, ext_soil = find_extinct_indices(plant_ab[:, -1], 
                                                           rem_plant, 
                                                           soil_ab[:, -1],
                                                           rem_soil)
                A = remove_extinctions_matrix(A, ext_plant)
                B = remove_extinctions_matrix(B, ext_soil)
                #Set initial conditions to the final state of previous 
                #integration
                z0 = list(np.hstack([rem_plant, rem_soil]))
                #Solve diferential equations again
                #Re-integrate system
                plant_ab, soil_ab, t = integrate_PSF(model, [1, 2000], z0, 
                                                     (A, B, n_plant))
                #Increase integration cycle counter
                n_int += 1
                #Check for convergence 
                convergence = check_convergence(plant_ab[:, -1], 
                                                soil_ab[:, -1], tol)
                if not convergence:
                    #Don't count this simulation
                    n_act -= 1
                    #End current simulation
                    break
                else:
                    #Remove extinctions after integration
                    rem_plant, rem_soil = remove_extinctions(plant_ab[:, -1], 
                                                             soil_ab[:, -1], 
                                                             tol) 
                    n_plant = len(rem_plant)
                    #Check equilibrium
                    equilibrium = check_equilibrium(rem_plant, rem_soil)
            #Check if while loop exited due to exceeding integration cycles or
            #if the integration diverged
            if not convergence:
                #If it did, flag as non-convergent integration
                df.loc[n_act + n_sim*n_vec_it, 'n_p_f'] = -1
                #Don't count this simulation
                n_act -= 1
            elif n_int > n_int_max:
                #couln't trim the system down to two species
                df.loc[n_act + n_sim*n_vec_it, 'n_p_f'] = -2
                #Don't count this simulation
                n_act -= 1
            else:
                df.loc[n_act + n_sim*n_vec_it, 'n_p_f'] = len(rem_plant)
            n_act += 1
            print(n_act, 'equilibria reached for ', n_vec[n_vec_it], 
                  'starting species', end = '\r')
        n_vec_it += 1
    df.to_csv('../data/results.csv')
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
