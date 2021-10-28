#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
import progressbar
from functions import *
import pandas as pd

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Number of species
    n = 4
    #Set extinct threshold
    tol = 1e-9
    like = 0
    n_int_max = 20
    #Perturbation strenght
    p = 0.1
    #Set to False convergence flags
    conv = False
    conv_per = False
    #Set number of simulations
    n_sim = 2
    sim_it = 0
    #Initialize emtpy lists containin variables of interest.
    sim_vec = []
    t_vec = []
    t_per_vec = []
    p_l = []
    p_ab = []
    p_ab_per = []
    s_l = []
    s_ab = []
    s_ab_per = []
    cyc_vec = []
    while sim_it < n_sim:
        print('\nSimulation number: ', sim_it)
        #Keep sampling and integrating until convergence
        #Initialize list to keep track of extinct plants and soils
        ext_ind = []
        while not conv or not conv_per:
            #Sample matrices A and B with H-S constraints
            A, B = sampling_matrices(n)
            A_per = A + np.random.uniform(1-p, 1+p, size = np.shape(A))
            B_per = B
            #Integrate unperurbed system
            plant_ab, soil_ab, t = integrate_PSF(model, [0, 2000], 2*n*[1/n], 
                                                 (A, B, n))
            #Check convergence of unperurbed system
            conv = check_convergence(plant_ab[:, -1], soil_ab[:, -1], tol)
            #Integrate perurbed system
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [0, 2000], 
                                                             2*n*[1/n], 
                                                             (A_per, B_per, n))
            #Check convergence of perurbed system
            conv_per = check_convergence(plant_ab_per[:, -1], 
                                         soil_ab_per[:, -1], tol)
        #Remove extinct species from perturbed system
        rem_plant, rem_soil = remove_extinctions(plant_ab_per[:, -1], 
                                                 soil_ab_per[:, -1], 
                                                 tol) 
        #Check if the convergent perturbed system has arrived to equilibrium
        equilibrium = check_equilibrium(rem_plant, 
                                        rem_soil)
        #Update number of plants
        n_plant = len(rem_plant)
        #Remove extinct species from payoff matrices
        ext_plant, ext_soil = find_extinct_indices(plant_ab_per[:, -1], 
                                                   rem_plant, 
                                                   soil_ab_per[:, -1],
                                                   rem_soil)
        #Record extinctions
        ext_ind += list(ext_plant)
        #Set number of integration cycles to zero
        n_int = 0
        #Store
        t_vec += list(np.tile(t_per, n))
        t_points = len(t_per)
        cyc_vec += t_points*n*[n_int]
        if ext_ind:
            plant_ab_store = np.insert(plant_ab, 0, ext_ind, axis = 0)
            soil_ab_store = np.insert(soil_ab, 0, ext_ind, axis = 0)
        #Set labels of plants and soils in long format by repeating each 
        #label t_points times
        p_l += list(np.repeat(np.arange(n), t_points))
        s_l += list(np.repeat(np.arange(n), t_points))
        #Set simulation vector by repeating each simulation the number of 
        #time points times the number of speciess)
        sim_vec += t_points*n*[sim_it]
        #Need to check that the three abovee vectors remain the same length
        p_ab += list(plant_ab.flatten())
        s_ab += list(soil_ab.flatten())
        p_ab_per += list(plant_ab_per.flatten())
        s_ab_per += list(soil_ab_per.flatten())
        #Keep integrating both systems until equilibrium in perturbed one 
        #is reached
        while not equilibrium:
            A_per = remove_extinctions_matrix(A_per, ext_plant)
            B_per = remove_extinctions_matrix(B_per, ext_soil)
            #Record extinctions
            ext_ind += list(ext_plant)
            #Set initial conditions to the final state of previous integration
            z0_per = list(np.hstack([rem_plant, rem_soil]))
            #Re-integrate systems
            import ipdb; ipdb.set_trace(context = 20)
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [0, 2000], 
                                                             z0_per, 
                                                             (A_per, B_per, 
                                                              n_plant))
            #Remove extinctions after integration
            rem_plant, rem_soil = remove_extinctions(plant_ab_per[:, -1], 
                                                     soil_ab_per[:, -1], tol) 
            #Set initial conditions of unperturbed system for reintegration
            z0 = list(np.hstack([plant_ab[:, -1], soil_ab[:, -1]]))
            #Set time span based on the time vector of the perturbed 
            #integration
            #Integrate unperturb system for the same time as the perturb system 
            print(len(z0))
            print(A.shape)
            print(B.shape)
            import ipdb; ipdb.set_trace(context = 20)
            sol = solve_ivp(model, [0, 2000], z0, args = (A, B, n), 
                            dense_output = True)
            #Check for convergence in both cases
            conv_per = check_convergence(plant_ab_per[:, -1], 
                                         soil_ab_per[:, -1], tol)
            conv = check_convergence(plant_ab[:, -1], soil_ab[:, -1], tol)
            if not conv_per or not conv:
                n_int += 1
                #Get bad indices
                ind_rem = np.where(np.array(sim_vec) == sim_it)[0]
                if sim_it > 0:
                    import ipdb; ipdb.set_trace(context = 20)
                    #Remove information of divergent system 
                    del p_ab[ind_rem]
                    del s_ab[ind_rem]
                    del t_vec[ind_rem]
                    del p_l[ind_rem]
                    del s_l[ind_rem]
                    del sim_vec[ind_rem]
                    del p_ab_per[ind_rem]
                    del s_ab_per[ind_rem]
                    del cyc_vec[ind_rem]
                #Don't count this simulation
                sim_it -= 1
                print('\nSimulation does not converge')
                #End current simulation
                break
            #Increase integration cycle counter
            n_int += 1
            #Evaluate unperturbed solution in time vector of perturbed one.
            z = sol.sol(t_per)
            plant_ab = z[0:n, :]
            soil_ab = z[n:2*n, :]
            #If it is the case, insert rows of 0 at extinct positions
            if ext_ind:
                plant_ab_store = np.insert(plant_ab, 0, ext_ind, axis = 0)
                soil_ab_store = np.insert(soil_ab, 0, ext_ind, axis = 0)
            #Store variables of interest
            #Note to add the end point of previous integration
            t_vec += list(np.tile(t_per+t_vec[-1], n))
            t_points = len(t_per)
            cyc_vec += t_points*n*[n_int]
            #Set labels of plants and soils in long format by repeating each 
            #label t_points times
            p_l += list(np.repeat(np.arange(n), t_points))
            s_l += list(np.repeat(np.arange(n), t_points))
            #Set simulation vector by repeating each simulation the number of 
            #time points times the number of speciess)
            sim_vec += t_points*n*[sim_it]
            #Need to check that the three abovee vectors remain the same length
            p_ab += list(plant_ab_store.flatten())
            s_ab += list(soil_ab_store.flatten())
            p_ab_per += list(plant_ab_per.flatten())
            s_ab_per += list(soil_ab_per.flatten())
            #Check equilibrium
            equilibrium = check_equilibrium(rem_plant, rem_soil)
            print('Number of integrations: ', n_int)
            if n_int > n_int_max:
                print('\nExceeded number of integration cycles')
                break

        #Create dataframe with lists
        #Store number of integration cycles at the end

        #Store dynamics of each integration in a matrix

        #Store matrix A and A_per if the system is convergent throughout
        #all the process of relaxation to a 2sp. equilibrium
        if sim_it > n_sim:
            print('\nCompleted number of simulations')
            break
        sim_it += 1
    #Set data frame column names vector
    names = ['n_sim', 't', 'p_l', 'p_ab', 'p_ab_per', 's_l', 's_ab', 
             's_ab_per', 'n_cyc']
    df = pd.DataFrame({'n_sim': sim_vec,
                       't': t_vec, 
                       'p_l': p_l, 
                       'p_ab': p_ab,
                       'p_ab_per': p_ab_per,
                       's_l': s_l,
                       's_ab': s_ab,
                       's_ab_per': s_ab_per,
                       'cyc_vec': cyc_vec
                      })
    import ipdb; ipdb.set_trace(context = 20)
    df.to_csv('../data/example_search.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

