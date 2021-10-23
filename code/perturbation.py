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
    n_it_max = 20
    #Perturbation strenght
    p = 0.01
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
    while sim_it < n_sim:
        #Keep sampling and integrating until convergence
        while not conv or not conv_per:
            #Sample matrices A and B with H-S constraints
            A, B = sampling_matrices(n)
            A_per = A + np.random.uniform(1-p, 1+p, size = np.shape(A))
            #Integrate unperurbed system
            plant_ab, soil_ab, t = integrate_PSF(model, [1, 2000], 2*n*[1/n], 
                                                 (A, B, n))
            #Check convergence of unperurbed system
            conv = check_convergence(plant_ab, soil_ab)
            #Integrate perurbed system
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [1, 2000], 
                                                             2*n*[1/n], 
                                                             (A_per, B, n))
            #Check convergence of perurbed system
            conv_per = check_convergence(plant_ab_per, soil_ab_per)
        #Remove extinct species from perturbed system
        rem_plant, rem_soil = remove_extinctions(plant_ab_per[:, -1], 
                                                 soil_ab_per[:, -1], 
                                                 tol) 
        #Check if the convergent perturbed system has arrived to equilibrium
        equilibrium = check_equilibrium(rem_plant, 
                                        rem_soil)
        #Set number of integration cycles to zero
        n_int = 0
        #Keep integrating both systems until equilibrium is reached
        while not equilibrium:
            #Remove extinct species from payoff matrices
            ext_plant, ext_soil = find_extinct_indices(plant_ab_per[:, -1], 
                                                       rem_plant, 
                                                       soil_ab_per[:, -1],
                                                       rem_soil)
            A_per = remove_extinctions_matrix(A_per, ext_plant)
            B = remove_extinctions_matrix(B, ext_soil)
            #Set initial conditions to the final state of previous integration
            z0_per = list(np.hstack([rem_plant, rem_soil]))
            #Re-integrate systems
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [1, 2000], 
                                                             z0_per, 
                                                             (A_per, B, 
                                                              n_plant))
            #Remove extinctions after integration
            rem_plant, rem_soil = remove_extinctions(plant_ab_per[:, -1], 
                                                     soil_ab_per[:, -1], tol) 
            #Set initial conditions of unperturbed system for reintegration
            z_0 = list(np.hstack([plant_ab[:, -1], soil_ab[:, -1]]))
            #Set time span based on the time vector of the perturbed 
            #integration
            t_span = [0, max(t_per) - min(t_per)]
            #Integrate unperturb system for the same time as the perturb system 
            plant_ab, soil_ab, t = integrate_PSF(model, 
                                                 t_span,
                                                 z_0, (A, B, n))
            #Check for convergence in both cases
            conv_per = check_convergence(plant_ab_per[:, -1], 
                                            soil_ab_per[:, -1], tol)
            conv = check_convergence(plant_ab[:, -1], soil_ab[:, -1], tol)
            if not conv_per or not conv:
                #Don't count this simulation
                sim_it -= 1
                #End current simulation
                break
            #Increase integration cycle counter
            n_int += 1
            #Store variables of interest
            t_vec = t_vec.append(list(t_per))
            t_points = len(t_per)
            #Set labels in long format by repeating each label t_points times
            p_l = p_l.append(np.repeat(np.arange(n), t_points))
            s_l = s_l.append(np.repeat(np.arange(n), t_points))
            #Set simulation vector by repeating each simulation the number of 
            #time points times the number of speciess)
            sim_vec = sim_vec.append(t_points*n*[sim_it])
            #Need to check that the three abovee vectors remain the same length
            #Need to create a vector of abundances including the 0s of extinct
            #species
            p_ab = []
            s_ab = []
            p_ab_per = []
            s_ab_per = []
            #Check equilibrium
            equilibrium = check_equilibrium(rem_plant, rem_soil)

        #Store number of integration cycles at the end

        #Store dynamics of each integration in a matrix

        #Store matrix A and A_per if the system is convergent throughout
        #all the process of relaxation to a 2sp. equilibrium

        sim_it += 1
        if sim_it > n_sim:
            break
        #Set data frame column names vector
        names = ['n_sim', 't', 'p_l', 'p_ab', 's_l', 's_ab']
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

