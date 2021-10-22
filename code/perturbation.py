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
    n = 6
    like = 0
    n_it_max = 20
    n_it = 0
    #Perturbation strenght
    p = 0.01
    #Set to False convergence flags
    conv = False
    conv_per = False
    #Set number of simulations
    n_sim = 2
    sim_it = 0
    while sim_it < n_sim:
        #Keep sampling and integrating until convergence
        while not conv or not conv_per:
            #Sample matrices A and B with H-S constraints
            A, B = sampling_matrices(n)
            A_per = A + np.random.uniform(1-p, 1+p, size = np.shape(A))
            #Integrate unperurbed system
            plant_ab, soil_ab, t = integrate_PSF(model, [1, 2000], 
                                                 2*n*[1/n], 
                                                 (A, B, n))
            #Check convergence of unperurbed system
            conv = check_convergence(plant_ab, soil_ab)
            #Integrate perurbed system
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [1, 2000], 
                                                             2*n*[1/n], 
                                                             (A_per, B, n))
            #Check convergence of perurbed system
            conv_per = check_convergence(plant_ab_per, soil_ab_per)
        #Check if the convergent perturbed system has arrived to equilibrium
        equilibrium = check_equilibrium(plan_ab_per[:, -1], 
                                        soil_ab_per[:, -1])
        #Keep integrating both systems until equilibrium is reached
        while not equilibrium:
            #Organize new initial conditions and re-integrate, evaluating
            #equilibrium at the end of each simulation

            #Store number of integration cycles at the end

            #Store dynamics of each integration in a matrix

            #Store matrix A and A_pert if the system is convergent throughout
            #all the process of relaxation to a 2sp. equilibrium

        sim_it += 1
        if sim_it > n_sim:
            break
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

