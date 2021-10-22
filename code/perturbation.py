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
    while not like : 
        #Sample matrices A and B with H-S constraints
        A, B = sampling_matrices(n)
        #Integrate system
        plant_ab, soil_ab, t = integrate_PSF(model, [1, 20000], 2*n*[1/n], 
                                             (A, B, n))
        #Check convergence
        conv = check_convergence(plant_ab, soil_ab)
        #Perturb system
        p = 0.01
        A_per = A + np.random.uniform(1-p, 1+p, size = np.shape(A))
        plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [1, 20000], 
                                                         2*n*[1/n], 
                                                         (A_per, B, n))
        #Check convergence
        conv_per = check_convergence(plant_ab_per, soil_ab_per)
        while not conv or not conv_per:
            A, B = sampling_matrices(n)
            A_per = A + np.random.uniform(1-p, 1+p, size = np.shape(A))
            #Integrate unperurbed system
            plant_ab, soil_ab, t = integrate_PSF(model, [1, 20000], 
                                                 2*n*[1/n], 
                                                 (A, B, n))
            #Check convergence of unperurbed system
            conv = check_convergence(plant_ab, soil_ab)
            #Integrate perurbed system
            plant_ab_per, soil_ab_per, t_per = integrate_PSF(model, [1, 20000], 
                                                             2*n*[1/n], 
                                                             (A_per, B, n))
            #Check convergence of perurbed system
            conv_per = check_convergence(plant_ab_per, soil_ab_per)
        #Plot 
        ax1 = plt.subplot(121)
        plot_solution(t, plant_ab)
        ax2 = plt.subplot(122)
        plot_solution(t_per, plant_ab_per)
        plt.show()
        n_it += 1
        if n_it > n_it_max:
            break
    #Ask input of record or not
    #Record A, plant_ab, soil_ab for unperturbed and perturbed if yes

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

