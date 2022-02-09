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
    #vector of timescales
    epsilon = np.array([1, 10, 100])
    #Number of species 
    n=3
    #Extinct threshold
    tol = 1e-9
    #Set maximum number of integrations before skip to next simulation
    n_int_max = 100
    #set number of final species
    n_final = 0
    #Search for parameters 3->2 and 3->1
    #Run n_sim simulations for n starting plant and soil species
    while n_final != 2:
        #Sample random matrix A 
        A_cand = np.random.random(size = (n, n))
        #Ensure that A is not singular and feasible
        singular = check_singularity(A_cand)
        feasible = check_feasibility(A_cand, n)
        while singular or not feasible:
            A_cand = np.random.random(size = (n, n))
            singular = check_singularity(A_cand)
            feasible = check_feasibility(A_cand, n)
        #Sample random matrix A 
        A = A_cand
        #Sample matrix B
        B = np.diag(np.random.random(size = n))
        n_final = integrate_n(A, B, n, tol, n_int_max)
    #Run dynamics with different timescales
    for e in epsilon:
        p, q, t = integrate_PSF(model_timescale, [1, 2000], 2*n*[1/n], 
                                (A, B, n, e)) 
        plot_solution(t, np.vstack([p, q]), n_final)
        plt.show()





    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

