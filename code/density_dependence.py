#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from functions import *
import pandas as pd

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #set tolerance
    tol=1e-3
    #number of species 
    n=3
    #run n_sim simulations for n starting plant and soil species
    running=True
    #find first community
    pos = 0
    c = 0.1
    while running:
        #sample random matrix A 
        A = sample_A(n)
        #sample matrix B
        B = np.diag(np.random.random(size = n))
        #solve system numerically
        sol = solve_ivp(model_self_regulation, t_span = [1, 2000], y0 = 2*n*[1/n], 
                        method = 'BDF', args =(A, B, n, 0))
        sol_dens = solve_ivp(model_self_regulation, t_span = [1, 2000], 
                             y0 = 2*n*[1/n], method = 'BDF',
                             args =(A, B, n, c)) 
        #check for divergence
        if np.any(sol.y[:, -1] > 1+tol) or np.any(sol_dens.y[:, -1] > 1+tol):
            continue
        end_ab = sol.y[:n, -1]
        n_sp = len(np.where(end_ab > tol)[0])
        ax1 = plt.subplot(321)
        plot_solution(sol.t, sol.y, n)
        ax2 = plt.subplot(322)
        plot_solution(sol_dens.t, sol_dens.y, n)
        #next goal, find community with n_rem species
        n_rem = search_next(n_sp)
        #search for dynamics with n_rem ending species
        while n_sp != n_rem: 
            #sample random matrix A 
            A = sample_A(n)
            #sample matrix B
            B = np.diag(np.random.random(size = n))
            #solve system numerically
            sol = solve_ivp(model_self_regulation, t_span = [1, 2000], y0 = 2*n*[1/n], 
                            method = 'BDF', args=(A, B, n, 0))
            sol_dens = solve_ivp(model_self_regulation, t_span = [1, 2000], 
                                 y0 = 2*n*[1/n], method = 'BDF',
                                 args=(A, B, n, c)) 
            #check for divergence
            end_ab = sol.y[:n, -1]
            n_sp = len(np.where(end_ab > tol)[0])
            if np.any(sol.y[:, -1] > 1+tol) or np.any(sol_dens.y[:, -1] > 1+tol):
                continue
        ax3= plt.subplot(323)
        plot_solution(sol.t, sol.y, n)
        ax4 = plt.subplot(324)
        plot_solution(sol_dens.t, sol_dens.y, n)
        #generate dynamics wiht special parametrization 
        A, B = sampling_matrices(n)
        #solve system numerically
        sol = solve_ivp(model_self_regulation, t_span = [1, 2000], y0 = 2*n*[1/n], 
                        method = 'BDF', args =(A, B, n, 0))
        sol_dens = solve_ivp(model_self_regulation, t_span = [1, 2000], 
                             y0 = 2*n*[1/n], method = 'BDF',
                             args =(A, B, n, c)) 
        running=False
        if np.any(sol.y[:, -1] > 1+tol) or np.any(sol_dens.y[:, -1] > 1+tol):
            running = True
        ax5 = plt.subplot(325)
        plot_solution(sol.t, sol.y, n)
        ax6 = plt.subplot(326)
        plot_solution(sol_dens.t, sol_dens.y, n)
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
