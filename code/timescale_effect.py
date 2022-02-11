#!/usr/bin/env python3

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from functions import *

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #vector of timescales
    epsilon = np.array([1, 5, 10])
    #set tolerance
    tol=1e-9
    #Number of species 
    n=3
    #Run n_sim simulations for n starting plant and soil species
    divergent=True
    while divergent:
        #Sample random matrix A 
        A = sample_A(n)
        #Sample matrix B
        B = np.diag(np.random.random(size = n))
        #Loop over timescales
        sim=0
        fig, axs = plt.subplots(3)
        fig.subplots_adjust(hspace=0.5)
        for i in epsilon:
            sol = solve_ivp(model_timescale, t_span = [1, 500], y0 = 2*n*[1/n], 
                            method = 'BDF', args =(A, B, n, i))
            #Check for divergence
            if np.any(sol.y[:, -1] > 1+tol):
                plt.close()
                break
            for j in range(len(sol.y)):
                if j < n:
                    if j == n-2:
                        axs[sim].plot(sol.t, sol.y[j,:], 'g', label='plant')
                    axs[sim].plot(sol.t, sol.y[j,:], 'g')
                    axs[sim].set_title('$\epsilon = %i$' %i)
                else:
                    if j == n:
                        axs[sim].plot(sol.t, sol.y[j,:], 'red', alpha=0.5, 
                                      label = 'soil', linestyle = 'dashed')
                    else: 
                        axs[sim].plot(sol.t, sol.y[j,:], 'red', alpha=0.5,
                                      linestyle = 'dashed') 
            sim += 1
        if sim == len(epsilon):
            divergent = False
    plt.legend()
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):

    status = main(sys.argv)
    sys.exit(status)
