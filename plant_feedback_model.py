#!/usr/bin/env python3

__appname__ = '[model.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pylab as plt
import progressbar

## CONSTANTS ##


## FUNCTIONS ##

def model(t, z, params):

    '''
    Diferential equations of plant-feedback model.
    '''
    #Unpack parameter values
    A, B, n, n = map(params.get, ('A', 'B', 'n', 'n'))
    #Separate plant and soil vecotrs and reshape
    p = np.array(z[0:n]).reshape(n, 1)
    q = np.array(z[n:2*n]).reshape(n, 1)
    #Create column vector of ones
    Ip = np.ones(n).reshape(n, 1)
    Iq = np.ones(n).reshape(n, 1)
    #Model equations
    dpdt = np.diag(p.transpose()[0]) @ (A @ q - \
            (p.transpose()[0] @ A @ q) * Iq) 
    dqdt = np.diag(q.transpose()[0]) @ (B @ p - \
            (q.transpose()[0] @ B @ p) * Ip)
    return(list(dpdt.reshape(n)) + list(dqdt.reshape(n)))

def sampling_matrices(n):
    '''
    Sampling parameters to construct A and B supporting full coexistence
    '''
    #Random (positive) diagonal matrix
    B = np.diag(np.random.rand(n)) 
    #Random negative constant
    c = -np.random.rand(1)    
    #Random vector by which to shift cols of A (+1 to ensure A > 0)
    col_shifts = np.ones(n) + np.random.random(n) 
    #First create matrix of constant columns
    A = np.ones(shape = (n, n)) @ np.diag(col_shifts)
    #Make A a rescaled version of B
    A = A + c * B 
    return(A, B)


def main(argv):
    '''Main function'''
    
    #Number of plants (and soils)
    n = 10
    n_sim = 100
    #Extinct threshold
    tol = 1e-3
    #Preallocate
    n_plants = np.zeros(n_sim)
    n_soils = np.zeros(n_sim)

    for i in progressbar.progressbar(range(n_sim)):
        #Set matrices
        #A, B = sampling_matrices(n, n)
        A = np.random.random(size = (10, 10))
        B = np.diag(np.random.random(size = 10))
        params = {'A':A,
                  'B':B, 
                  'n':n,
                  'n':n
                 }
        #Create time vector
        tspan = tuple([1, 2000])
        #Set initial conditions
        z0 = list(np.ones(2*n)/n)
        #Solve diferential equations
        sol = solve_ivp(lambda t,z: model(t,z, params),
                        tspan, z0,
                        method = 'BDF', atol = 0.0001 )
        #Check for convergence
        if not all(sol.y[:,-1] <= 1):
            #Skip this iteration
            continue
        #Get plant and soil abundances
        plant_ab = sol.y[0:n, -1]
        soil_ab = sol.y[n:2*n, -1]
        #Number of plant and soil survivors
        n_plants[i] = len(plant_ab[(plant_ab > tol) & (plant_ab <= 1)])
        n_soils[i] = len(soil_ab[(soil_ab > tol) & (soil_ab <= 1)])
        if (n_plants[i] > 2) | (n_soils[i] > 2):
            #Plot 
            for i in range(np.shape(sol.y)[0]):
                linestyle = 'solid'
                color = 'green'
                if i > n:
                    linestyle = 'dashed'
                    color = 'black'
                plt.plot(sol.t, sol.y[i,:], linestyle = linestyle, 
                         color = color)
            plt.show()
    plt.hist(n_plants, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    plt.hist(n_soils, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
    
