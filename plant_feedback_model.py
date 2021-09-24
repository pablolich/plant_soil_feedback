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
    A, B, n_p, n_m = map(params.get, ('A', 'B', 'n_p', 'n_m'))
    #Separate plant and soil vecotrs and reshape
    p = np.array(z[0:n_p]).reshape(n_p, 1)
    q = np.array(z[n_p:n_p+n_m]).reshape(n_m, 1)
    #Create column vector of ones
    Ip = np.ones(n_p).reshape(n_p, 1)
    Iq = np.ones(n_m).reshape(n_m, 1)
    #Model equations
    dpdt = np.diag(p.transpose()[0]) @ (A @ q - \
            (p.transpose()[0] @ A @ q) * Iq) 
    dqdt = np.diag(q.transpose()[0]) @ (B @ p - \
            (q.transpose()[0] @ B @ p) * Ip)
    return(list(dpdt.reshape(n_p)) + list(dqdt.reshape(n_m)))

def sampling_matrices(n_p, n_m):
    '''
    Sampling parameters to construct A and B supporting full coexistence
    '''
    #Sample d_j's
    d_vec = np.random.rand(n_p)
    #Sample c_j's
    c_vec = np.random.rand(n_p)
    #Sample a constant c (less than zero)
    c = -1*abs(np.random.rand(1))
    #Construct off-diagonal elements of matrix of the rescaled game C
    C_vec =  -d_vec/c
    C = np.diag(C_vec) @ np.ones(shape = (n_p, n_p))
    #Sample diagonals of C 
    c_diag = np.random.rand(n_p)
    np.fill_diagonal(C, c_diag)
    #Create matrices
    A = C + c_vec.reshape(n_p, 1)
    B = c*C.T + d_vec
    return(A, B)


def main(argv):
    '''Main function'''
    
    #Number of plants
    n_p = 10
    #Number of soil communities
    n_m = 10
    n_sim = 10000
    #Extinct threshold
    tol = 1e-3
    #Preallocate
    n_plants = np.zeros(n_sim)
    n_soils = np.zeros(n_sim)

    for i in progressbar.progressbar(range(n_sim)):
        #Set matrices
        A, B = sampling_matrices(n_p, n_m)
        #A = np.random.randint(10, size = (10, 10))
        #B = np.diag(np.random.randint(10, size = 10))
        params = {'A':A,
                  'B':B, 
                  'n_p':n_p,
                  'n_m':n_m
                 }
        #Create time vector
        tspan = tuple([1, 200])
        #Set initial conditions
        z0 = list(0.1*np.ones(n_p)) + list(0.1*np.ones(n_m))
        #Solve diferential equations
        sol = solve_ivp(lambda t,z: model(t,z, params),
                        tspan, z0,
                        method = 'BDF', atol = 0.0001 )
        #Check for convergence
        if not all(sol.y[:,-1] <= 1):
            #Skip this iteration
            continue
        #Get plant and soil abundances
        plant_ab = sol.y[0:n_p, -1]
        soil_ab = sol.y[n_p:n_p + n_m, -1]
        #Number of plant and soil survivors
        n_plants[i] = len(plant_ab[(plant_ab > tol) & (plant_ab <= 1)])
        n_soils[i] = len(soil_ab[(soil_ab > tol) & (soil_ab <= 1)])
        if (n_plants[i] < 10) | (n_soils[i] < 10):
            #Plot 
            for i in range(np.shape(sol.y)[0]):
                linestyle = 'solid'
                color = 'green'
                if i > n_p:
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
    
