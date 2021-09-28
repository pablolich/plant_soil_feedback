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
    A, B, n = map(params.get, ('A', 'B', 'n'))
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

def check_equilibrium(plant_ab, soil_ab, tol):
    '''
    Check wether the obtain equilibrium is a real one by checking if either:
        1. The abundances have remained constant for a while
        2. The abundances are periodic.

    Parameters: 
        plant_ab (txn): Plant abundance matrix across integrated t_span
        soil_ab (txn): Soil abundance matrix across integrated t_span

    Returns:
        equil (bool): Boolean stating whether the equilibrium was real (True) or
                      not (False).
    '''

    #Get rid of species below the tolerance threshold
    plant_ab = plant_ab[plant_ab[:,-1] > tol, :]
    #Check if the equilibrium vector remains constant for a while.
    equilibrium = False

    return equilibrium

def plant_soil_extinction(plant_ab, soil_ab, tol):
    '''
    Determine if extinct soils also have their corresponding plant extinct
    '''
    #Get indices of extinct soils
    ext_soil_ind = set(np.where(soil_ab < tol)[0])
    #Get indices of extinct plants
    ext_plants_ind = set(np.where(plant_ab < tol)[0])
    #Make sure that plants are contained in soil (that is, all numbers in 
    #ext_soil_ind are also in ext_plants_ind
    contained = ext_plants_ind.issubset(ext_soil_ind)
    return contained



def main(argv):
    '''Main function'''
    
    #Number of plants (and soils)
    n = 5
    n_sim = 100
    #Extinct threshold
    tol = 1e-3
    #Preallocate
    n_plants = np.zeros(n_sim)
    n_soils = np.zeros(n_sim)

    for i in progressbar.progressbar(range(n_sim)):
        #Set matrices
        A, B = sampling_matrices(n)
        A = np.random.random(size = (n, n))
        #Ensure that A is not singular and that it is feasible
        sing_unfeas = True
        while sing_unfeas:
            try:
                #Check for non-singularity
                Ainv = np.linalg.inv(A)
                #Calculate equilibrium signs to check for feasibility
                eq_prop = np.linalg.inv(A) @ np.ones(n).reshape(n, 1)
                if np.all(eq_prop > 0):
                    #If both conditions hold use this A matrix
                    sing_unfeas = False
                else:
                    #Otherwise, raise an exception and sample another A matrix
                    raise Exception('Matrix A has unfeasible equilibrium')
            except:
                A = np.random.random(size = (n, n))
        B = np.diag(np.random.random(size = n))
        params = {'A':A,
                  'B':B, 
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
        #Get plant and soil abundances
        plant_ab = sol.y[0:n, -1]
        soil_ab = sol.y[n:2*n, -1]
        #Check for convergence and also check wether all extinct soils have 
        #also their corresponding plant extinct too
        correct = plant_soil_extinction(plant_ab, soil_ab, tol)
        if  (not all(sol.y[:,-1] <= 1)) | (not correct):
            #Skip this iteration
            continue
        #Check if this is true equilibrium
        check_equilibrium(sol.y[0:n,:], sol[n:2*n,:])
        if (n_plants[i] > 2) | (n_soils[i] > 2):
            #In the case that we have more than 2 plants or soils in the final
            #community, re-integrate for another 2000 timesteps
            #Set initial conditions
            z0 = list(plant_ab) + list(soil_ab)
            #Solve diferential equations
            sol = solve_ivp(lambda t,z: model(t,z, params),
                            tspan, z0,
                            method = 'BDF', atol = 0.0001 )
            for i in range(np.shape(sol.y)[0]):
                linestyle = 'solid'
                color = 'green'
                if i > n:
                    linestyle = 'dashed'
                    color = 'black'
                plt.plot(sol.t, sol.y[i,:], linestyle = linestyle, 
                         color = color)
            plt.show()
            import ipdb; ipdb.set_trace(context = 20)
            #Check for convergence
            if not all(sol.y[:,-1] <= 1):
                #Skip this iteration
                continue
        else:
            #Number of plant and soil survivors
            n_plants[i] = len(plant_ab[(plant_ab > tol) & (plant_ab <= 1)])
            n_soils[i] = len(soil_ab[(soil_ab > tol) & (soil_ab <= 1)])
    plt.hist(n_plants, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    plt.hist(n_soils, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
    
