#!/usr/bin/env python3

__appname__ = '[model.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
import progressbar
from functions import *

def main(argv):
    '''Main function'''
    
    #Number of plants (and soils)
    n = 5
    #Number of simulations
    n_sim = 10
    #Current simulation number
    n_act = 0
    #Extinct threshold
    tol = 1e-9
    #Preallocate number of coexisting plants/soils
    n_plants = np.zeros(n_sim)
    n_soils = np.zeros(n_sim)
    #Set maximum number of integrations before declaring not equilibrium
    n_int_max = 10
    #Start integration
    while n_act <= n_sim:
        #Sample matrices A and B with H-S constraints
        #A, B = sampling_matrices(n)
        #Sample A without H-S constraints
        A = np.random.random(size = (n, n))
        #Ensure that A is not singular and that it is feasible
        sing_unfeas = True
        while sing_unfeas:
            #Keep sampling A matrices until feasible and non-singular
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
                #A, B = sampling_matrices(n)
                A = np.random.random(size = (n, n))
        #Sample matrix B
        B = np.diag(np.random.random(size = n))
        #Create dictionary of parameters
        params = {'A':A,
                  'B':B, 
                  'n':n}
        #Create time vector
        tspan = tuple([1, 2000])
        #Set initial conditions
        z0 = list(np.ones(2*n)/n)
        #Solve diferential equations
        sol = solve_ivp(lambda t,z: model(t,z, params),
                        tspan, z0,
                        method = 'BDF', atol = tol)
        #Get plant and soil abundances
        plant_ab = sol.y[0:n, :]
        soil_ab = sol.y[n:2*n, :]
        #Check for convergence and also check wether all extinct soils have 
        #their corresponding plant extinct too
        correct = plant_soil_extinction(plant_ab, soil_ab, tol)
        if  (not all(sol.y[:,-1] <= 1)) | (not correct):
            #Skip this iteration
            continue
        #Find rows of extinct species
        rows_rem_plant = np.where(plant_ab[:, -1] < tol)[0]
        rows_rem_soil = np.where(soil_ab[:, -1] < tol)[0]
        #Get rid of rows with species at equilibrium with abundance below the 
        #tolerance threshold
        plant_ab_rem = np.delete(plant_ab, rows_rem_plant, axis = 0) 
        soil_ab_rem = np.delete(soil_ab, rows_rem_soil, axis = 0)
        #Plot either way
        for i in range(np.shape(sol.y)[0]):
            linestyle = 'solid'
            color = 'green'
            if i >= n:
                linestyle = 'dashed'
                color = 'black'
            plt.plot(sol.t, sol.y[i,:], linestyle = linestyle, 
                     color = color)
        plt.show()

        #Check if true equilibrium is reached
        equilibrium_val = check_equilibrium(plant_ab_rem, soil_ab_rem, 
                                            tol, n, equilibrium = False, 
                                            tol_float = 1e-3)
        print('Test of equilibrium says that it is: ', equilibrium_val)
        n_int = 0
        while (not equilibrium_val) & (n_int <= n_int_max):
            #Keep integrating system until (1) it diverges, or (2) it converges
            #to a valid equilibrium as determined with the check_equilibrium 
            #function

            #Plot either way
            for i in range(np.shape(sol.y)[0]):
                linestyle = 'solid'
                color = 'green'
                if i >= n:
                    linestyle = 'dashed'
                    color = 'black'
                plt.plot(sol.t, sol.y[i,:], linestyle = linestyle, 
                         color = color)
            plt.show()

            #If not, update parameters for reintegration
            n_r = len(plant_ab_rem)
            n = len(soil_ab_rem)
            if n_r != n:
                #The number of soils is not the same as the number of plants, 
                #so we keep integrating until 1 or 2 (see above) are satisfied
                #Get initial conditions of endpoints at previous integration
                z0 = list(np.hstack([plant_ab[:,-1], soil_ab[:, -1]]))
                #Solve diferential equations again
                sol = solve_ivp(lambda t,z: model(t,z, params),
                                tspan, z0,
                                method = 'BDF', atol = 0.0001 )
                #Get plant and soil abundances
                plant_ab = sol.y[0:n, :]
                soil_ab = sol.y[n:2*n, :]
                #Check for convergence 
                if not all(sol.y[:,-1] <= 1):
                    #Terminate while loop trying to find equilibruim for this 
                    #system
                    break
                else:
                    #Check if true equilibrium is reached again
                    equilibrium_val = check_equilibrium(sol.y[0:n,:], 
                                                        sol.y[n:2*n,:], 
                                                        tol, n, 
                                                        equilibrium = False, 
                                                        tol_float = 1e-3)
                    print('Keep integrating...: ', n_int, end = '\r')
                    n_int += 1
            else:
                #The number of soils is the same as the number of plants, but 
                #equilibrium has not yet been reached. Then we keep integrating
                #but we update the system by removing the extinct species.
                #Update number of species and soils
                params['n'] = n_r
                params['n'] = n_r
                params['n'] = n
                A_r = np.delete(A, rows_rem_plant, 0)
                A_r = np.delete(A_r, rows_rem_plant, 1)
                B_r = np.delete(B, rows_rem_soil, 0)
                B_r = np.delete(B_r, rows_rem_soil, 1)
                #Update matrices A and B
                params['A'] = A_r
                params['B'] = B_r
                #Set initial conditions to the final state of previous 
                #integration
                z0 = list(np.hstack([plant_ab_rem[:,-1], soil_ab_rem[:, -1]]))
                #Solve diferential equations again
                sol = solve_ivp(lambda t,z: model(t,z, params),
                                tspan, z0,
                                method = 'BDF', atol = 0.0001 )
                #Get plant and soil abundances
                plant_ab = sol.y[0:n, :]
                soil_ab = sol.y[n:2*n, :]
                #Check for convergence and also check wether all extinct soils
                #have their corresponding plant extinct too
                #Check for convergence 
                if not all(sol.y[:,-1] <= 1):
                    #Terminate while loop trying to find equilibruim for this 
                    #system
                    break
                else:
                    #Check if true equilibrium is reached again
                    equilibrium_val = check_equilibrium(sol.y[0:n_r,:], 
                                                        sol.y[n_r:2*n_r,:], 
                                                        tol, n_r,
                                                        equilibrium = False, 
                                                        tol_float = 1e-3)
                    print('Keep integrating...: ', n_int, end = '\r')
                    n_int += 1
        for i in range(np.shape(sol.y)[0]):
            linestyle = 'solid'
            color = 'green'
            if i >= params['n']:
                linestyle = 'dashed'
                color = 'black'
            plt.plot(sol.t, sol.y[i,:], linestyle = linestyle, 
                     color = color)
        plt.show()
        n_act += 1
        print('Number of equilibria reached: ', n_act, end = '\r')
    plt.hist(n_plants, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    plt.hist(n_soils, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
    
