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
    q = np.array(z[n:n + n]).reshape(n, 1)
    #Create column vector of ones
    Ip = np.ones(shape = n).reshape(n, 1)
    Iq = np.ones(shape = n).reshape(n, 1)
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

def find_sign_change(array):
    '''
    Detect sign changes in array
    Output:
        array: where 0 means there is no sign change with respect the previous
               element, and 1 means that there is.
    '''
    asign = np.sign(array)
    #Consider 0 as 'positive'
    asign[np.where(asign == 0)[0]] = 1
    return ((np.roll(asign, 1) - asign) != 0).astype(int)[1:]

def find_peaks(time_series):
    '''
    Find the peaks of a periodic time series

    Parameters: 
        time_series (1xt): Vector of time series data.

    Output:
        peaks (dict): Dictionary with indices and values of peaks
    '''
    #Calculate vector of diferences
    diff_vec = time_series[1:] - time_series[0:-1]
    #Get a vector indicating where are sign changes 
    sign_changes = find_sign_change(diff_vec)
    #Get indices before and after the sign change 
    import ipdb; ipdb.set_trace(context = 20)
    ind_pairs = [(i+1, i+2) for i in sign_changes 
                 if sign_changes[i] != sign_changes[i+1]]
    #Initialize vectors to store peaks and indices where peaks are reached
    n_peaks = len(ind_pairs)
    peaks = np.zeros(n_peaks)
    inds = peaks
    for i in range(n_peaks):
        #Find values of time_series at these two peak candidate, and select the
        #maximum, because we are looking for peaks (not valleys)
        peaks[i]= max(time_series[ind_pairs[i][0]], time_series[ind_pairs[i][1]])
        #Store position where peak is reached (if there are various, take the
        #midle position)
        inds[i] = int(np.median(np.where(time_series == peaks[i])[0]))
    #Pack and return
    return zip(list(inds), list(peaks))
    
vector  = np.array([1, 2, 3, 3, 2, 1])
peaks = find_peaks(vector)

def check_equilibrium(plant_ab, soil_ab, tol, n, equilibrium, tol_float):
    '''
    Check wether the obtain equilibrium is a real one by checking if either:
        1. The abundances have remained constant for a while
        2. The abundances are periodic. To check if the solution is periodic
           we: 
           2.1. Check that there are values of solution repeat
           2.2. Check that value of solution at the peaks remains constant
           2.3. Check that time between peaks remains constant

    Parameters: 
        plant_ab (txn): Plant abundance matrix across integrated t_span
        soil_ab (txn): Soil abundance matrix across integrated t_span
        tol (float): Tolerance under which we consider a species extinct
        n (int): Number of plants/soils
        tol_float (float): Tolerance under which we consider two peaks the same

    Returns:
        equil (bool): Boolean stating whether the equilibrium was real (True) or
                      not (False).
    '''
    #Check that after extinctions removal we have equal number of 
    #plants and soils
    if len(plant_ab) != len(soil_ab):
        return False
    #Get number of remaining species
    n_r = len(plant_ab)
    #Check if the equilibrium vector remains constant for a while.
    #Get counts of last element of abundance vector of each species 
    last_count_plant = np.zeros(n_r)
    last_count_soil = last_count_plant 
    for i in range(n_r):
        u, plant_counts = np.unique(np.round(plant_ab[i, :], 
                                             decimals = int(-np.log10(tol))), 
                                    return_counts = True)
        last_count_plant[i] = plant_counts[-1] 
        u, soil_counts = np.unique(np.round(soil_ab[i, :], 
                                            decimals = int(-np.log10(tol))), 
                                    return_counts = True)
        last_count_soil[i] = soil_counts[-1]
    #Get last point of soil and plant abundances
    endpoint_plant = plant_ab[:, -1]
    endpoint_soil = soil_ab[:, -1]
    #Get difference between time-series and endpoints
    diff_plant = abs(endpoint_plant.reshape(n_r, 1) - plant_ab)
    diff_soil = abs(endpoint_soil.reshape(n_r, 1) - soil_ab) 
    #Check if there is any fixations
    if (np.any(plant_ab[:, -1] == 1)) | (np.any(soil_ab[:, -1] == 1)):
        return True
    #Check if last element of all rows is (that is, all elements of the vectors
    #last_count_plant and last_count_soil are greater than 1)
    elif (np.all(last_count_plant > 1)) & (np.all(last_count_soil > 1)):
        return True
    #Check for periodicity
    #First, check for 2.1
    elif (np.all(np.sum((diff_plant < tol_float), axis = 1) > 2)) & \
         (np.all(np.sum((diff_soil < tol_float), axis = 1) > 2)): 
        ##If 2.1 is satisfied, check find the peaks of the periodic solution
        #for i in range(n_r):
        #    #For the plants
        #    peaks_plants = find_peaks(
        return True
    #There are neither static equilibria nor periodic oscilations
    else:
        return False 

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
    nim = 100
    n_act = 0
    #Extinct threshold
    tol = 1e-9
    #Preallocate
    nlants = np.zeros(nim)
    noils = np.zeros(nim)
    n_int_max = 10

    while n_act <= nim:
        #Sample matrices A and B with H-S constraints
        #A, B = sampling_matrices(n)
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
    plt.hist(nlants, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    plt.hist(noils, bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
    
