import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import solve_ivp

def model(t, z, A, B, n):
    '''
    Diferential equations of plant-feedback model.
    '''
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

def model_timescale(t, z, A, B, n, epsilon):
    '''
    Diferential equations of plant-feedback model.
    '''
    #Separate plant and soil vecotrs and reshape
    p = np.array(z[0:n]).reshape(n, 1)
    q = np.array(z[n:n + n]).reshape(n, 1)
    #Create column vector of ones
    Ip = np.ones(shape = n).reshape(n, 1)
    Iq = np.ones(shape = n).reshape(n, 1)
    #Model equations
    dpdt = np.diag(p.transpose()[0]) @ (A @ q - \
            (p.transpose()[0] @ A @ q) * Iq) 
    dqdt = epsilon*np.diag(q.transpose()[0]) @ (B @ p - \
            (q.transpose()[0] @ B @ p) * Ip)
    return(list(dpdt.reshape(n)) + list(dqdt.reshape(n)))

def integrate_PSF(fun, t_span, z0, args, method = 'BDF', max_step = np.inf):
    '''
    Wrapper for integrator
    '''
    #Solve diferential equations
    sol = solve_ivp(fun, t_span, z0, method = method, args = args)

    #Get number of species 
    n = args[-1]
    #Get plant and soil abundances
    plant_ab = sol.y[0:n, :]
    soil_ab = sol.y[n:2*n, :]
    t = sol.t
    return (plant_ab, soil_ab, t)

def plot_solution(t, f_t, n_p):
    '''
    Plot solution
    '''
    #Get number of species 
    n = len(f_t)
    for i in range(n):
        if i < n_p:
            plt.plot(t, f_t[i,:])
        else:
            plt.plot(t, f_t[i,:], 'white')
    return 0

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

def check_feasibility(A, n):
    '''
    Check the feasibility of matrix A
    '''
    #Calculate equilibrium signs to check for feasibility
    eq_prop = np.linalg.inv(A) @ np.ones(n).reshape(n, 1)
    #Check if all have the same sign
    feasibility = np.all(eq_prop>0) if eq_prop[0] > 0 else np.all(eq_prop < 0)
    return feasibility 

def check_singularity(A):
    '''
    Check if matrix A is singular
    '''
    try:
        #Check for non-singularity
        Ainv = np.linalg.inv(A)
        singular = False
    except:
        #If error occurs when trying to invert, declare matrix as singular
        singular = True
    return singular

def remove_extinctions(plants, soils, tol):
    '''
    Remove those elements less than the tolerance in both vectors 
    simultaneously
    '''
    #Find indices of elements less than tolerance in each vector
    ext_ind_plants = np.where(plants < tol)[0]
    ext_ind_soils = np.where(soils < tol)[0]
    #Find index of elements less than tolerance in both vectors simultaneously
    ext_ind_both = np.intersect1d(ext_ind_plants, ext_ind_soils, 
                                  assume_unique = True)
    new_plants = np.delete(plants, ext_ind_both, axis = 0)
    new_soils = np.delete(soils, ext_ind_both, axis = 0)
    return(new_plants, new_soils)

def remove_extinctions_matrix(matrix, extinctions):
    '''
    Remove columns and rows corresponding to extinct plants/soils
    '''
    matrix_row = np.delete(matrix, extinctions, 0)
    matrix_new = np.delete(matrix_row, extinctions, 1)
    return matrix_new

def check_equilibrium(plant_vector, soil_vector):
    '''
    Check if the number of plants are soils are less or equal than 2
    '''
    n_plants = len(plant_vector)
    n_soils = len(soil_vector)
    if n_plants > 2 or n_soils > 2:
        #If either plant/soil have more than 2, we have not reached
        #equilibnrium
        equilibrium = False
    else:
        #Otherwise (neither have more than 2), declare equilibrium
        equilibrium = True
    return equilibrium

def check_equilibrium_bis(plant_vector, soil_vector):
    '''
    Check if only 2 species are greater than the tolerance
    '''
    n_plants = sum(plant_vector > 0)
    n_soils = sum(soil_vector > 0)
    if n_plants > 2:
        #If either plant/soil have more than 2, we have not reached
        #equilibnrium
        equilibrium = False
    else:
        #Otherwise (neither have more than 2), declare equilibrium
        equilibrium = True
    return equilibrium


def check_convergence(plants, soils, tol):
    '''
    Check that integration has converged by verifying that all abundances at
    equilibrium are less than 1
    '''
    conv_plants = np.all(plants <= 1+tol)
    conv_soils = np.all(soils <= 1+tol)
    if conv_plants and conv_soils:
        convergence = True
    else:
        convergence = False
    return convergence

def find_extinct_indices(plants, plants_rem, soils, soils_rem):
    '''
    Find extinct indices (those elements that are not in x, but not in x_rem)
    in both plants and soils.
    '''
    #Get elements in plants that are not in plants_rem
    ext_plant = list(set(plants) - set(plants_rem))
    #Find indices of extinct species in original vector
    o, ext_plant_ind, o = np.intersect1d(plants, ext_plant, 
                                         return_indices = True)
    ext_soil = list(set(soils) - set(soils_rem))
    o, ext_soil_ind, o = np.intersect1d(soils, ext_soil, 
                                        return_indices = True)

    return np.intersect1d(ext_plant_ind, ext_soil_ind)

def integrate_n(A, B, n, tol, n_int_max, epsilon):
    '''
    Integrate system as many times as needed until. Keep integrating system 
    until (1) it diverges, (2) it converges to an equilibrium with 2 or less 
    than 2 species, or (3) the number of integration cycles surpasses our 
    limit. 
    '''
    #Integrate
    plant_ab, soil_ab, t = integrate_PSF(model_timescale, [1, 2000], 2*n*[1/n], 
                                         (A, B, n, epsilon))
    #Check for convergence 
    convergence = check_convergence(plant_ab[:, -1], soil_ab[:, -1],
                                    tol)
    if not convergence:
        #If integration doesn't converge, we skip it.
        return 0
    #Set extinctions to 0
    plant_ab[plant_ab[:, -1] < tol, :] = 0
    #Check if we have reached equilibrium
    equilibrium = check_equilibrium_bis(plant_ab[:, -1], soil_ab[:, -1])
    #Initialize counter for number of integrations
    n_int = 0
    while (not equilibrium) & (n_int <= n_int_max):
        #Set initial conditions to the final state of previous 
        #integration
        z0 = list(np.hstack([plant_ab[:, -1], soil_ab[:, -1]]))
        #Re-integrate system
        plant_ab, soil_ab, t = integrate_PSF(model_timescale, [1, 2000], z0, 
                                             (A, B, n, epsilon))
        #Increase integration cycle counter
        n_int += 1
        #Check for convergence 
        convergence = check_convergence(plant_ab[:, -1], 
                                        soil_ab[:, -1], tol)
        if not convergence:
            return 0
        else:
            #Check equilibrium
            equilibrium = check_equilibrium_bis(plant_ab[:, -1], 
                                            soil_ab[:, -1])
    return(plant_ab, soil_ab, t)

def sample_A(n):
    '''
    Sample a feasible and non-singular A
    '''
    #Sample random matrix A 
    A_cand = np.random.random(size = (n, n))
    #Ensure that A is not singular and feasible
    singular = check_singularity(A_cand)
    feasible = check_feasibility(A_cand, n)
    while singular or not feasible:
        A_cand = np.random.random(size = (n, n))
        singular = check_singularity(A_cand)
        feasible = check_feasibility(A_cand, n)
    return A_cand
