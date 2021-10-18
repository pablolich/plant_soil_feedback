import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import solve_ivp

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
    ind_pairs = [(i+1, i+2) for i in range(len(sign_changes)-1)
                 if sign_changes[i] != sign_changes[i+1]]
    #Initialize vectors to store peaks and indices where peaks are reached
    n_peaks = len(ind_pairs)
    peaks = np.zeros(n_peaks)
    inds = np.zeros(n_peaks)
    for i in range(n_peaks):
        #Find values of time_series at these two peak candidate, and select the
        #maximum, because we are looking for peaks (not valleys)
        peaks[i]= max(time_series[ind_pairs[i][0]], time_series[ind_pairs[i][1]])
        #Store position where peak is reached (if there are various, take the
        #midle position)
        inds[i] = int(np.median(np.where(time_series == peaks[i])[0]))
    #Pack and return
    return zip(list(inds), list(peaks))

def cost(means):
    '''
    Calculates the cost of a certain clustering
    '''
    #Calculate total distance between all ordered pairs of group averages
    dist = np.mean(abs(means[:, np.newaxis] - means))/2
    #Number of groups
    n_groups = len(means)
    return n_groups/dist
    
def cluster(vector):
    '''
    Cluster 1-d vector by sorting it, and finding its largest gaps
    '''
    #Sort vector
    sorted_vec = np.sort(vector)
    #Get indicers that would sort the vector
    sorted_ind = np.argsort(vector)
    #Get where are the sorted elements in the original vector
    com, ind_original, ind2 = np.intersect1d(vector, sorted_vec, 
                                             return_indices = True)
    #Get differences between elements
    gaps = sorted_vec[1:] - sorted_vec[0:-1]
    #Count the number of 0 (that is, number of peaks with the same height)
    n_zeros = len(np.where(gaps == 0)[0])
    if n_zeros > 1:
        #Equilibrium has been reached
        return True
    #Initialize vector of cost of each partition (there are as many partitions 
    #as there are points in the time series)
    cost_vec = np.zeros(len(sorted_vec))
    #While switch and counter...
    go = True
    count = 0
    #Set minimum number of groups for the first partition
    n_groups = 2
    while go:
        #Store previous candidate partition
        try: prev_candidate = candidate_partition
        #If this is the first iteration, do nothing
        except: pass
        #Find the largest n_groups-1 gaps (note that one gap yields two groups)
        max_gaps = np.flip(np.sort(gaps))[0:n_groups-1]
        #Find indices of these gaps
        gap, gap_ind, max_gap_ind = np.intersect1d(gaps, max_gaps, 
                                                   return_indices = True)
        #Get indices of these gaps in the original vector
        gap_ind_origin = sorted_ind[gap_ind]
        #Re-group the points in n_groups groups, divided by gap_ind
        candidate_partition = np.split(sorted_vec, np.sort(gap_ind+1))
        #Calculate mean of each group
        mean_groups = [np.mean(i) for i in candidate_partition]
        #Obtain cost of this partition
        cost_vec[count] = cost(np.array(mean_groups))
        #If this is the first iteration, do nothing, but increase counters
        if count == 0: 
            count += 1
            n_groups += 1
        #Check if the cost is getting better or worse in the next iteration 
        #(we are okay with getting stuck in a local minima)
        else: 
            if cost_vec[count - 1] < cost_vec[count]:
                #If previous partition was better, terminate the loop
                go = False
            else:
                #If previous partition was worse, partition further 
                count += 1
                n_groups += 1
    #Get number of elements in each group
    n_each = [len(i) for i in prev_candidate]
    #Create vector of membersips in the sorted vector
    memberships_sorted = np.repeat(np.arange(len(n_each)), n_each)
    #Assign memberships to each element of the original vector
    memberships = memberships_sorted[ind_original]
    return np.array([[vector],[memberships]])

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
        print('Ones are reached')
        return True
    #Check if last element of all rows is (that is, all elements of the vectors
    #last_count_plant and last_count_soil are greater than 1)
    elif (np.all(last_count_plant > 1)) & (np.all(last_count_soil > 1)):
        print('Constant solution is reached')
        return True
    #Check for periodicity
    #First, check for 2.1
    elif (np.all(np.sum((diff_plant < tol_float), axis = 1) > 2)) & \
         (np.all(np.sum((diff_soil < tol_float), axis = 1) > 2)): 
        print('Checking for periodicity')
        ##If 2.1 is satisfied, find the peak values and position of the 
        #periodic solution 
        bool_peak_heights = np.zeros(2*n_r)
        #Put plants and soils together
        sol = np.vstack([plant_ab, soil_ab])
        for i in range(2*n_r):
            peaks = find_peaks(sol[i,:])
            #Cluster them
            peak_vals = np.array(list(peaks))[:,1]
            peak_groups  = cluster(peak_vals)
            if peak_groups == True:
                return True
            #Check if the heights are bounded, increasing, or decreasing
            else:
                import ipdb; ipdb.set_trace(context = 20)
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

