import numpy as np
import matplotlib

def extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, xo = 500.0):

    if theta < -10.0:
        theta = -10.0

    Ao = np.power(xo, 2.0)
    
    if ld_list['area'] >= Ao:
  
        elevation_f = np.array(ld_list['elevation'] - base_elevation)
        elevation = np.append(elevation, np.array(elevation_f))
        index = ld_list['index']
        chi_f = chi_o + (1 / np.array(ld_list['area']))**theta[0] * np.array(ld_list['distance_scale']) * de[index[0], index[1]]
        chi = chi + [chi_f]
        if ld_list.get('next') is not None:
            for next_list in ld_list['next']:
                elevation, chi = extract_chi_elevation_values(next_list, de, theta, chi_f, elevation, chi, base_elevation, xo = xo)
    
    return elevation, chi    

def extract_dA_elevation_values(ld_list):
  
    all_elevation = np.array(ld_list['elevation'])
    all_dA = np.array(ld_list['dA'])
    
    if ld_list.get('next') is not None:
        for next_list in ld_list['next']:
            elevation, dA = extract_dA_elevation_values(next_list)
            all_elevation = np.append(all_elevation,elevation)
            all_dA = np.append(all_dA,dA)
    
    return all_elevation, all_dA  

def extract_profile_values(ld_list, xo = 500.0, items = ()):
    
    Ao = np.power(xo, 2.0)
    return_list = list()
    
    if ld_list.get('area') is None or ld_list.get('area') >= Ao:
        this_list = list()
        for arg in items:
            this_list += [ld_list[arg]]
        return_list.append(this_list)
        if ld_list.get('next', None) is not None:
            for next_list in ld_list['next']:
                return_list += extract_profile_values(next_list, xo = xo, items = items)
    
    return return_list

def chi_elevation(ld_list, de, theta, xo = 500.0):
    
    chi_o = 0
    elevation = []
    chi = []
    base_elevation = ld_list['elevation']
    e, c = extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, xo=xo)
    return np.array(e), np.array(c)

def best_ks_with_wrss_list(ld_list, de, theta, xo = 500):
    
    e, c = chi_elevation(ld_list, de, theta, xo=xo)
    A = np.vstack([c]).T
    sol = np.linalg.lstsq(A, e)
    m = sol[0]
    WRSS = sol[1]    
    
    return (m, WRSS)

def best_ks_with_r2_list(ld_list, de, theta, xo = 500):
    (ks, WRSS) = best_ks_with_wrss_list(ld_list, de, theta, xo =xo)
    SS = uninformative_SS_list(ld_list, de, xo)
    return (ks, 1 - (WRSS / SS))
    
def uninformative_SS_list(ld_list, de, xo = 500):
    
    e, c = chi_elevation(ld_list, de, [0.5], xo=xo)
    mean_elevation = np.mean(e)
    return np.sum(np.power(e-mean_elevation, 2))

def best_ks_and_theta_with_wrss_list(ld_list, de, xo = 500, maxiter = 100, maxfun = 200):
    
    import scipy.optimize
    chi_ks = lambda theta: best_ks_with_wrss_list(ld_list, de, theta, xo=xo)[1]
    if len(chi_ks([0.5])) == 0:
        return (0, 0, 0)
    (xopt, funval, iter, funcalls, warnflag) = scipy.optimize.fmin(chi_ks, np.array([0.5]), (), 1E-5, 1E-5, 100, 200, True, True, 0, None)
    (m, WRSS) = best_ks_with_wrss_list(ld_list, de, [xopt[0]], xo)
    SS = uninformative_SS_list(ld_list, de, xo)
    
    R2 = 1 - (WRSS / SS)
    if warnflag == 1 or warnflag == 2:
        R2 = 0.0
    return (m, xopt[0], R2)
    
def best_ks_and_theta_with_wrss(elevation, flow_direction_or_length, area, outlet, xo = 500):
    
    ld_list = flow_direction_or_length.map_values_to_recursive_list(outlet, area = area, elevation = elevation)
    de = area._mean_pixel_dimension()
    
    return best_ks_and_theta_with_wrss_list(ld_list, de, xo = xo)

def hi(elevation, flow_direction, dA, outlet):
    ld_list = flow_direction.map_values_to_recursive_list(outlet, dA = dA, elevation = elevation)
    return hi_list(ld_list)

def hi_list(ld_list):
    
    elevation, dA = extract_dA_elevation_values(ld_list)
   
    if np.nan in elevation:
        return 0.0

    mean_elevation = np.mean(elevation)
    max_elevation = np.max(elevation)
    min_elevation = np.min(elevation) 

    return (mean_elevation - min_elevation) / (max_elevation - min_elevation) 

def area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, minimum_area = 1.0E7, min_pathlength_in_cells=100):
    
    import dem as d
    mean_pixel_dimension = d.BaseSpatialGrid()
    mean_pixel_dimension._copy_info_from_grid(area, True)
    mean_pixel_dimension._griddata = area._mean_pixel_dimension()
    
    ld_list = flow_direction.map_values_to_recursive_list(outlet, elevation = elevation, area = area, de = mean_pixel_dimension)
    
    area = [ld_list['area']]
    elevation = [ld_list['elevation']]
    de = [ld_list['de']*ld_list['distance_scale']]
    tributary_ld = []
    #first we gather information of area, elevation, and de from our flow direction, which should map out our mainstem and tributary
    def get_elevations_and_areas(ld_list, area, elevation, de, tributary_ld, minimum_area_to_consider):
        maximum_area = 0.0
        if ld_list.get('next', None) is not None:
            for next in ld_list['next']:
                if (next['area'] > minimum_area_to_consider) and (next['area'] > maximum_area):
                    maximum_area = next['area'] #for loop that iteratively adds to the area if the next area from ld list is greater than the previous
            
            for next in ld_list['next']:
                if next['area'] == maximum_area:
                    area += [next['area']]
                    elevation += [next['elevation']]
                    de += [next['de'] * next['distance_scale']]
                    (area, elevation, de, tributary_ld) = get_elevations_and_areas(next, area, elevation, de, tributary_ld, minimum_area_to_consider)
                elif next['area'] > minimum_area_to_consider:
                    tributary_ld.append(next)   # if your area is equal to the max area, then you must be the mainstem, and add up all those values. else if, add those values to a tributary list 
        return (area, elevation, de, tributary_ld)
    
    (area, elevation, de, tributary_ld) = get_elevations_and_areas(ld_list, area, elevation, de, tributary_ld, minimum_area)
    # get elevations and areas gets the elevation and area for the mainstem
    area = [area]
    elevation = [elevation]
    de = [de]
    
    while len(tributary_ld) > 0:
        next_tributary_ld = []
        for trb_ld in tributary_ld:
            this_area = [trb_ld['area']]
            this_elevation = [trb_ld['elevation']]
            this_de = [trb_ld['de']+trb_ld['distance_scale']]
            (this_area, this_elevation, this_de, next_tributary_ld) = get_elevations_and_areas(trb_ld, this_area, this_elevation, this_de, next_tributary_ld, minimum_area)
            if(len(this_area) >= min_pathlength_in_cells):
                area.append(this_area)
                elevation.append(this_elevation)
                de.append(this_de)
        tributary_ld = next_tributary_ld

    return_area = []
    return_elevation = []
    return_de = []
    
    for (this_area, this_elevation, this_de) in zip(area, elevation, de):
        this_return_elevation = []
        this_return_de = []
        this_return_area = []
        for (elevation_value, de_value, area_value) in zip(this_elevation, this_de, this_area):
            if elevation_value is not None and de_value is not None and area_value is not None:
                this_return_elevation += [elevation_value - this_elevation[0]]
                this_return_de += [de_value]
                this_return_area += [area_value]
        
        if len(this_return_elevation) > 4:
            return_elevation.append(this_return_elevation)
            return_de.append(this_return_de)
            return_area.append(this_return_area)

    return (return_area, return_elevation, return_de)

def indexes_for_mainstem_and_tributaries(outlet, flow_direction, area, minimum_area = 1.0E7):
        
    ld_list = flow_direction.map_values_to_recursive_list(outlet, area = area)
    
    area = [ld_list['area']]
    indexes = [ld_list['index']]
    tributary_ld = []
    
    def get_indexes(ld_list, indexes, area, tributary_ld, minimum_area_to_consider):
        maximum_area = 0.0
        for next_ld in ld_list['next']:
            if (next_ld['area'] > minimum_area_to_consider) and (next_ld['area'] > maximum_area):
                maximum_area = next_ld['area']
        
        for next_ld in ld_list['next']:
            if next_ld['area'] == maximum_area:
                indexes += [next_ld['index']]
                (indexes, tributary_ld) = get_indexes(next_ld, indexes, area, tributary_ld, minimum_area_to_consider)
            elif next_ld['area'] > minimum_area_to_consider:
                tributary_ld.append(next_ld)   
        return (indexes, tributary_ld)
    
    (indexes, tributary_ld) = get_indexes(ld_list, indexes, area, tributary_ld, minimum_area)
    
    indexes = [indexes]
    
    while len(tributary_ld) > 0:
        next_tributary_ld = []
        for trb_ld in tributary_ld:
            this_index = [trb_ld['index']]
            this_area = [trb_ld['area']]
            (this_index, next_tributary_ld) = get_indexes(trb_ld, this_index, this_area, next_tributary_ld, minimum_area)
            indexes.append(this_index)
        tributary_ld = next_tributary_ld

    return_index = []
    
    for index in indexes:
        this_return_index = []
        for index_value in index:
            if index_value is not None:
                this_return_index += [index_value]
        
        if len(this_return_index) > 4:
            return_index.append(this_return_index)

    return return_index

    
def best_ks_theta_wrss_for_outlet(outlet, flow_direction, elevation, area, minimum_area = 1E7, min_pathlength_in_cells = 100):

    def chi_for_profile(area, de, theta):
        chi = []
        chi_value = 0.0
        for (a, d) in zip(area, de):
            chi += [chi_value]
            chi_value += (1/a)**theta[0] * d
        return chi
    
    def best_ks_with_wrss(chi, elevation):
        A = np.vstack([chi]).T
        sol = np.linalg.lstsq(A, elevation)
        m = sol[0]
        WRSS = sol[1]
        
        return (m, WRSS)
    
    def steves_best_guess_for_fit_chi(outlet,flow_direction, elevation, area, minimum_area=1E7, min_pathlength_in_cells = min_pathlength_in_cells):
        (area, elevation, de) = area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, minimum_area=minimum_area, min_pathlength_in_cells = min_pathlength_in_cells)
        
        def fit_by_making_chi(area, elevation, de, theta):
            chi = chi_for_profile(area ,de , theta)
            mean_elevation = np.mean(elevation)
            SS = np.sum(np.power(elevation - mean_elevation, 2))
            # find best fit ks and WRSS
            ks_wrss = (best_ks_with_wrss(chi, elevation), SS)
            # return these values for this tributary
            return ks_wrss
            # This is my lists of lists
        rez = []
        rez_list=[]
        import scipy.optimize
        for (area_all, elevation_all, de_all) in zip(area, elevation, de):    
            chi_ks_all = lambda theta: fit_by_making_chi(area_all, elevation_all, de_all, theta)[0][1]
            try:
                (xopt, _, _, _, warnflag) = scipy.optimize.fmin(chi_ks_all, np.array([0.5]), (), 1E-5, 1E-5, 100, 200, True, True, 0, None)
                ((ks, WRSS), SS) = fit_by_making_chi(area_all, elevation_all, de_all, np.array([xopt[0]]))
                ks = ks[0]
                WRSS = WRSS[0] 
            except:
                SS = 1.0
                WRSS = 0.0
                xopt = [0.0]
                ks = 0.0
                warnflag = 1
            R2 = 1 - (WRSS / SS) #calculate R2 value
            if warnflag == 1 or warnflag == 2:
                R2 = 0.0
                            
            theta_all = xopt[0]
            rez.append({'theta': theta_all,
                            'R2': R2,
                            'ks': ks,
                            'area': np.max(area_all)}) #LIST OF DICTIONARIES GOES HERE
            rez_list.append([theta_all,R2,ks,np.max(area_all)])

        return rez
    return steves_best_guess_for_fit_chi(outlet,flow_direction, elevation, area, minimum_area=minimum_area, min_pathlength_in_cells=min_pathlength_in_cells)

def best_ks_wrss_for_outlet(outlet, flow_direction, elevation, area, minimum_area = 1E7, theta = 0.4):

    def chi_for_profile(area, de, theta):
        chi = []
        chi_value = 0.0
        for (a, d) in zip(area, de):
            chi += [chi_value]
            chi_value += (1/a)**theta[0] * d
        return chi
    
    def best_ks_with_wrss(chi, elevation):
        A = np.vstack([chi]).T
        sol = np.linalg.lstsq(A, elevation)
        m = sol[0]
        WRSS = sol[1]
        
        return (m, WRSS)
    
    def steves_best_guess_for_fit_chi(outlet,flow_direction, elevation, area, minimum_area=1E7):
        (area, elevation, de) = area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, minimum_area=minimum_area)
        
        def fit_by_making_chi(area, elevation, de, theta):
            chi = chi_for_profile(area ,de , theta)
            mean_elevation = np.mean(elevation)
            SS = np.sum(np.power(elevation - mean_elevation, 2))
            # find best fit ks and WRSS
            ks_wrss = (best_ks_with_wrss(chi, elevation), SS)
            # return these values for this tributary
            return ks_wrss
            # This is my lists of lists
        rez = []
        rez_list=[]
        for (area_all, elevation_all, de_all) in zip(area, elevation, de):    
            try:
                ((ks, WRSS), SS) = fit_by_making_chi(area_all, elevation_all, de_all, np.array([theta]))
                ks = ks[0]
                WRSS = WRSS[0] 
            except:
                SS = 1.0
                WRSS = 0.0
                ks = 0.0
            
            R2 = 1 - (WRSS / SS) #calculate R2 value
                            
            rez.append({'R2': R2,
                        'ks': ks,
                        'area': np.max(area_all)}) #LIST OF DICTIONARIES GOES HERE
            rez_list.append([R2,ks,np.max(area_all)])

        return rez
    return steves_best_guess_for_fit_chi(outlet,flow_direction, elevation, area, minimum_area=minimum_area)

def best_ks_theta(outlet, flow_direction, elevation, area, minimum_area):
 #def best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, theta, minimum_area):
#         (area, elevation, de) = area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, theta, minimum_area)
#         area = area[0]
#         elevation = elevation[0]
#         de = de[0]
#         chi = chi_for_profile(area, de, theta)
#         mean_elevation = np.mean(elevation)
#         SS = np.sum(np.power(elevation - mean_elevation, 2))
#         return (best_ks_with_wrss(chi, elevation), SS)
    
#     def best_ks_theta_wrss_for_tribs(outlet, flow_direction, elevation, area, theta, minimum_area):
#         (area, elevation, de) = area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, theta, minimum_area)
#         area = area[1:] #assumes that all the rest of the list 
#         # are the tributaries and thus puts all the areas, eleveations, and de into a list
#         elevation = elevation[1:]
#         de = de[1:]
#         chi = []
#         for (area_profile, de_profile) in zip(area, de): #PLEASE EXPLAIN
#             chi.append(chi_for_profile(area_profile, de_profile, theta))
#         #chi = [c for sublist in chi for c in sublist] # puts together all the separate lists into one list for chi
#         #elevation = [e for sublist in elevation for e in sublist] # puts together all the separate lists into one list for elevation
#         mean_elevation = np.mean(elevation)
#         SS = np.sum(np.power(elevation - mean_elevation, 2))
#         dobbs = (best_ks_with_wrss(chi, elevation), SS)
#         return dobbs  # calculates the best ks for the tributaries
    # Currently, the chi and elevation data are put into the same list, and then theta and ks values are
    # approximated from there. I need to figure out how to keep these values separated so that I can calculate
    # individual ks and theta values for each tributary. 
    # WILLS IDEA: instead of running dobbs, loop each set of chis for each tributary, run dobbs and store it in a list and return that set of lists. loop through that list 
    # using the optimizatyion and store in a new dictionary
    
     #import scipy.optimize
#     chi_ks_mainstem = lambda theta: best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, theta, minimum_area)[0][1]
#     chi_ks_tribs = lambda theta: best_ks_theta_wrss_for_tribs(outlet, flow_direction, elevation, area, theta, minimum_area)[0][1]
#     # we use lambda when we require a nameless function for a short period of time
#     try:
#         (xopt, _, _, _, warnflag) = scipy.optimize.fmin(chi_ks_mainstem, np.array([0.5]), (), 1E-5, 1E-5, 100, 200, True, True, 0, None)
#         ((ks, WRSS), SS) = best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, np.array([xopt[0]]), minimum_area)
#         ks = ks[0]
#         WRSS = WRSS[0]
#         # fmin minimizes a function using the downhill simplex algorithm 
#     except:
#         SS = 1.0
#         WRSS = 0.0
#         xopt = [0.0]
#         ks = 0.0
#         warnflag = 1
#  
#     R2 = 1 - (WRSS / SS) #calculate R2 value
#     if warnflag == 1 or warnflag == 2:
#         R2 = 0.0
#     theta_mainstem = xopt[0]
#     R2_mainstem = R2
#     ks_mainstem = ks
#     
#     try:
#         (xopt, _, _, _, warnflag) = scipy.optimize.fmin(chi_ks_tribs, np.array([0.5]), (), 1E-5, 1E-5, 100, 200, True, True, 0, None)
#         ((ks, WRSS), SS) = best_ks_theta_wrss_for_tribs(outlet, flow_direction, elevation, area, np.array([xopt[0]]), minimum_area)
#         ks = ks[0]
#         WRSS = WRSS[0]
#     except:
#         SS = 1.0
#         WRSS = 0.0
#         xopt[0] = 0.0
#         ks = 0.0
#         
#     R2 = 1 - (WRSS / SS)
#     if warnflag == 1 or warnflag == 2:
#         R2 = 0.0
#     theta_tribs = xopt[0]
#     R2_tribs = R2
#     ks_tribs = ks
#     
#     return{'mainstem': {'theta': theta_mainstem,
#                         'R2': R2_mainstem,
#                         'ks': ks_mainstem},
#            'tributaries': {'theta': theta_tribs,
#                            'R2': R2_tribs,
#                            'ks': ks_tribs}
#            }

    def chi_for_profile(area, de, theta):
        # area, de,defined from area_elevation_for_mainstem_and_tributaries
        chi = [] # chi initially defined as an empty list
        chi_value = 0.0
        for (a, d) in zip(area, de):
            chi += [chi_value]
            chi_value += (1/a)**theta[0] * d #calculates chi for every area and de value previously calculated
        return chi

    def best_ks_with_wrss(chi, elevation):
        A = np.vstack([chi]).T # stacks the array of the chi list vertically
        sol = np.linalg.lstsq(A, elevation) # returns the least-squares solution to a linear matrix equation of A and elevation
        m = sol[0] # m is the solution for the mainstem
        WRSS = sol[1]  # not sure what WRSS is...      
        return (m, WRSS)

    def best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, theta, minimum_area):
        (area, elevation, de) = area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, theta, minimum_area)
        # finds area, elevation, and de using area_elevation_for_mainstem_and_tributaries
        area = area[0] # the mainstem values will be the first values of the list so it defines area, elevation, and de as [0] 
        elevation = elevation[0]
        de = de[0]
        chi = chi_for_profile(area, de, theta) # calculates chi values for the mainstem
        mean_elevation = np.mean(elevation) 
        SS = np.sum(np.power(elevation - mean_elevation, 2))
        return (best_ks_with_wrss(chi, elevation), SS)

    import scipy.optimize
    chi_ks_mainstem = lambda theta: best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, theta, minimum_area)[0][1]
    try:
        seed = np.random.rand(1)
        (xopt, _, _, _, warnflag, sols) = scipy.optimize.fmin(chi_ks_mainstem, seed, (), 1E-5, 1E-5, 100, 200, True, True, True, None)
        ((ks, WRSS), SS) = best_ks_theta_wrss_for_mainstem(outlet, flow_direction, elevation, area, np.array([xopt[0]]), minimum_area)
        ks = ks[0]
        WRSS = WRSS[0]
        print(sols)
    except OSError:
        SS = 1.0
        WRSS = 0.0
        xopt = [0.0]
        ks = 0.0
        warnflag = 1
 
    R2 = 1 - (WRSS / SS)
    if warnflag == 1 or warnflag == 2:
        R2 = 0.0
    theta = xopt[0]
    return {'theta': theta,
            'R2': R2,
            'ks': ks }
        
def map_chi_profiles(elevation, flow_direction, area, outlet, minimum_area = 1.0E6, theta = 0.5, start_at = 0.0, downstream = True, Ao = 1.0E6):
    
    ((row, col),) = elevation._xy_to_rowscols((outlet,))
    base_elevation = elevation[row,col]
    
    return_map = {}
    
    def map_ld_link(current_chi, ld_list, downstream_sign, minimum_area):
        (current_row, current_column) = ld_list['index']
        
        if ld_list.get('next') is None:
            return
        for next_list in ld_list['next']:
            
            if next_list['area'] >= minimum_area:
                
                index = next_list['index']
                (next_row, next_column) = index
                if (current_row != next_row) & (current_column != next_column):
                    next_chi = current_chi + (1 / np.array(next_list['area']))**theta * np.array(ld_list['de']*1.414*downstream_sign)
                else:
                    next_chi = current_chi + (1 / np.array(next_list['area']))**theta * np.array(ld_list['de']*downstream_sign)
                next_elevation = next_list['elevation'] - base_elevation
                return_map[index] = (next_chi, next_elevation)
                map_ld_link(next_chi, next_list, downstream_sign, minimum_area)
    
    import dem as d                            
    mean_pixel_dimension = d.BaseSpatialGrid()
    mean_pixel_dimension._copy_info_from_grid(area, True)
    mean_pixel_dimension._griddata = area._mean_pixel_dimension()
                
    ld_list = flow_direction.map_values_to_recursive_list(outlet, elevation = elevation, area = area, de = mean_pixel_dimension)
    current_chi = start_at
    if downstream:
        downstream_sign = 1.0
    else:
        downstream_sign = -1.0
    
    map_ld_link(current_chi, ld_list, downstream_sign, minimum_area)
    
    return return_map

