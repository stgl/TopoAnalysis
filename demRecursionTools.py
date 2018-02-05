import numpy as np

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

def chi_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, theta = 0.5, minimum_area = 1.0E7):
    
    import dem as d
    mean_pixel_dimension = d.BaseSpatialGrid()
    mean_pixel_dimension._copy_info_from_grid(area, True)
    mean_pixel_dimension._griddata = area._mean_pixel_dimension()
    
    ld_list = flow_direction.map_values_to_recursive_list(outlet, elevation = elevation, area = area, de = mean_pixel_dimension)
    
    area = [ld_list['area']]
    elevation = [ld_list['elevation']]
    de = [ld_list['de']*ld_list['distance_scale']]
    tributary_ld = []
    
    def get_elevations_and_areas(ld_list, area, elevation, de, tributary_ld, minimum_area_to_consider):
        maximum_area = 0.0
        for next in ld_list['next']:
            if (next['area'] > minimum_area_to_consider) and (next['area'] > maximum_area):
                maximum_area = next['area']
        
        for next in ld_list['next']:
            if next['area'] == maximum_area:
                area += [next['area']]
                elevation += [next['elevation']]
                de += [next['de'] * next['distance_scale']]
                (area, elevation, de, tributary_ld) = get_elevations_and_areas(next, area, elevation, de, tributary_ld, minimum_area_to_consider)
            elif next['area'] > minimum_area_to_consider:
                tributary_ld.append(next)   
        return (area, elevation, de, tributary_ld)
    
    (area, elevation, de, tributary_ld) = get_elevations_and_areas(ld_list, area, elevation, de, tributary_ld, minimum_area)
    
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
            area.append(this_area)
            elevation.append(this_elevation)
        tributary_ld = next_tributary_ld

    return_chi = []
    return_elevation = []
    for (this_area, this_elevation, this_de) in zip(area, elevation, de):
        print(this_area)
        print(this_elevation)
        this_return_elevation = []
        for elevation_value in this_elevation:
            this_return_elevation += [elevation_value - this_elevation[0]]
        return_elevation.append(this_return_elevation)
        this_chi_value = 0.0
        this_chi = []
        for (area_value, de_value) in zip(this_area, this_de):
            if area_value == this_area[0]:
                this_chi += [0.0]
            else:
                this_chi_value += (1 / area_value)**theta * de_value
                this_chi += [this_chi_value]
        return_chi.append(this_chi)    
        
    return (return_chi, elevation)
    
                
        
        