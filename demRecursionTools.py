import numpy as np

def extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, xo = 500.0):

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

def extract_profile_values(ld_list, xo = 500.0, items = ()):
    
    Ao = np.power(xo, 2.0)
    return_list = list()
    
    if ld_list['area'] >= Ao:
        this_list = list()
        for arg in items:
            this_list += [ld_list[arg]]
        print(this_list)
        return_list.append(this_list)
        if ld_list.get('next', None) is not None:
            return_list += extract_profile_values(ld_list['next'][0], xo = xo, items = items)
    
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

def uninformative_SS_list(ld_list, de, xo = 500):
    
    e, c = chi_elevation(ld_list, de, [0.5], xo=xo)
    mean_elevation = np.mean(e)
    return np.sum(np.power(e-mean_elevation, 2))

def best_ks_and_theta_with_wrss_list(ld_list, de, xo = 500):
    
    import scipy.optimize
    chi_ks = lambda theta: best_ks_with_wrss_list(ld_list, de, theta, xo=xo)[1]
    if len(chi_ks([0.5])) == 0:
        return (0, 0, 0)
    xopt = scipy.optimize.fmin(func=chi_ks, x0=np.array([0.5]))
    (m, WRSS) = best_ks_with_wrss_list(ld_list, de, [xopt[0]], xo)
    SS = uninformative_SS_list(ld_list, de, xo)
    
    return (m, xopt[0], 1 - (WRSS / SS))
    
def best_ks_and_theta_with_wrss(elevation, flow_direction_or_length, area, outlet, xo = 500):
    
    ld_list = flow_direction_or_length.map_values_to_recursive_list(outlet, area = area, elevation = elevation)
    de = area._mean_pixel_dimension()
    
    return best_ks_and_theta_with_wrss_list(ld_list, de, xo = xo)


    