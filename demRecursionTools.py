import numpy as np

def extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, xo = 500.0):

    Ao = np.power(xo, 2.0)
    if ld_list['area'] >= Ao:
        elevation_f = ld_list['elevation'] - base_elevation
        elevation = elevation + [elevation_f]
        index = ld_list['index']
        chi_f = chi_o + (1 / ld_list['area'])**theta * ld_list['distance_scale'] * de[index[0], index[1]]
        chi = chi + [chi_f]
        if ld_list.get('next') is not None:
            for next_list in ld_list['next']:
                elevation, chi = extract_chi_elevation_values(next_list, de, theta, chi_f, elevation, chi, base_elevation)
    
    return elevation, chi    

def chi_elevation(ld_list, de, theta, xo = 500.0):
    
    chi_o = 0
    elevation = []
    chi = []
    base_elevation = ld_list['elevation']
    
    return extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, xo)

def best_ks_with_wrss_list(ld_list, de, theta, xo = 500):
    
    e, c = chi_elevation(ld_list, de, theta, xo)
    A = np.vstack([c]).T
    sol = np.linalg.lstsq(A, e)
    m = sol[0]
    WRSS = sol[1]    
    
    return (m, WRSS)

def uninformative_SS_list(ld_list, de, xo = 500):
    
    e, c = chi_elevation(ld_list, de, 0.5, xo)
    mean_elevation = np.mean(e)
    return np.sum(np.power(e-mean_elevation, 2))

def best_ks_and_theta_with_wrss(elevation, flow_direction, area, outlet, xo = 500):
    
    ld_list = flow_direction.map_values_to_recursive_list(outlet, area = area, elevation = elevation)
    de = area._mean_pixel_dimension()
    
    import scipy.optimize
    chi_ks = lambda theta: best_ks_with_wrss_list(ld_list, de, theta, xo)[1]
    if len(chi_ks(0.5)) == 0:
        return (0, 0, 0)
    xopt = scipy.optimize.fmin(func=chi_ks, x0=np.array([0.5]))
    (m, WRSS) = best_ks_with_wrss_list(ld_list, de, xopt[0], xo)
    SS = uninformative_SS_list(ld_list, de, xo)
    
    return (m, xopt[0], 1 - (WRSS / SS))