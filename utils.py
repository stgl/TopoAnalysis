def calc_ks_for_outlet(outlet, theta, **kwargs):
    fd = kwargs['flow_direction']
    xo = kwargs.get('xo') if kwargs.get('xo') is not None else 500.0
    kwargs.pop('xo')
    plot = kwargs.pop('plot', None)
    
    de = fd._mean_pixel_dimension()
    ld_list = fd.map_values_to_recursive_list(outlet, **kwargs)
    
    from demRecursionTools import best_ks_with_r2_list
    from numpy import array
    ret = best_ks_with_r2_list(ld_list, de, array([theta]), xo = xo)

    if plot is not None and plot == theta:
        import matplotlib.pylab as plt
        from demRecursionTools import chi_elevation
        import numpy as np
        e, c = chi_elevation(ld_list, de, array([theta]), xo=xo)
        plt.figure()
        plt.plot(c, e, 'k.')
        plt.hold(True)
        (ks, _) = ret
        plt.plot([0, np.max(c)],[0, ks*np.max(c)], 'k-')
        
    return ret
        
    