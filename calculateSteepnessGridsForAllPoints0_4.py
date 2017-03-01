import dem as d
import numpy as np

def calculate_ks(prefix, outlet_distance, suffix):
    chi = d.Chi.load(prefix + '_chi_' + outlet_distance + suffix)
    relief = d.ScaledRelief.load(prefix + '_relief_' + outlet_distance + suffix)
    ks = d.BaseSpatialGrid()
    ks._copy_info_from_grid(relief)
    ks._griddata = np.divide(ks._griddata, chi._griddata)
    ks.save(prefix + '_ks_' + outlet_distance + suffix)
    
suffix = '_0_4_1000000'
prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']
outlet_distances = ['50000', '100000', '200000', '400000']

for prefix in prefixes:
    print('Prefix: ' + prefix)
    for outlet_distance in outlet_distances:
        print('Outlet distance: ' + outlet_distance)
        calculate_ks(prefix, outlet_distance, suffix)


        