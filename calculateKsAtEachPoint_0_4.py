#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import dem as d
import numpy as np

suffix = '0_4'

for prefix in prefixes:
    
    import copy
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_' + suffix)
    relief = d.ScaledRelief.load(prefix + '_relief_250000_' + suffix)
    
    ks = copy.deepcopy(ksi)
    ks._griddata = np.divide(relief._griddata, ksi._griddata)
    
    ks.save(prefix + '_ks_250000_' + suffix)
    

        

    
