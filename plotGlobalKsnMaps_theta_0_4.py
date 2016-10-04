#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')
threshold = 100

import matplotlib
matplotlib.use('Agg')

import dem as d
import copy

a = [0, 35000, 0, 1200000]
suffix = '0_4'

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_' + suffix)
    relief = d.ScaledRelief.load(prefix + '_relief_250000_' + suffix)
    ks = copy.deepcopy(ksi)
    ks._griddata = ((relief._griddata / ksi._griddata) >= threshold).astype(int)
    
    ks.save(prefix + '_ks_mask_250000_' + suffix + "_" + str(threshold))
        
    
