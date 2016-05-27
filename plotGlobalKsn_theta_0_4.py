#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import matplotlib
matplotlib.use('Agg')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids
import numpy as np

a = [0, 35000, 0, 1200000]
suffix = '0_4'

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_' + suffix)
    relief = d.ScaledRelief.load(prefix + '_relief_250000_' + suffix)
    
    fig = plt.figure()
    ksi_vec, relief_vec = plotGrids(ksi, relief, 'k.', rasterized = True, markersize=5.0)
    plt.axis(a) 
    plt.savefig(prefix + suffix + '.png',dpi=600)
    plt.close()
    
    if 'all_ksi_vec' not in locals():
        all_ksi_vec = ksi_vec
        all_relief_vec = relief_vec
    else:
        all_ksi_vec = np.append(all_ksi_vec,ksi_vec)
        all_relief_vec = np.append(all_relief_vec, relief_vec)
    print('loaded ' + prefix)
        
fig = plt.figure()

plt.plot(all_ksi_vec, all_relief_vec, 'k.', rasterized = True, markersize=5.0)
plt.axis(a)
plt.savefig('all_' + suffix + '.png',dpi=600)
    
