#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import matplotlib
matplotlib.use('Agg')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids
import numpy as np

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_0_4')
    relief = d.ScaledRelief.load(prefix + '_relief_250000_0_4')
    
    fig = plt.figure()
    ksi_vec, relief_vec = plotGrids(ksi, relief, 'k.', rasterized = True)
    plt.savefig(prefix + '0_4.png',dpi=600)
    plt.close()
    
    if 'all_ksi_vec' not in locals():
        all_ksi_vec = ksi_vec
        all_relief_vec = relief_vec
    else:
        all_ksi_vec = np.append(all_ksi_vec,ksi_vec)
        all_relief_vec = np.append(all_relief_vec, relief_vec)
    print('loaded ' + prefix)
        
fig = plt.figure()

plt.plot(all_ksi_vec, all_relief_vec, 'k.', rasterized = True)
plt.savefig('all_0_4.png',dpi=600)
    
