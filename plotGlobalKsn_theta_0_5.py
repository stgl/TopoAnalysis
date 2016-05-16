#! /usr/bin/env python

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids
import numpy as np

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_0_0_5')
    relief = d.ScaledRelief.load(prefix + '_relief_250000_0_0_5')
    
    fig = plt.figure()
    ksi_vec, relief_vec = plotGrids(ksi, relief, 'k.', rasterized = True)
    #plt.savefig(prefix + '0_5.png',dpi=600)
    #plt.close()
    
    if 'all_ksi_vec' not in locals():
        all_ksi_vec = ksi_vec
        all_relief_vec = relief_vec
    else:
        all_ksi_vec = np.append(all_ksi_vec,ksi_vec)
        all_relief_vec = np.append(all_relief_vec, relief_vec)
    print('loaded ' + prefix)
        
fig = plt.figure()

plt.plot(all_ksi_vec, all_relief_vec, 'k.', rasterized = True)
plt.savefig('all_0_5.png',dpi=600)
    