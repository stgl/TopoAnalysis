#! /usr/bin/env python2.7

import matplotlib
matplotlib.use('Agg')
import denudationRateAnalysis as dra

portenga_data = dra.read_csv('portengadata.csv')
del(portenga_data[0])

dr = [r[3] for r in portenga_data]
dr_sig = [r[4] for r in portenga_data]

from matplotlib import pyplot as plt
import numpy as np

suffixes = ['0_4', '0_5', '0_6']

for suffix in suffixes:
    data = np.load('ks_area_data_' + suffix + '.npz')
    ks = data['ksn_vec']
    plt.figure()
    plt.loglog(np.float64(dr)/1000.0,ks,'k.')
    for (m, s, k) in zip(np.float64(dr)/1000.0, np.float64(dr_sig)/1000.0, ks):
        plt.loglog((m-s, m+s), (k,k), 'k-')
    
    
    plt.axis('image')
    plt.axis([1e-4,1e1,1e-1,1e4])

    plt.grid()
    plt.savefig('dr_steepness_' + suffix + '.eps')
    