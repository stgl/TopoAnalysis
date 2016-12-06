#! /usr/bin/env python2.7

import pickle

prefixes = ('brazil', 'costarica', 'guatemala', 'taiwan', 'venezuela')

colors = {'brazil': 'b',
          'costarica': 'c',
          'guatemala': 'g',
          'taiwan': 'm',
          'venezuela': 'r'}
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import numpy as np
import denudationRateAnalysis as dra

suffixes = ('04', '05', '06')

sm = 1.36
b = 3.29
# R2 is 0.96

ksdict = pickle.load(open('ks.p','rb'))

for suffix in suffixes:
    for prefix in prefixes:
        ks = ksdict[prefix][suffix[0]+'_'+suffix[1]]
        #ks = np.load(prefix + '_ks' + suffix + '.npy')
        dr = np.load(prefix + '/' + prefix + '_dr.npy')
        dr_sig = np.load(prefix + '/' + prefix + '_drstd.npy')
        
        plt.loglog(dr,ks,colors[prefix] + '.')
        for (m,s,k) in zip(dr,dr_sig,ks):
            plt.loglog([m-s, m+s],[k, k], colors[prefix] + '-')
    
    if suffix == '04':
        dra.plot_stock_and_montgomery()
        x = (1.0e-4, 1.0e1)
        y = (10**b * x[0]**sm, 10**b * x[1]**sm)
        plt.loglog(x,y,'k-')
        
    plt.grid()
    plt.axis('image')
    plt.axis([1e-4,1e1,1e-1,1e4])    
    plt.savefig('alldata_dr_steepness_' + suffix + '.eps')
    plt.close()
            
