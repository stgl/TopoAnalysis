#! /usr/bin/env python2.7

import matplotlib
matplotlib.use('Agg')
import denudationRateAnalysis as dra
from matplotlib import pyplot as plt
import numpy as np

portenga_data = dra.read_csv('portenga_filtered_steepness.csv')
 
dr = [r[3] for r in portenga_data]
dr_sig = [r[4] for r in portenga_data]

ks = dict()
r2 = dict()

ks['0.4'] = [np.float64(r[10].replace('[','').replace(']','')) for r in portenga_data]
r2['0.4'] = [np.float64(r[11].replace('[','').replace(']','')) for r in portenga_data]

ks['0.5'] = [np.float64(r[12].replace('[','').replace(']','')) for r in portenga_data]
r2['0.5'] = [np.float64(r[13].replace('[','').replace(']','')) for r in portenga_data]

ks['0.6'] = [np.float64(r[14].replace('[','').replace(']','')) for r in portenga_data]
r2['0.6'] = [np.float64(r[15].replace('[','').replace(']','')) for r in portenga_data]


keys = ['0.4', '0.5', '0.6']

maxks = 0

for key in keys:
    
    plt.figure()
    for (m, s, k, r) in zip(np.float64(dr)/1000.0, np.float64(dr_sig)/1000.0, ks[key], r2[key]):
        
        m = np.array(m)
        s = np.array(s)
        k = np.array(k)
        r = np.array(r)

	if key == '0.4' and k > maxks:
		maxks = k

        if r >= 0.9:
		plt.loglog((m-s, m+s), (k,k), 'k-')
        	plt.loglog(m,k,'k.') 
        
        if (r < 0.9) & (r >= 0.7):
		plt.loglog((m-s, m+s), (k,k), 'b-')
        	plt.loglog(m,k,'b.')
        
        if r < 0.7:
        	plt.loglog((m-s, m+s), (k,k), 'r-')
        	plt.loglog(m,k,'r.')
    
    if key == '0.4':
        dra.plot_stock_and_montgomery()
        
    plt.axis('image')
    plt.axis([1e-4,1e1,1e-1,1e4])

    plt.grid()
    plt.savefig('dr_steepness_' + key.replace('.','_') + '.eps')

print(maxks)
    
