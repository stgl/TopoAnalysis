#! /usr/bin/env python2.7

import pickle

ksdict = pickle.load(open('ks.p','rb'))

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

# Calculate fit to data using Venezuela and Brazil:

dr_regress = list()
ks_regress = list()
dr = np.load('costarica/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('guatemala/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('taiwan/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
ks_regress = ks_regress + [np.log10(x[0]) for x in ksdict['costarica']['0_4']] + [np.log10(y[0]) for y in ksdict['guatemala']['0_4']] + [np.log10(z[0]) for z in ksdict['taiwan']['0_4']]

dr_regress_check = list()
ks_regress_check = list()

for (dr_el, ks_el) in zip(dr_regress, ks_regress):
     if np.isfinite(dr_el) and np.isfinite(ks_el):
         dr_regress_check.append(dr_el)
         ks_regress_check.append(ks_el)

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check,2))

sm = m

print('slope: ' + str(sm))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))

# Calculate fit to rest of data:

dr_regress = list()
ks_regress = list()
dr = np.load('brazil/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('venezuela/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
ks_regress = ks_regress + [np.log10(x[0]) for x in ksdict['brazil']['0_4']] + [np.log10(y[0]) for y in ksdict['venezuela']['0_4']]

dr_regress_check = list()
ks_regress_check = list()

for (dr_el, ks_el) in zip(dr_regress, ks_regress):
     if np.isfinite(dr_el) and np.isfinite(ks_el):
         dr_regress_check.append(dr_el)
         ks_regress_check.append(ks_el)

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check,2))

sm = m

print('slope: ' + str(sm))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))

for suffix in suffixes:
    for prefix in prefixes:
        this_ks = ksdict[prefix][suffix[0]+'_'+suffix[1]]
        ks = [ks[0] for ks in this_ks]
        dr = np.load(prefix + '/'  + 'dr.npy')
        dr_sig = np.load(prefix + '/' + 'drstd.npy')
        
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
            
