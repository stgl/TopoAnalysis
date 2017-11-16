#! /usr/bin/env python2.6

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

print(np.power(10,ks_regress_check))
print(np.power(10,dr_regress_check))

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check-np.mean(ks_regress_check),2))

sm = m

syx = np.sqrt(R / (len(dr_regress_check)-2))
sb = syx / np.sqrt(np.power(np.sum(dr_regress_check - np.mean(dr_regress_check)),2))

print('slope: ' + str(sm))
print('slope uncertainty: ' + str(sb))
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

print(np.power(10,ks_regress_check))
print(np.power(10,dr_regress_check))

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check-np.mean(ks_regress_check),2))

sm = m

syx = np.sqrt(R / (len(dr_regress_check)-2))
sb = syx / np.sqrt(np.power(np.sum(dr_regress_check - np.mean(dr_regress_check)),2))

print('slope: ' + str(sm))
print('slope uncertainty: ' + str(sb))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))

print('regression for ALL data: 0.4')

dr_regress = list()
ks_regress = list()
dr = np.load('brazil/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('venezuela/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('costarica/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('guatemala/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('taiwan/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
ks_regress = ks_regress + [np.log10(x[0]) for x in ksdict['brazil']['0_4']] + [np.log10(y[0]) for y in ksdict['venezuela']['0_4']] + [np.log10(x[0]) for x in ksdict['costarica']['0_4']] + [np.log10(y[0]) for y in ksdict['guatemala']['0_4']] + [np.log10(z[0]) for z in ksdict['taiwan']['0_4']]

dr_regress_check = list()
ks_regress_check = list()

for (dr_el, ks_el) in zip(dr_regress, ks_regress):
     if np.isfinite(dr_el) and np.isfinite(ks_el):
         dr_regress_check.append(dr_el)
         ks_regress_check.append(ks_el)

print(np.power(10,ks_regress_check))
print(np.power(10,dr_regress_check))

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check-np.mean(ks_regress_check),2))

sm = m

syx = np.sqrt(R / (len(dr_regress_check)-2))
sb = syx / np.sqrt(np.power(np.sum(dr_regress_check - np.mean(dr_regress_check)),2))

print('slope: ' + str(sm))
print('slope uncertainty: ' + str(sb))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))

print('regression for ALL data: 0.5')

dr_regress = list()
ks_regress = list()
dr = np.load('brazil/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('venezuela/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('costarica/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('guatemala/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('taiwan/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
ks_regress = ks_regress + [np.log10(x[0]) for x in ksdict['brazil']['0_5']] + [np.log10(y[0]) for y in ksdict['venezuela']['0_5']] + [np.log10(x[0]) for x in ksdict['costarica']['0_5']] + [np.log10(y[0]) for y in ksdict['guatemala']['0_5']] + [np.log10(z[0]) for z in ksdict['taiwan']['0_5']]

dr_regress_check = list()
ks_regress_check = list()

for (dr_el, ks_el) in zip(dr_regress, ks_regress):
     if np.isfinite(dr_el) and np.isfinite(ks_el):
         dr_regress_check.append(dr_el)
         ks_regress_check.append(ks_el)

print(np.power(10,ks_regress_check))
print(np.power(10,dr_regress_check))

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check-np.mean(ks_regress_check),2))

sm = m

syx = np.sqrt(R / (len(dr_regress_check)-2))
sb = syx / np.sqrt(np.power(np.sum(dr_regress_check - np.mean(dr_regress_check)),2))

print('slope: ' + str(sm))
print('slope uncertainty: ' + str(sb))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))



print('regression for ALL data: 0.6')

dr_regress = list()
ks_regress = list()
dr = np.load('brazil/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('venezuela/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('costarica/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('guatemala/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
dr = np.load('taiwan/dr.npy')
dr_regress = dr_regress + [np.log10(x) for x in dr]
ks_regress = ks_regress + [np.log10(x[0]) for x in ksdict['brazil']['0_6']] + [np.log10(y[0]) for y in ksdict['venezuela']['0_6']] + [np.log10(x[0]) for x in ksdict['costarica']['0_6']] + [np.log10(y[0]) for y in ksdict['guatemala']['0_6']] + [np.log10(z[0]) for z in ksdict['taiwan']['0_6']]

dr_regress_check = list()
ks_regress_check = list()

for (dr_el, ks_el) in zip(dr_regress, ks_regress):
     if np.isfinite(dr_el) and np.isfinite(ks_el):
         dr_regress_check.append(dr_el)
         ks_regress_check.append(ks_el)

print(np.power(10,ks_regress_check))
print(np.power(10,dr_regress_check))

A = np.vstack([dr_regress_check, np.ones(len(dr_regress_check))]).T

ks_regress_check = np.array(ks_regress_check)

result = np.linalg.lstsq(A, ks_regress_check)
(m, b) = result[0]
R = result[1]
print(R)
R2 = 1 - R / np.sum(np.power(ks_regress_check-np.mean(ks_regress_check),2))

sm = m

syx = np.sqrt(R / (len(dr_regress_check)-2))
sb = syx / np.sqrt(np.power(np.sum(dr_regress_check - np.mean(dr_regress_check)),2))

print('slope: ' + str(sm))
print('slope uncertainty: ' + str(sb))
print('intercept: ' + str(b))
print('R2 = ' + str(R2))



harel_dr = [0.51563227, 0.42379423, 1.6724333, 1.0277396, 0.76055958, 0.96075007, 0.95488516, 0.84357639, 0.4529031, 0.83152162, 0.97846715, 1.4010328, 0.67417162, 3.8464178, 0.39482937, 0.00764078, 0.012289158, 0.007522207, 0.009929992, 0.00527313, 0.009494532, 0.009138029, 0.011390154, 0.0108523, 0.010885823, 0.010237483, 0.008462175, 0.012757005, 0.009944075, 0.005874143]
harel_ks_05 = [403.1063345, 374.5994588, 357.7528502, 426.6259243, 449.4054637, 469.4254082, 454.7907542, 452.6106149, 441.1653558, 680.8821015, 434.0378082, 473.2591487, 309.9089028, 490.3473985, 354.8220746, 56.30693648, 75.89848786, 45.69839069, 37.12863245, 41.94146584, 53.7260818, 44.42326007, 90.44146107, 81.11600757, 115.9466899, 71.39814714, 65.2297482, 77.62707199, 76.65422994, 112.040401]


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
    
    if suffix == '05':
        plt.loglog(harel_dr, harel_ks_05, 'k.')
            
    plt.grid()
    plt.axis('image')
    plt.axis([1e-4,1e1,1e-1,1e4])    
    plt.savefig('alldata_dr_steepness_' + suffix + '.eps')
    plt.close()
            
