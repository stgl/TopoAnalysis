#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import matplotlib
matplotlib.use('Agg')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids
import numpy as np
import demMethods as dm

a = [0, 70000, 0, 2400000]
suffix = '0_4'

dx = 100.0
dy = 10000.0

contours = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

x_bins = np.arange(a[0],a[1],dx)
y_bins = np.arange(a[2],a[3],dy)

for prefix in prefixes:
    print(prefix)    
    ksi = d.Ksi.load(prefix + '_ksi_2000000_' + suffix)
    relief = d.ScaledRelief.load(prefix + '_relief_2000000_' + suffix)
    
    H, xedges, yedges = dm.create_density(ksi,relief,x_bins,y_bins)
    
    if 'H_tot' not in locals():
        H_tot = H.copy()
    else:
        H_tot = H_tot + H

H_tot = np.flipud(H_tot.T)
H_tot = H_tot / np.sum(H_tot) / dx / dy
v = np.ndarray.flatten(H_tot)
v = np.sort(v)

vc = np.cumsum(v) * dx * dy

plt.figure()
plt.imshow(np.log10(H_tot), extent = a)
plt.colorbar()
plt.ion()
plt.axis('normal')

for contour in contours:
    i = np.where(vc >= contour)
    i = np.min(i)
    contour_value = v[i]
    plt.contour(np.flipud(H_tot) > contour_value, levels = [0], extent = a)

plt.savefig('density' + suffix + '.eps')


    
