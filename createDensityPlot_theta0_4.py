#! /usr/bin/env python2.7

prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import matplotlib
matplotlib.use('Agg')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids
import numpy as np
import demMethods as dm

kss_to_report = [100.0, 150.0, 200.0]

a = [0, 20000, 0, 2000000]
suffix = '0_4'

dx = 100.0
dy = 10000.0

contours = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
basin_lengths = [50000, 100000, 200000, 400000]

x_bins = np.arange(a[0],a[1],dx)
y_bins = np.arange(a[2],a[3],dy)

chi_vec = np.array([])
relief_vec = np.array([])

f = open('stats_for_concavity0_4.txt','w')

for basin_length in basin_lengths:
    for prefix in prefixes:
        print(prefix)    
        chi = d.GeographicChi.load(prefix + '_chi_' + str(basin_length) + '_' + suffix + '_1000000')
        relief = d.ChiScaledRelief.load(prefix + '_relief_' + str(basin_length) + '_' + suffix + '_1000000')
        this_chi, this_relief = dm.extract_values_from_grid(chi, relief, ignore_zeros=True)
        chi_vec = np.concatenate((chi_vec, this_chi))
        relief_vec = np.concatenate((relief_vec, this_relief))
        ks_vec = this_relief / this_chi
        for ks_to_report in kss_to_report:
            i = np.where(ks_vec > ks_to_report)
            f.write('Fraction of points exceeding threshold: ' + str(ks_to_report) + '; for dataset: ' + prefix + '; basin length: ' + str(basin_length) + '; concavity: ' + suffix + ': ' + str(float(len(i[0]))/float(len(this_chi))) + '\n')
        H, xedges, yedges = dm.create_density(this_chi,this_relief,x_bins,y_bins)
        H = np.flipud(H.T)
        H = H / np.sum(H) / dx / dy
        v = np.ndarray.flatten(H)
        v = np.sort(v)
        vc = np.cumsum(v) * dx * dy
        plt.figure()
        plt.imshow(np.log10(H), extent = a)
        plt.colorbar()
        plt.ion()
        plt.axis('normal')
        for contour in contours:
            i = np.where(vc >= contour)
            i = np.min(i)
            contour_value = v[i]
            plt.contour(np.flipud(H) > contour_value, levels = [0], extent = a)
    
            plt.savefig('density_' + prefix + '_' + str(basin_length) + '_' + suffix + '.eps')
    
    ks_vec = relief_vec / chi_vec
    for ks_to_report in kss_to_report:
        i = np.where(ks_vec > ks_to_report)
        f.write('Fraction of points exceeding threshold: ' + str(ks_to_report) + '; basin length: ' + str(basin_length) + '; concavity: ' + suffix + ': ' + str(len(i[0])/len(this_chi)) + '\n')

    H, xedges, yedges = dm.create_density(chi_vec,relief_vec,x_bins,y_bins)
    H = np.flipud(H.T)
    H = H / np.sum(H) / dx / dy
    v = np.ndarray.flatten(H)
    v = np.sort(v)
    
    vc = np.cumsum(v) * dx * dy
    
    plt.figure()
    plt.imshow(np.log10(H), extent = a)
    plt.colorbar()
    plt.ion()
    plt.axis('normal')
    
    for contour in contours:
        i = np.where(vc >= contour)
        i = np.min(i)
        contour_value = v[i]
        plt.contour(np.flipud(H) > contour_value, levels = [0], extent = a)
    
    plt.savefig('density_' + str(basin_length) + '_'+ suffix + '.eps')


    
