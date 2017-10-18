#! /usr/bin/env python2.6

import denudationRateAnalysis as dra
import dem as d
import numpy as np
import pickle

prefixes = ['brazil', 'costarica', 'guatemala', 'taiwan', 'venezuela']
suffixes = ['0_4', '0_5', '0_6']
theta_values = [0.4, 0.5, 0.6]

ks = dict()

Ao = 1000000.0
A_cutoff = 500000.0 

for prefix in prefixes:
    print(prefix)
    d8 = d.FlowDirectionD8.load(prefix + "/" + prefix + '_flow_direction')
    area = d.Area.load(prefix + '/' + prefix + '_area')
    locs = np.load(prefix + '/' + prefix +  '_sample_locations.npy')
    local_dict = dict()
    elevation = d.Elevation.load(prefix + '/' + prefix + '_elevation')
    filled = d.FilledElevation(elevation = elevation)
    d8.sort()
    d8._sort_indexes = filled.sort(reverse = False, force = True)
    elevation = d.Elevation.load(prefix + "/" + prefix + '_elevation')
    for theta in theta_values:
        chi = d.Chi(area = area, flow_direction = d8, theta = theta, Ao = Ao, outlets = locs)
        relief = d.ChiScaledRelief(elevation = elevation, flow_direction = d8, theta = theta, Ao = Ao, outlets = locs)
        local_dict[str(theta).replace('.','_')] = dra.calculate_ks_for_sample(locs, d8, chi, relief, area, Ao, A_cutoff=A_cutoff)
    ks[prefix] = local_dict
    print(local_dict['0_4'])

pickle.dump( ks, open( "ks.p", "wb" ) )
