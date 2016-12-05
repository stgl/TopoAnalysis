import denudationRateAnalysis as dra
import dem as d
import numpy as np
import pickle

prefixes = ['brazil', 'costarica', 'guatemala', 'taiwan', 'venezuela']
suffixes = ['0_4', '0_5', '0_6']

ks = dict()

Ao = 2000000.0

for prefix in prefixes:
    print(prefix)

    d8 = d.FlowDirectionD8.load(prefix + "/" + prefix + '_flow_direction')
    area = d.Area.load(prefix + '/' + prefix + '_area')
    locs = np.load(prefix + '/' + prefix +  '_sample_locations.npy')
    local_dict = dict()
    filled = d.FilledElevation.load(prefix + "/" + prefix + "_filled") 
    for suffix in suffixes:    
#        ksi = d.Ksi.load(prefix + "/" + prefix + '_ksi' + suffix)
#        relief = d.ScaledRelief.load(prefix + '/' + prefix + '_relief' + suffix)
	d8.sort()
	d8._sort_indexes = filled.sort(reverse = False, force = True)
	elevation = d.Elevation.load(full_prefix + '_elevation')
	ksi = d.Ksi(area = area, flow_direction = d8, theta = theta, Ao = Ao, flow_length = length)
	relief = d.ScaledRelief(flow_direction = d8, elevation = elevation, flow_length = length, Ao = Ao, theta = theta, area = area)
	relief.save(prefix + "/" + prefix + '_' + str(A0) + '_' + str(theta) + '_relief')
	ksi.save(prefix + "/" + prefix + '_' + str(Ao) + '_' + str(theta) + '_ksi')
        theta = float(suffix.replace('_','.'))
        print(theta)
        local_dict[suffix] = dra.calculate_ks_for_sample(locs, d8, ksi, relief, area, Ao, theta)
    
    	ks[prefix] = local_dict

pickle.dump( ks, open( "ks.p", "wb" ) )
