import denudationRateAnalysis as dra
import dem as d
import numpy as np
import pickle

prefixes = ['brazil', 'costarica', 'guatemala', 'taiwan', 'venezuela']
suffixes = ['0_4', '0_5', '0_6']
theta_values = [0.4, 0.5, 0.6]

area_dic = dict()

for prefix in prefixes:
    print(prefix)
    area = d.Area.load(prefix + '/' + prefix + '_area')
    locs = np.load(prefix + '/' + prefix +  '_sample_locations.npy')
    locs_rowscols = area._xy_to_rowscols(locs)
    area_vec = list()
    for loc in locs_rowscols:
        area_vec.append(area[loc[0],loc[1]])
    area_dic[prefix] = tuple(area_vec)
    
pickle.dump( area_dic, open( "basin_areas.p", "wb" ) )