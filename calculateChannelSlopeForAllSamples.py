import denudationRateAnalysis as dra
import dem as d
import numpy as np
import pickle

prefixes = ['brazil', 'costarica', 'guatemala', 'taiwan', 'venezuela']

slp = dict()

for prefix in prefixes:
    print(prefix)

    d8 = d.FlowDirectionD8.load(prefix + "/" + prefix + '_flow_direction')
    elevation = d.Elevation.load(prefix + "/" + prefix + '_elevation')
    area = d.Area.load(prefix + '/' + prefix + '_area')
    slope = d.ChannelSlope(elevation = elevation, flow_direction = d8)
    slope.save(prefix + "/" + prefix + "_channelslope")
    
    locs = np.load(prefix + '/' + prefix +  '_sample_locations.npy')
    
    slp[prefix] = dra.calculate_slope_fraction_for_sample(locs, d8, area, slope)
    
pickle.dump( slp, open( "slp.p", "wb" ) )