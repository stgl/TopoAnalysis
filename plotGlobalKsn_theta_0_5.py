prefixes = ( 'au', 'ca', 'na', 'sa')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_0_0_5')
    relief = d.ScaledRelief.load(prefix + '_relief_250000_0_0_5')
    
    fig = plt.figure()
    ksi_vec, relief_vec = plotGrids(ksi, relief, 'k.', rasterized = True)
    plt.savefig(prefix + '0_5.png',600)
    plt.close()
    
    if 'all_ksi_vec' not in locals():
        all_ksi_vec = ksi_vec
        all_relief_vec = relief_vec
    else:
        all_ksi_vec.append(ksi_vec)
        all_relief_vec.append(relief_vec)
        
    
fig = plt.figure()

plt.plot(all_ksi_vec, all_relief_vec, 'k.', rasterize = True)
plt.savefig('all_0_5.png',600)
    