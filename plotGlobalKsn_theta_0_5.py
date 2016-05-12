prefixes = ('af', 'as', 'au', 'ca', 'eu', 'na', 'sa')

import dem as d
from matplotlib import pyplot as plt
from demMethods import plotGrids

for prefix in prefixes:
    
    ksi = d.Ksi.load(prefix + '_ksi_250000_0_0_5')
    relief = d.ScaledRelief.load(prefix + '_relief_250000_0_0_5')
    
    fig = plt.figure()
    plotGrids(ksi, relief, 'k.', rasterized = True)
    