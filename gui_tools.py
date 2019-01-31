import matplotlib
matplotlib.use('TKAgg')

import dem as d
import matplotlib.pylab as plt
import numpy as np
import pickle as p
import matplotlib.cm as cm
import sys

if sys.version_info[0] >= 3:
    raw_input = input
    
def plot_dem(dem, hs):
    plt.close('all')

    dem.plot()
    hs.plot(alpha = 0.5, cmap = plt.cm.gray)

def select_outlets(dem, fd, prefix, hs = None, outlet_filename = 'outlets.p', color = 'r'):
    if hs is None:
        hs = d.Hillshade(elevation = dem, azimuth = 320, inclination = 20)
        
    plot_dem(dem, hs)
    
    keep_going = True
    
    while keep_going:
        zoom_ok = False
        print('\nZoom or pan to view, \npress spacebar when ready to select point upstream of outlet:\n')
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()
        print('\nClick on point upstream of outlet.')
        xy = plt.ginput(1)[0]
        xy_path = fd.search_down_flow_direction_from_xy_location(xy)
        plt.plot(xy)
        plt.plot(xy_path)
        xy_path_plot = list(zip(*xy_path))
        path = plt.plot(xy_path_plot[0],xy_path_plot[1], color+'-')
        print('\nClick on the outlet location.')
        xy = plt.ginput(1)[0]
    
        min_distance = 1E12
    
        for loc in xy_path:
        
            if (np.sqrt( np.power(loc[0] - xy[0], 2) + np.power(loc[1] - xy[1], 2)) < min_distance) or min_distance == 1E12:
                outlet_loc = loc
                min_distance = np.sqrt( np.power(loc[0] - xy[0], 2) + np.power(loc[1] - xy[1], 2))
    
        plt.figure(1)
        plt.plot(outlet_loc[0], outlet_loc[1], color+'o')
        path.pop(0).remove()
    
        outlet_prefix = raw_input('Type a name for this outlet (leaving this blank will prevent outlet from being saved and will complete the selection process: ')
    
        if outlet_prefix != '':
            import os.path
            if os.path.isfile(outlet_filename):
                outlets = p.load(open(outlet_filename, 'rb'))
            else:
                outlets = dict()
        
            this_tile = outlets.get(prefix, dict())
            this_tile[outlet_prefix] = outlet_loc
    
            outlets[prefix] = this_tile
        
            p.dump(outlets, open(outlet_filename, 'wb'))
        else:
            keep_going = False
    
