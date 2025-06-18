from . import dem as d
import matplotlib.pylab as plt
import numpy as np

def plot_downstream_profile(elevation, flow_direction, outlet, plot_code, downstream = True, start_at = 0.0, mean_pixel_dimension = None, figure = None):
    
    rc = flow_direction.search_down_flow_direction(outlet)
    
    rows = np.array(list(zip(*rc))[0], dtype = np.int)
    cols = np.array(list(zip(*rc))[1], dtype = np.int)
    
    if downstream:
        direction = 1.0
    else:
        direction = -1.0
        
    if mean_pixel_dimension is None:
        pixel_dimension = np.ones_like(rows) * flow_direction._georef_info.dx
    else:
        pixel_dimension = mean_pixel_dimension[rows, cols]
        
    for i in range(rows.shape[0]):
        if i == 0:
            length = np.array([start_at])
        else:
            if (rows[i] != rows[i-1]) & (cols[i] != cols[i-1]):
                length = np.append(length, np.array([length[i-1] + pixel_dimension[i-1] * 1.414 * direction]))
            else:
                length = np.append(length, np.array([length[i-1] + pixel_dimension[i-1] * direction]))
    
    elevation_profile = elevation._griddata[rows, cols]
    
    if figure is None:
        figure = plt.figure()
            
    plt.figure(figure.number)
    plt.plot(length, elevation_profile, plot_code)
    return rows, cols, length, elevation_profile
    
def plot_recursive_upstream_profiles(elevation, flow_direction, area, outlet, plot_code, downstream = False, start_at = 0.0, figure = None, minimum_area = 1.0E6):
    
    def plot_ld_link(current_length, ld_list, plot_code, downstream_sign, minimum_area):
        if len(ld_list['index']) == 1:
            (current_row, current_column) = ld_list['index'][0]
        else:
            (current_row, current_column) = ld_list['index']

        
        if ld_list.get('next') is None:
            return
        for next_list in ld_list['next']:
            
            if next_list['area'] >= minimum_area:
                
                (next_row, next_column) = next_list['index']
                if (current_row != next_row) & (current_column != next_column):
                    next_length = current_length + (ld_list['de'] * 1.414 * downstream_sign)
                else:
                    next_length = current_length + (ld_list['de'] * downstream_sign)
                plt.plot([current_length, next_length], [ld_list['elevation'], next_list['elevation']], plot_code)
                plot_ld_link(next_length, next_list, plot_code, downstream_sign, minimum_area)
                
    mean_pixel_dimension = d.BaseSpatialGrid()
    mean_pixel_dimension._copy_info_from_grid(area, True)
    mean_pixel_dimension._griddata = area._mean_pixel_dimension()
                
    ld_list = flow_direction.map_values_to_recursive_list(outlet, elevation = elevation, area = area, de = mean_pixel_dimension)
    current_length = start_at
    if downstream:
        downstream_sign = 1.0
    else:
        downstream_sign = -1.0
        
    #if figure is None:
    #    figure = plt.figure()
    
    #plt.figure(figure.number)
    plot_ld_link(current_length, ld_list, plot_code, downstream_sign, minimum_area)
    

                        
def plot_profiles(elevation, flow_direction, area, outlet, plot_code, minimum_area = 1.0E6, figure = None):
    
    if figure is None:
        figure = plt.figure()
        
    plot_downstream_profile(elevation, flow_direction, outlet, plot_code, mean_pixel_dimension= area._mean_pixel_dimension(), figure = figure)        
    plot_recursive_upstream_profiles(elevation, flow_direction, area, outlet, plot_code, figure = figure, minimum_area = minimum_area)
    
def plot_profiles_with_outlet_code(prefix, code, plot_code, dem, fd, area, minimum_area = 1.0E7):
    
    import pickle
    outlets = pickle.load(open('outlets.p', 'rb'))
    outlet = outlets[prefix][code]
    plot_profiles(dem, fd, area, outlet, plot_code, minimum_area)
    
def plot_chi_profiles(elevation, flow_direction, area, outlet, plot_code, minimum_area = 1.0E6, figure = None, theta = 0.5, start_at = 0.0, downstream = True, Ao = 1.0E6):
    
    ((row, col),) = elevation._xy_to_rowscols((outlet,))
    base_elevation = elevation[row,col]
    
    def plot_ld_link(current_chi, ld_list, plot_code, downstream_sign, minimum_area):
        (current_row, current_column) = ld_list['index']
        
        if ld_list.get('next') is None:
            return
        for next_list in ld_list['next']:
            
            if next_list['area'] >= minimum_area:
                
                (next_row, next_column) = next_list['index']
                if (current_row != next_row) & (current_column != next_column):
                    next_chi = current_chi + (1 / np.array(next_list['area']))**theta * np.array(ld_list['de']*1.414*downstream_sign)
                                        
                else:
                    next_chi = current_chi + (1 / np.array(next_list['area']))**theta * np.array(ld_list['de']*downstream_sign)
                plt.plot([current_chi, next_chi], [(ld_list['elevation'] - base_elevation) * np.power(Ao,theta), (next_list['elevation'] - base_elevation)*np.power(Ao,theta)], plot_code)
                plot_ld_link(next_chi, next_list, plot_code, downstream_sign, minimum_area)
                
    mean_pixel_dimension = d.BaseSpatialGrid()
    mean_pixel_dimension._copy_info_from_grid(area, True)
    mean_pixel_dimension._griddata = area._mean_pixel_dimension()
                
    ld_list = flow_direction.map_values_to_recursive_list(outlet, elevation = elevation, area = area, de = mean_pixel_dimension)
    current_chi = start_at
    if downstream:
        downstream_sign = 1.0
    else:
        downstream_sign = -1.0
        
    if figure is None:
        figure = plt.figure()
    
    plt.figure(figure.number)
    plot_ld_link(current_chi, ld_list, plot_code, downstream_sign, minimum_area)  

def plot_chi_profiles_with_outlet_code(prefix, code, plot_code, dem, fd, area, minimum_area = 1.0E7, theta = 0.5):
    
    import pickle
    outlets = pickle.load(open('outlets.p', 'rb'))
    outlet = outlets[prefix][code]
    plot_chi_profiles(dem, fd, area, outlet, plot_code, minimum_area = minimum_area, theta = theta)


def interactive_chi_profiles_and_map_view(prefix, code, plot_code, dem, fd, area, hillshade, minimum_area = 1.0E7, theta = 0.5, Ao = 1.0E6):
        
    
    import pickle
    outlets = pickle.load(open('outlets.p', 'rb'))
    outlet = outlets[prefix][code]
    
    from demRecursionTools import map_chi_profiles
    
    chi_map = map_chi_profiles(dem, fd, area, outlet, minimum_area = minimum_area, theta = theta, Ao = Ao)
    indexes = chi_map.keys()
    import operator
    (chi, _) = list(zip(*operator.itemgetter(*indexes)(chi_map)))
    coordinates = list(zip(*hillshade._rowscols_to_xy(indexes)))
    
    fig1 = plt.figure(1)
    
    hillshade.plot(cmap = plt.cm.gray)
    plt.scatter(coordinates[0], coordinates[1], c=chi, s=0.5)
    ax = plt.gca()
    
    fig2 = plt.figure(2)
    plot_chi_profiles(dem, fd, area, outlet, plot_code, minimum_area = minimum_area, theta = theta, figure=fig2, Ao = Ao)
    
    current_marker = plt.plot(0,0,'bo');
    
    def hover(current_marker, Ao, theta):
        def hoverwrapper(event):
            if event.inaxes == ax:
                x = event.xdata
                y = event.ydata
                (i, ) = hillshade._xy_to_rowscols(((x, y), ))
                chi_elevation = chi_map.get(i)
                if chi_elevation is not None:
                    (chi, elevation) = chi_elevation
                    current_marker[0].set_data([chi], [elevation*np.power(Ao, theta)])
            
        return hoverwrapper
    
    fig1.canvas.mpl_connect('motion_notify_event', hover(current_marker, Ao, theta))
