import dem as d
import csv
import numpy as np


def read_csv(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        return list(reader)
    
def best_ksn(ksi, scaled_relief, xo = 500):
    
    A = np.vstack([ksi-np.ones(len(ksi))*float(xo)]).T 
    return np.linalg.lstsq(A, scaled_relief)

def best_chi_ksn(chi, scaled_relief, xo = 500):
    
    A = np.vstack([chi]).T
    return np.linalg.lstsq(A, scaled_relief)

def best_ks_and_theta_with_wrss(elevation, flow_direction, area, outlet, xo = 500):
    
    import scipy.optimize
    chi_ks = lambda theta: best_ks_with_wrss(elevation, flow_direction, area, theta, outlet, xo)[1]
    if len(chi_ks(0.5)) == 0:
        return (0, 0, 0)
    xopt = scipy.optimize.fmin(func=chi_ks, x0=np.array([0.5]))
    (m, WRSS) = best_ks_with_wrss(elevation, flow_direction, area, xopt[0], outlet, xo)
    SS = uninformative_SS(elevation, flow_direction, area, outlet, xopt[0], xo)
    
    return (m, xopt[0], 1 - (WRSS / SS))
    
def best_ks_and_theta_with_clip(elevation, flow_direction, area, outlet, area_measured, xo = 500, pixel_radius = 5):
    
    bounds = flow_direction.bounds_of_basin_for_outlet(outlet);
    elev_clip = elevation.clip_to_bounds(bounds)
    area_clip = area.clip_to_bounds(bounds)
    fd_clip = flow_direction.clip_to_bounds(bounds)
    outlet_xy = (outlet[1], outlet[0])
    locations_snap = area_clip.snap_locations_to_closest_value((outlet_xy, ), (area_measured, ), pixel_radius = pixel_radius)
    return best_ks_and_theta_with_wrss(elev_clip, fd_clip, area_clip, locations_snap, xo)

def best_ks_with_wrss(elevation, flow_direction, area, theta, outlet, xo = 500):
    
    chi = d.GeographicChi(area = area, flow_direction = flow_direction, theta = theta, outlets = (outlet[0], ), Ao = np.power(xo,2))
    valid_indexes = np.where((chi._griddata != np.NAN) & (elevation._griddata != np.NAN) & (chi._griddata != 0.0))
    outlet_rowscols = flow_direction._xy_to_rowscols((outlet[0], ))[0]
    base_elevation = elevation[outlet_rowscols[0], outlet_rowscols[1]]
    sol = best_chi_ksn(chi._griddata[valid_indexes], elevation._griddata[valid_indexes]-base_elevation)
    m = sol[0]
    WRSS = sol[1]    
    
    return (m, WRSS)

def uninformative_SS(elevation, flow_direction, area, outlet, theta, xo = 500):
    
    chi = d.GeographicChi(area = area, flow_direction = flow_direction, theta = theta, outlets = (outlet[0], ), Ao = np.power(xo,2))
    valid_indexes = np.where((chi._griddata != np.NAN) & (elevation._griddata != np.NAN) & (chi._griddata != 0.0))
    mean_elevation = np.mean(elevation._griddata[valid_indexes])
    return np.sum(np.power(elevation._griddata[valid_indexes]-mean_elevation, 2))
    
def find_ksi_scaled_relief(lat, lon, area, ksi, relief, d8, A_measured, pixel_radius = 5):
    
    index = area._xy_to_rowscols(((lon,lat),))[0]
    if index[0] is None:
        return None, None, None
    #row, col = area.find_nearest_cell_with_value(index, A_measured, pixel_radius)
    A_calculated = area[index[0], index[1]]
    indexes_of_area = d8.get_indexes_of_upstream_cells(index[0], index[1])
    ksi_values = list()
    relief_values = list()
    for (row, col) in indexes_of_area:
        if ksi[row,col] is None or relief[row,col] is None or np.isnan(ksi[row, col]) or np.isnan(relief[row,col]):
            return None, None, None
        ksi_values.append(ksi[row,col])
        relief_values.append(relief[row,col])
    return ksi_values, relief_values, A_calculated
    
def calculate_ksn_for_data(data, Ao = 250000, theta = 0.5):
   
    import sys
    sys.setrecursionlimit(1000000)
    
    pixel_radius = 5
    reject_fraction = 0.1
    
    prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']
    lats = list()
    lons = list()
    areas = list()
    for sample_name, lat, lon, dr, dr_sig, a in data:
        lats.append(float(lat))
        lons.append(float(lon))
        areas.append(float(a)*1.0E6)
    
    locations = zip(lons,lats)
    ksn_vec = np.zeros(len(areas), dtype = np.float64)
    a_calc_vec = np.zeros(len(areas), dtype = np.float64)
    
    for prefix in prefixes:
        print('Loading prefix: ' + prefix)
        area = d.GeographicArea.load(prefix + '_area')
        #ksi = d.GeographicKsi.load(prefix + '_ksi_' + suffix)
        #relief = d.ScaledRelief.load(prefix + '_relief_' + suffix)
        d8 = d.FlowDirectionD8.load(prefix + '_flow_direction')
        elevation = d.Elevation.load(prefix + '_elevation')
        print('Done loading prefix: ' + prefix)
        counter = 0
        xo = np.mean(d8._mean_pixel_dimension(flow_direction = d8) * d8.pixel_scale())
        locs = tuple()
        areas_for_valid_points = tuple()
        for (lon, lat), area_m in zip(locations, areas):
            if elevation.location_in_grid((lon, lat)) and area_m > 5.0*Ao:
                locs = locs + ((counter, (lon, lat)),)
                areas_for_valid_points = areas_for_valid_points + (area_m, )
            counter = counter + 1
        counter = tuple([el[0] for el in locs]);
        locat = tuple([el[1] for el in locs]);
        locations_snap = area.snap_locations_to_closest_value(locat, areas_for_valid_points, pixel_radius = pixel_radius)
        locations_snap_indexes = area._xy_to_rowscols(locations_snap)
        dem_derived_areas = [area[i,j] for (i,j) in locations_snap_indexes]
        for (lon, lat), target_area, measured_area in zip(locat, areas_for_valid_points, dem_derived_areas):
            print('Longitude: ' + str(lon) + '; Latitude: ' + str(lat) + '; target area: ' + str(target_area) + '; measured area: ' + str(measured_area) + '; fractional difference: ' + str(np.abs(measured_area - target_area)/target_area))

        fraction_difference = [np.abs(derived_area - measured_area) / measured_area for (derived_area, measured_area) in zip(dem_derived_areas, areas_for_valid_points)]
        chi = d.GeographicChi(area = area, flow_direction = d8, theta = theta, Ao = Ao, outlets = (locations_snap[0],))
        scaled_relief = d.ChiScaledRelief(elevation = elevation, flow_direction = d8, theta = theta, Ao = Ao, outlets = (locations_snap[0],)) 
        first_element = True
            
        for (lon, lat), areas_m, counter_v, areas_dem, sample_fraction_difference in zip(locations_snap, areas_for_valid_points, counter, dem_derived_areas, fraction_difference):
            if first_element:
                first_element = False
            else:
                chi._create_from_inputs(area = area, flow_direction = d8, theta = theta, Ao = Ao, outlets = ((lon, lat),))
                scaled_relief._create_from_inputs(elevation = elevation, flow_direction = d8, theta = theta, Ao = Ao, outlets = ((lon, lat),)) 
            chi_vec, scaled_relief_vec, a_calc = find_ksi_scaled_relief(lat, lon, area, chi, scaled_relief, d8, area_m, pixel_radius)
            if chi_vec is not None and sample_fraction_difference < reject_fraction:
                best_fit, residuals, rank, s = best_ksn(chi_vec, scaled_relief_vec)
                best_ks = best_fit[0]
                ksn_vec[counter_v] = best_ks
                a_calc_vec[counter_v] = a_calc
                print('lat = {0}, long = {1}, ksn = {2}, np = {3}, reported area = {4}, measured area = {5}, np/4 = {6}'.format(lat,lon,best_ks, len(chi_vec), areas_m/1.0E6, areas_dem/1.0E6, len(chi_vec)/4))
        
        area = None
        d8 = None
        elevation = None
        
    return ksn_vec, a_calc_vec


def extract_all_ksi_relief_values_for_position(position, d8, area, ksi, relief, Ao, mask=None, A_cutoff = None):
    (row, col) = area._xy_to_rowscols((position, ))[0]
        
    indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
    ksi_values = list()
    relief_values = list()
    if A_cutoff is None:
        A_cutoff = Ao
    for (row, col) in indexes_of_area:
        if ksi[row,col] is None or relief[row,col] is None or np.isnan(ksi[row, col]) or np.isnan(relief[row,col]):
            return None, None, None
        if area[row, col] >= A_cutoff:
            if mask is None:
                ksi_values.append(ksi[row,col])
                relief_values.append(relief[row,col])
            elif mask[row, col] == 1:
                ksi_values.append(ksi[row,col])
                relief_values.append(relief[row,col])
    return ksi_values, relief_values

def create_chi_grid_for_geographic_prefix(prefix, thetas, Ao, basin_lengths):
    
    area = d.GeographicArea.load(prefix + "_area")
    d8 = d.FlowDirectionD8.load(prefix + "_flow_direction")
    flow_length = d.GeographicFlowLength.load(prefix + "_flow_length")
    elevation = d.Elevation.load(prefix + "_elevation")
    for theta in thetas:
        for basin_length in basin_lengths:
            chi = d.GeographicChi(area = area, flow_direction=d8, flow_length=flow_length, theta=theta, Ao=Ao, basin_length=basin_length)
            chi.save(prefix + '_chi_' + str(basin_length) + '_' + str(theta).replace('.', '_') + '_' + str(Ao))
            chi = None
            relief = d.ChiScaledRelief(elevation = elevation, area = area, flow_direction=d8, flow_length=flow_length, theta=theta, Ao=Ao, basin_length=basin_length)
            relief.save(prefix + '_relief_' + str(basin_length) + '_' + str(theta).replace('.', '_') + '_' + str(Ao))
            relief = None
                                 
def calculate_ks_for_sample(v, d8, ksi, relief, area, Ao = 250000, mask = None, xo = None, A_cutoff = None):
        
    ks = list()
    if xo is None:
        xo = np.mean(d8._mean_pixel_dimension(flow_direction = d8) * d8.pixel_scale())
    for position in v:
        ksi_values, relief_values = extract_all_ksi_relief_values_for_position(position, d8, area, ksi, relief, Ao, mask, A_cutoff=A_cutoff)
        from matplotlib import pyplot as plt
        try:
            best_fit, residuals, rank, s = best_chi_ksn(ksi_values, relief_values, xo)
            best_ks = best_fit[0] 
            model_residuals = residuals[0] 
            relief_array = np.array(relief_values)
            relief_mean = np.mean(relief_array)
            total_residuals = np.sum((relief_array - relief_mean)**2)
            R2 = 1 - model_residuals / total_residuals        
        except:
            best_ks = 0
            R2 = 0    
        ks.append((best_ks, R2))

    return ks

def plot_relief_and_ksi(v, d8, ksi, relief, area, Ao = 250000, mask = None):
    from matplotlib import pyplot as plt
    for position in v:
	ksi_values, relief_values = extract_all_ksi_relief_values_for_position(position, d8, area, ksi, relief, Ao, mask = mask)
        plt.figure()
        plt.plot(ksi_values, relief_values, 'k.', rasterized = True)


def calculate_slope_fraction_for_sample(v, d8, area, slope, cutoff = 0.2):
        
    fraction = list()
    
    for position in v:
        
        (row, col) = area._xy_to_rowscols((position, ))[0]
        
        indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
        total_number_of_points_in_basin = len(indexes_of_area)
        number_of_points_in_basin_above_cutoff = 0
        for (row, col) in indexes_of_area:
            print(slope[row,col])
            if slope[row,col] > cutoff:
                number_of_points_in_basin_above_cutoff += 1
                    
        fraction.append(number_of_points_in_basin_above_cutoff / total_number_of_points_in_basin)

    return fraction

def plot_stock_and_montgomery():
    
    K = {'granitoids': (4.4e-7, 4.3e-6),
         'volcaniclastics': (4.8e-5, 3.0e-4),
         'mudstones': (4.7e-4, 7.0e-3),
         'basalt': (3.8e-6, 7.3e-6)}
    
    colors = {'granitoids': 'r-',
         'volcaniclastics': 'b-',
         'mudstones': 'g-',
         'basalt': 'm-'}
    
    from matplotlib import pyplot as plt
    
    for key in K.keys():
        
        U = (1e-4, 1e1)
        ks = (U[0] / 1000.0 / K[key][0], U[1] / 1000.0 / K[key][0])
        plt.plot(U,ks,colors[key])
        ks = (U[0] / 1000.0 / K[key][1], U[1] / 1000.0 / K[key][1])
        plt.plot(U,ks,colors[key])

        

    
