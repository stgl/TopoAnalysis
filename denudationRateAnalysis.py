import dem as d
import csv
import numpy as np

def read_csv(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        return list(reader)
    
def best_ksn(ksi, scaled_relief, xo = 500, theta = 0.5):
    
    A = np.vstack(ksi) - xo 
    return np.linalg.lstsq(A, scaled_relief)[0]

def find_ksi_scaled_relief(lat, lon, area, ksi, relief, d8, A_measured, pixel_radius = 5):
    
    index = area._xy_to_rowscols(((lon,lat),))[0]
    if index[0] is None:
        return None, None, None
    row, col = area.find_nearest_cell_with_value(index, A_measured, pixel_radius)
    A_calculated = area[row,col]
    indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
    ksi_values = list()
    relief_values = list()
    for (row, col) in indexes_of_area:
        if ksi[row,col] is None or relief[row,col] is None or np.isnan(ksi[row, col]) or np.isnan(relief[row,col]):
            return None, None, None
        ksi_values.append(ksi[row,col])
        relief_values.append(relief[row,col])
    return ksi_values, relief_values, A_calculated
    
def calculate_ksn_for_data(data, Ao = 250000, theta = 0.5, xo = 500):
    
    import sys
    sys.setrecursionlimit(1000000)
    
    prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']

    suffix = str(Ao) + '_' + str(theta).replace('.', '_')
    
    lats = list()
    lons = list()
    areas = list()
    for sample_name, lat, lon, dr, dr_sig, a in data:
        lats.append(float(lat))
        lons.append(float(lon))
        areas.append(float(a))
    
    locations = zip(lons,lats)
    ksn_vec = np.zeros(len(areas), dtype = np.float64)
    a_calc_vec = np.zeros(len(areas), dtype = np.float64)
    
    for prefix in prefixes:
        print('Loading prefix: ' + prefix)
        area = d.GeographicArea.load(prefix + '_area')
        ksi = d.GeographicKsi.load(prefix + '_ksi_' + suffix)
        relief = d.ScaledRelief.load(prefix + '_relief_' + suffix)
        d8 = d.FlowDirectionD8.load(prefix + '_flow_direction')

        print('Done loading prefix: ' + prefix)
        counter = 0
        
        for (lon, lat), area_m in zip(locations, areas):
            
            ksi_vec, relief_vec, a_calc = find_ksi_scaled_relief(lat, lon, area, ksi, relief, d8, area_m*1.0e6, 15)
            if ksi_vec is not None and (abs(area_m*1.0e6 - a_calc) < abs(area_m*1.0e6 - a_calc_vec[counter])):
                ksn = best_ksn(ksi_vec, relief_vec, xo, theta)[0]
                ksn_vec[counter] = ksn
                a_calc_vec[counter] = a_calc
                print 'lat = {0}, long = {1}, ksn = {2}'.format(lat,lon,ksn)

            counter = counter + 1
        
    return ksn_vec, a_calc_vec
                
def calculate_ks_for_sample(v, d8, ksi, relief, area, Ao = 250000, theta = 0.5):
        
    ks = list()
    
    for position in v:
        
        (row, col) = area._xy_to_rowscols((position, ))[0]
        
        indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
        ksi_values = list()
        relief_values = list()
        for (row, col) in indexes_of_area:
            if ksi[row,col] is None or relief[row,col] is None or np.isnan(ksi[row, col]) or np.isnan(relief[row,col]):
                return None, None, None
            if area[row, col] >= Ao:
                ksi_values.append(ksi[row,col])
                relief_values.append(relief[row,col])
        
        best_ks = best_ksn(ksi_values, relief_values, 90.0, theta)[0]
        
        ksi_array = np.array(ksi_values)
        relief_array = np.array(relief_values)
        relief_mean = np.mean(relief_array)
        total_residuals = np.sum((relief_array - relief_mean)**2)
        model_residuals = np.sum((best_ks * (ksi_array - 90.0) - relief_array )**2)
        R2 = 1 - model_residuals / total_residuals            
        ks.append((best_ks, R2))

    return ks

def calculate_slope_fraction_for_sample(v, d8, area, slope, cutoff = 0.2):
        
    fraction = list()
    
    for position in v:
        
        (row, col) = area._xy_to_rowscols((position, ))[0]
        
        indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
        total_number_of_points_in_basin = len(indexes_of_area)
        number_of_points_in_basin_above_cutoff = 0
        for (row, col) in indexes_of_area:
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

        

    
