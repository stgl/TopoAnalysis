import dem as d
import csv
import numpy as np

def read_csv(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        return list(reader)
    
def best_ksn(ksi, scaled_relief):
    
    A = np.vstack(ksi).T
    return np.linalg.lstsq(A, scaled_relief)[0]

def find_ksi_scaled_relief(lat, lon, area, ksi, relief, d8, A_measured, pixel_radius = 5):
    
    index = area._xy_to_rowscols(((lon,lat),))[0]
    row, col = area.find_nearest_cell_with_value(index, A_measured, pixel_radius)
    A_calculated = area[row,col]
    indexes_of_area = d8.get_indexes_of_upstream_cells(row, col)
    ksi_values = list()
    relief_values = list()
    for (row, col) in indexes_of_area:
        ksi_values.append(ksi[row,col])
        relief_values.append(relief[row,col])
    return ksi_values, relief_values, A_calculated
    
